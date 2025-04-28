#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <curl/curl.h>
#include"rapidjson/document.h"

using namespace std;
using namespace rapidjson;

// ——— Simple BlockingQueue<T> ———
template<typename T>
class BlockingQueue {
    queue<T>       q;
    mutex          m;
    condition_variable cv;
    bool           done{false};
public:
    void push(const T& v) {
        { lock_guard<mutex> lk(m); q.push(v); }
        cv.notify_one();
    }
    // pop(): returns false if queue empty & done()
    bool pop(T& out) {
        unique_lock<mutex> lk(m);
        cv.wait(lk, [&]{ return done || !q.empty(); });
        if (q.empty()) return false;
        out = q.front(); q.pop();
        return true;
    }
    void finish() {
        { lock_guard<mutex> lk(m); done = true; }
        cv.notify_all();
    }
};
// ————————————————————————

static const string BASE_URL = 
    "http://hollywood-graph-crawler.bridgesuncc.org/neighbors/";
static const int    MAX_THREADS = 4;

// libcurl write callback
static size_t WriteCb(void* data, size_t sz, size_t nm, string* out) {
    out->append((char*)data, sz * nm);
    return sz * nm;
}

// Fetch & parse neighbors in one go
vector<string> get_neighbors(const string& node) {
    vector<string> nbrs;
    CURL* curl = curl_easy_init();
    if (!curl) return nbrs;

    // build URL with escaping
    char* esc = curl_easy_escape(curl, node.c_str(), (int)node.size());
    string url = BASE_URL + esc;
    curl_free(esc);

    // perform HTTP GET
    string resp;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCb);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &resp);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "ParallelCrawler/1.0");
    curl_easy_perform(curl);
    curl_easy_cleanup(curl);

    // parse JSON
    Document doc;
    doc.Parse(resp.c_str());
    if (doc.HasMember("neighbors") && doc["neighbors"].IsArray()) {
        for (auto& v : doc["neighbors"].GetArray())
            if (v.IsString())
                nbrs.emplace_back(v.GetString());
    }
    return nbrs;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <start_node> <depth>\n";
        return 1;
    }
    string start = argv[1];
    int    maxD  = stoi(argv[2]);

    // ensure output folder exists
    filesystem::create_directories("output");

    // shared data
    BlockingQueue<pair<string,int>> queue;
    unordered_set<string>           visited;
    vector<string>                  output;
    mutex                           lock;
    atomic<int>                     inFlight{1};

    // seed BFS
    visited.insert(start);
    queue.push({start, 0});

    // worker threads
    auto worker = [&]() {
        pair<string,int> cur;
        while (queue.pop(cur)) {
            auto [node, depth] = cur;

            {   // record node
                lock_guard<mutex> lk(lock);
                output.push_back(node);
            }

            if (depth < maxD) {
                for (auto& nbr : get_neighbors(node)) {
                    lock_guard<mutex> lk(lock);
                    if (visited.insert(nbr).second) {
                        ++inFlight;
                        queue.push({nbr, depth + 1});
                    }
                }
            }

            if (--inFlight == 0)
                queue.finish();
        }
    };

    // launch & time
    auto t0 = chrono::steady_clock::now();
    vector<thread> thr;
    for (int i = 0; i < MAX_THREADS; ++i)
        thr.emplace_back(worker);
    for (auto& t : thr) t.join();
    auto t1 = chrono::steady_clock::now();
    double elapsed = chrono::duration<double>(t1 - t0).count();

    // console output
    for (auto& n : output)
        cout << "-> " << n << "\n";
    cout << "Timer: " << elapsed << " Seconds\n";

    // file dump
    ofstream ofs("output/" + start + "_depth" + to_string(maxD) + ".txt");
    for (auto& n : output) ofs << n << "\n";
    ofs << "It took: " << elapsed << " Seconds to run\n" ;

    return 0;
}
