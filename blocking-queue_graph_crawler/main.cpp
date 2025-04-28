#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <unordered_set>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <curl/curl.h>
#include <stdexcept>
#include "rapidjson/error/error.h"
#include "rapidjson/reader.h"
#include <rapidjson/document.h>

using namespace std;
using namespace rapidjson;

// ---- your existing JSON / curl helpers ----

struct ParseException : runtime_error, ParseResult {
    ParseException(ParseErrorCode c, const char* m, size_t o)
      : runtime_error(m), ParseResult(c, o) {}
};
#define RAPIDJSON_PARSE_ERROR_NORETURN(code, offset) \
    throw ParseException(code, #code, offset)

bool debug = false;
const string SERVICE_URL = "http://hollywood-graph-crawler.bridgesuncc.org/neighbors/";

string url_encode(CURL* curl, const string& input) {
    char* out = curl_easy_escape(curl, input.c_str(), (int)input.size());
    string s = out; curl_free(out);
    return s;
}

size_t WriteCallback(void* contents, size_t size, size_t nmemb, string* output) {
    size_t total = size * nmemb;
    output->append((char*)contents, total);
    return total;
}

string fetch_neighbors(CURL* curl, const string& node) {
    string url = SERVICE_URL + url_encode(curl, node);
    string response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    struct curl_slist* headers = nullptr;
    headers = curl_slist_append(headers, "User-Agent: ParallelCrawler/1.0");
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);

    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK)
        cerr << "CURL error: " << curl_easy_strerror(res) << endl;
    curl_slist_free_all(headers);
    return (res == CURLE_OK) ? response : "{}";
}

vector<string> get_neighbors(const string& json_str) {
    vector<string> neighbors;
    Document doc;
    doc.Parse(json_str.c_str());
    if (doc.HasParseError())
        RAPIDJSON_PARSE_ERROR_NORETURN(doc.GetParseError(), doc.GetErrorOffset());
    if (doc.HasMember("neighbors") && doc["neighbors"].IsArray()) {
        for (auto& v : doc["neighbors"].GetArray())
            neighbors.push_back(v.GetString());
    }
    return neighbors;
}

// ---- blocking queue & BFS ----

template<typename T>
class BlockingQueue {
public:
    void push(const T& x) {
        { lock_guard<mutex> lk(m); q.push(x); }
        cv.notify_one();
    }
    bool pop(T& out) {
        unique_lock<mutex> lk(m);
        cv.wait(lk, [&]{ return closed || !q.empty(); });
        if (q.empty()) return false;
        out = q.front(); q.pop();
        return true;
    }
    void close() {
        { lock_guard<mutex> lk(m); closed = true; }
        cv.notify_all();
    }
private:
    queue<T> q; mutex m; condition_variable cv; bool closed = false;
};

struct Node { string name; int depth; };

const int MAX_THREADS = 4;

vector<string> parallel_bfs(const string& start, int maxDepth) {
    BlockingQueue<Node> workQ;
    unordered_set<string> visited;
    mutex visitM, resultM;
    vector<string> result;
    atomic<int> tasks{1};

    visited.insert(start);
    workQ.push({start, 0});

    auto worker = [&]{
        CURL* curl = curl_easy_init();
        Node cur;
        while (workQ.pop(cur)) {
            // record
            {
                lock_guard<mutex> lk(resultM);
                result.push_back(cur.name);
            }
            // expand
            if (cur.depth < maxDepth) {
                auto js = fetch_neighbors(curl, cur.name);
                for (auto& n : get_neighbors(js)) {
                    lock_guard<mutex> lk(visitM);
                    if (!visited.count(n)) {
                        visited.insert(n);
                        tasks++;
                        workQ.push({n, cur.depth + 1});
                    }
                }
            }
            if (--tasks == 0) 
                workQ.close();
        }
        curl_easy_cleanup(curl);
    };

    vector<thread> threads;
    for (int i = 0; i < MAX_THREADS; ++i)
        threads.emplace_back(worker);
    for (auto& t : threads) 
        t.join();

    return result;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <start_node> <depth>\n";
        return 1;
    }
    string start = argv[1];
    int depth = stoi(argv[2]);

    // ensure output dir
    filesystem::create_directories("output");

    auto t0 = chrono::steady_clock::now();
    auto nodes = parallel_bfs(start, depth);
    auto t1 = chrono::steady_clock::now();
    double elapsed = chrono::duration<double>(t1 - t0).count();

    // stdout
    for (auto& n : nodes)
        cout << "- " << n << "\n";
    cout << "Time to crawl: " << elapsed << "s\n";

    // write to file
    string fname = "output/" + start + "_depth" + to_string(depth) + ".txt";
    ofstream ofs(fname);
    for (auto& n : nodes) 
        ofs << n << "\n";
    ofs << "Time: " << elapsed << "s\n";
    ofs.close();

    return 0;
}
