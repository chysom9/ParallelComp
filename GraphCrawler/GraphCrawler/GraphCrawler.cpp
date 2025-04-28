#include <rapidjson/document.h>
#include <rapidjson/reader.h>
#include <iostream>
#include <string>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <curl/curl.h>
#include <chrono>
using namespace std;

string replaceSpaces(const string &s)
{
    string result;
    for (char c : s)
        result += (c == ' ') ? string("%20") : string(1, c);
    return result;
}
// take the string passes in
// if the s`tring has spaces in it, change it to '%20', if not we gucci.

// storing API output

// Callback function for handling API response data
// This function is called when data is received from the API
size_t apiRequest(void *contents, size_t size, size_t nmemb, std::string *output)
{
    // Calculate the total number of bytes of data received
    size_t total = size * nmemb;

    // Append the  data content to the  string.
    // contents is a pointer to the raw data received, and total is the number of bytes.
    output->append(static_cast<char *>(contents), total);

    // Return the total number of bytes processed.
    return total;
}

// get neighbors from API
std::vector<std::string> getNeighbors(const std::string &node)
{
    // Encode the node name by replacing spaces with URL encoding.
    std::string encodedNode = replaceSpaces(node);

    // Construct the URL using the encoded node name.
    std::string url = "http://hollywood-graph-crawler.bridgesuncc.org/neighbors/" + encodedNode;

    // Initialize a CURL session.
    CURL *curl = curl_easy_init();
    std::vector<std::string> neighbors;
    // I saw this next 4 lines on stack overflow because i was kinda lost
    // but the rest was mainly ME and some youtube tutorials.
    if (curl)
    {
        std::string response;

        // Set the URL and the write callback function for CURL.
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, apiRequest);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);

        // Perform the CURL request and check for errors.
        CURLcode res = curl_easy_perform(curl);
        if (res != CURLE_OK)
        {
            // If there's an error, print it.
            std::cerr << "CURL Error: " << curl_easy_strerror(res) << std::endl;
        }
        else
        {
            // If the request is successful, parse the JSON response.
            rapidjson::Document doc;
            doc.Parse(response.c_str());

            // Check if the JSON is valid and contains the "neighbors" field.
            if (!doc.HasParseError() && doc.HasMember("neighbors"))
            {
                // Loop through the array of neighbors and add them to the neighbors vector.
                for (const auto &neighbor : doc["neighbors"].GetArray())
                {
                    neighbors.push_back(neighbor.GetString());
                }
            }
            else
            {
                // If JSON is invalid or doesn't contain the expected data, print an error.
                std::cerr << "Errors: Invalid JSON response" << std::endl;
            }
        }

        // Clean up the CURL session.
        curl_easy_cleanup(curl);
    }

    // Return the list of neighbors.
    return neighbors;
}
// BFS
void crawl(const string &rootNode, int maxDepth)
{
    auto startClock = chrono::high_resolution_clock::now();

    queue<pair<string, int>> q;
    unordered_set<string> exploredNodes;

    q.push({rootNode, 0});
    exploredNodes.insert(rootNode);

    while (!q.empty())
    {
        auto [currentNode, currentDepth] = q.front();
        q.pop();

        cout << " Node: " << currentNode << " at Depth " << currentDepth << endl;

        if (currentDepth < maxDepth)
        {
            vector<string> adjacentNodes = getNeighbors(currentNode);
            for (const auto &adjacentNode : adjacentNodes)
            {
                if (exploredNodes.find(adjacentNode) == exploredNodes.end())
                {
                    exploredNodes.insert(adjacentNode);
                    q.push({adjacentNode, currentDepth + 1}); // Queue the adjacent node at a deeper level
                }
            }
        }
    }

    auto endClock = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsedTime = endClock - startClock;
    cout << "Total time for traversal: " << elapsedTime.count() << " seconds." << endl;
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cerr << "Incorrect usage. Expected format: ./crawler <starting_node> <maximum_depth>\n";
        return 1;
    }

    string startingNode = argv[1];
    int maxDepth = stoi(argv[2]);

    crawl(startingNode, maxDepth);

    return 0;
}
