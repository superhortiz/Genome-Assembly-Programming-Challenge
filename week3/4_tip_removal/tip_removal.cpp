#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <unordered_map>
#include <unordered_set>

using namespace std;

vector<string> get_k_mers(vector<string>& reads, int k) {
    // This function breaks the reads down by the value k
    vector<string> k_mers;

    for (const string& read : reads) {
        int n = read.size();

        for (int i = 0; i <= n - k; ++i) {
            k_mers.push_back(read.substr(i, k));
        }
    }
    return k_mers;
}

unordered_map<string, int> get_k_minus_1_to_vertex(const vector<string>& k_mers, const int k) {
    // Generates a mapping from unique (k-1)-mers to unique integer vertex IDs.
    unordered_map<string, int> k_minus_1_to_vertex;
    int vertex = -1;

    for (const string& k_mer : k_mers) {
        string prefix = k_mer.substr(0, k - 1);
        string suffix = k_mer.substr(1, k - 1);

        // Assign a new unique ID to the prefix if not already encountered
        if (k_minus_1_to_vertex.count(prefix) == 0) {
            k_minus_1_to_vertex[prefix] = ++vertex;
        }

        // Assign a new unique ID to the suffix if not already encountered
        if (k_minus_1_to_vertex.count(suffix) == 0) {
            k_minus_1_to_vertex[suffix] = ++vertex;
        }
    }

    return k_minus_1_to_vertex;
}

unordered_map<int, string> get_vertex_to_k_minus_1(const unordered_map<string, int>& k_minus_1_to_vertex) {
    // Generates the inverse mapping: from unique integer vertex IDs back to their
    // corresponding (k-1)-mer string representations.
    unordered_map<int, string> vertex_to_k_minus_1;

    for (const auto& el : k_minus_1_to_vertex) {
        int vertex = el.second;
        string k_minus_1_mer = el.first;
        vertex_to_k_minus_1[vertex] = k_minus_1_mer;
    }

    return vertex_to_k_minus_1;
}

vector<unordered_set<int>> bruijn_graph(vector<string>& k_mers, unordered_map<string, int>& k_minus_1_to_vertex, vector<int>& incoming_edges, vector<int>& outgoing_edges, int& num_edges, const int k) {
    // Constructs the de Bruijn graph from a given set of k-mers and (k-1)-mer mappings

    // Determine the total number of vertices in the graph
    int n_vert = k_minus_1_to_vertex.size();

    // Initialize the adjacency list for the graph
    vector<unordered_set<int>> adj(n_vert);

    // Initialize edge count to 0 (since it's passed by reference to be updated)
    num_edges = 0; 

    // Iterate through each k-mer to build the graph edges
    for (const string& k_mer : k_mers) {
        // Extract the (k-1)-mer prefix and suffix of the current k-mer
        string prefix = k_mer.substr(0, k - 1);
        string suffix = k_mer.substr(1, k - 1);

        // Get the integer IDs for the prefix and suffix nodes
        int prefix_vertex = k_minus_1_to_vertex[prefix];
        int suffix_vertex = k_minus_1_to_vertex[suffix];

        // Add a directed edge from the prefix vertex to the suffix vertex
        if (adj[prefix_vertex].count(suffix_vertex) == 0) {
            adj[prefix_vertex].emplace(suffix_vertex);
            outgoing_edges[prefix_vertex]++;
            incoming_edges[suffix_vertex]++;

            // Increment the total count of edge
            num_edges++;
        }
    }

    // Return the constructed adjacency list
    return adj;
}

int remove_tips(vector<unordered_set<int>>& graph, vector<int>& incoming_edges, vector<int>& outgoing_edges) {
    // Function to remove all the tip from the given graph

    vector<unordered_set<int>> reversed_graph(graph.size());
    vector<pair<int, int>> starting_tips;
    vector<pair<int, int>> ending_tips;

    // Iterate over all the vertices of the graph to reverse the graph
    // and identify tips.
    for (int v = 0; v < graph.size(); ++v) {
        for (const auto& w : graph[v]) {
            reversed_graph[w].emplace(v);

            if (incoming_edges[v] == 0) {
                starting_tips.push_back({v, w});
            }

            if (outgoing_edges[w] == 0) {
                ending_tips.push_back({v, w});
            }
        }
    }

    int tips_deleted = 0;

    // Iterate until we remove all the starting tips and their branches
    while (!starting_tips.empty()) {
        pair<int, int> edge = starting_tips.back();
        int u = edge.first;
        int v = edge.second;
        starting_tips.pop_back();

        if (graph[u].count(v) != 0) {
            graph[u].erase(v);
            incoming_edges[v]--;
            outgoing_edges[u]--;
            tips_deleted++;

            // Add new tips to the list to continue the process
            for (const int& w : graph[v]) {
                if (incoming_edges[v] == 0)
                    starting_tips.push_back({v, w});
            }
        }
    }

    // Iterate until we remove all the ending tips and their branches
    while (!ending_tips.empty()) {
        pair<int, int> edge = ending_tips.back();
        int v = edge.first;
        int w = edge.second;
        ending_tips.pop_back();

        if (reversed_graph[w].count(v) != 0) {
            graph[v].erase(w);
            reversed_graph[w].erase(v);
            incoming_edges[w]--;
            outgoing_edges[v]--;
            tips_deleted++;

            // Add new tips to the list to continue the process
            for (const int& u : reversed_graph[v]) {
                if (outgoing_edges[v] == 0)
                    ending_tips.push_back({u, v});
            }
        }
    }

    // Return the total number of tips deleted
    return tips_deleted;
}

int main() {
    // Define k value
    int k = 15;

    vector<string> reads;
    string read;

    while (cin >> read) {
        reads.push_back(read);
    }

    // Break the reads down by the value k
    vector<string> k_mers = get_k_mers(reads, k);

    // Create a mapping from unique (k-1)-mers to integer vertex IDs
    unordered_map<string, int> k_minus_1_to_vertex = get_k_minus_1_to_vertex(k_mers, k);

    // Create the inverse mapping (vertex ID back to (k-1)-mer string)
    unordered_map<int, string> vertex_to_k_minus_1 = get_vertex_to_k_minus_1(k_minus_1_to_vertex);

    // Construct the de Bruijn graph using the k-mers and (k-1)-mer mappings
    int num_edges;
    int n = k_minus_1_to_vertex.size();
    vector<int> incoming_edges(n, 0);
    vector<int> outgoing_edges(n, 0);
    vector<unordered_set<int>> graph = bruijn_graph(k_mers, k_minus_1_to_vertex, incoming_edges, outgoing_edges, num_edges, k);

    cout << remove_tips(graph, incoming_edges, outgoing_edges);

    return 0;
}