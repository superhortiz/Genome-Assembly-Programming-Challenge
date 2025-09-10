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

void dfs(const vector<unordered_set<int>>& graph, const vector<int>& incoming_edges, const vector<int>& outgoing_edges, unordered_map<int, vector<vector<int>>>& paths, vector<int> path, const int length, const int max_length, const int start_vertex, const int vertex) {
    // Perform DFS until the maximum length is reached

    // Add to the current path the current vertex
    path.push_back(vertex);

    // Add to the dictionary the path for the pair (start_vertex, vertex)
    if (incoming_edges[vertex] > 1 && vertex != start_vertex)
        paths[vertex].push_back(path);

    // If the maximum length is reached or there are no more vertices to visit, return
    if (length >= max_length || outgoing_edges[vertex] == 0)
        return;

    // Recursively explore all the adjacent vertices of the current vertex
    for (const int next_vertex : graph[vertex]) {
        // Discard self loops
        if (next_vertex != vertex)
            dfs(graph, incoming_edges, outgoing_edges, paths, path, length + 1, max_length, start_vertex, next_vertex);
    }
}

bool are_disjoint(const vector<int>& path1, const vector<int>& path2) {
    // Function to check if two paths are disjoint
    if (path1.size() > path2.size()) {
        return are_disjoint(path2, path1);
    }

    unordered_set<int> set_path1(path1.begin(), path1.end());

    for (int i = 1; i < path2.size() - 1; ++i) {
        int vertex = path2[i];
        if (set_path1.count(vertex) > 0)
            return false;
    }
    return true;
}

int find_bubbles(const vector<unordered_set<int>>& graph, const vector<int>& incoming_edges, const vector<int>& outgoing_edges, const int max_length) {
    int num_bubbles = 0;
    int n = graph.size();

    // Iterate over all the vertices in the graph
    for (int start_vertex = 0; start_vertex < n; ++start_vertex) {
        int num_paths = outgoing_edges[start_vertex];

        // Discard the vertices that don't have branching
        if (num_paths <= 1) continue;

        // Perform DFS to find all the pairs (start_vertex, w), where w is a 
        // vertex such incoming_edges[w] > 1
        int length = 0;
        vector<int> path;
        unordered_map<int, vector<vector<int>>> paths;
        dfs(graph, incoming_edges, outgoing_edges, paths, path, length, max_length, start_vertex, start_vertex);

        // Iterate over all the pairs (start_vertex, w) found
        for (const auto& v_w_pair : paths) {
            vector<vector<int>> directed_paths = v_w_pair.second;
            int num_paths = directed_paths.size();

            // Compare all the possible combination of paths and check if they are disjoint
            for (int i = 0; i < num_paths; ++i) {
                for (int j = i + 1; j < num_paths; ++j) {
                    if (are_disjoint(directed_paths[i], directed_paths[j]))
                        // If they are disjoint, increment the count of bubbles found
                        num_bubbles++;
                }
            }
        }
    }

    // Return the number of bubbles found in the graph
    return num_bubbles;
}

int main() {
    // Read the input data
    int k, t;
    cin >> k >> t;

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

    // Get the number of bubbles present in the graph
    cout << find_bubbles(graph, incoming_edges, outgoing_edges, t) << endl;

    return 0;
}