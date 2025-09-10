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

bool eulerian_cycle_exist(int n, vector<int>& incoming_edges, vector<int>& outgoing_edges) {
    // Function to verify if the graph satisfies Euler's Theorem.
    // It is guaranteed that the graph is strongly connected,
    // then this function just checks if the graph is balanced.
    for (int i = 0; i < n; ++i) {
        if (incoming_edges[i] != outgoing_edges[i]) {
            return false;
        }
    }
    return true;
}

int find_optimal_k(vector<string>& reads) {
    // This function uses binary search to find the optimal value of k
    int lo = 2;
    int hi = reads[0].size();
    int optimal_k = -1;

    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;

        vector<string> k_mers = get_k_mers(reads, mid);

        // Create a mapping from unique (k-1)-mers to integer vertex IDs
        unordered_map<string, int> k_minus_1_to_vertex = get_k_minus_1_to_vertex(k_mers, mid);

        // Construct the de Bruijn graph using the k-mers and (k-1)-mer mappings
        int num_edges;
        int n = k_minus_1_to_vertex.size();
        vector<int> incoming_edges(n, 0);
        vector<int> outgoing_edges(n, 0);
        vector<unordered_set<int>> graph = bruijn_graph(k_mers, k_minus_1_to_vertex, incoming_edges, outgoing_edges, num_edges, mid);

        // If an Eulerian cycle exists for the current 'k' (mid),
        // this 'k' is a potential optimal value. Try a larger 'k' in the upper half.
        if (eulerian_cycle_exist(n, incoming_edges, outgoing_edges)) {
            optimal_k = mid;
            lo = mid + 1;
        } else {
            // If no Eulerian cycle exists, 'k' is too large. Search in the lower half.
            hi = mid - 1;
        }
    }

    return optimal_k;
}

int main() {
    // Read the inputs
    vector<string> reads;
    string read;

    while (cin >> read) {
        reads.push_back(read);
    }

    // Find the optimal value of k
    cout << find_optimal_k(reads);
    return 0;
}