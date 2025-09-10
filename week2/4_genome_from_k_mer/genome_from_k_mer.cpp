#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <unordered_map>

using namespace std;

list<int> get_cycle(vector<vector<int>>& adj, int& num_edges, const int start_vertex) {
    // This function returns a cycle, that is found when it returns
    // to the starting vertex. Every edge used in the cycle is removed from
    // the adjacency list and the number of edges present in the adjacency list is updated.

    // Use a double linked list to store the cycle
    list<int> cycle;

    // Start traversing the graph from start_vertex
    int curr_vertex = start_vertex;

    // Loop to traverse the cycle
    while (true) {
        // Extract the next vertex to visit and update the
        // adjacency list and number of edges.
        int next_vertex = adj[curr_vertex].back();
        adj[curr_vertex].pop_back();
        cycle.emplace_back(next_vertex);
        num_edges--;

        // The loop finishes when we return to start_vertex
        if (next_vertex == start_vertex) break;

        // Update current vertex
        curr_vertex = next_vertex;
    }

    // Return the found cycle
    return cycle;
}

list<int> get_eulerian_cycle(vector<vector<int>>& adj, int num_edges) {
    // Function to find the Eulerian cycle using the Hierholzer's algorithm

    // Start in vertex 0 and get a cycle
    int start_vertex = 0;
    list<int> cycle = get_cycle(adj, num_edges, start_vertex);

    // Iterate until we do not have more edges to visit
    while (num_edges > 0) {
        // Iterate over the vertices present in the current cycle
        auto it = cycle.begin();
        while (it != cycle.end()) {
            int vertex = *it;

            // Check if we have a vertices with unvisited outgoing edges
            if (adj[vertex].size() != 0) {
                // Get a new cycle starting from the vertex we have found
                list<int> new_cycle = get_cycle(adj, num_edges, vertex);

                // Join the new cycle with the current cycle
                cycle.splice(next(it), new_cycle, new_cycle.begin(), new_cycle.end());
                break;
            }
            ++it;
        }
    }
    // Return the Eulerian cycle
    return cycle;
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

vector<vector<int>> bruijn_graph(vector<string>& k_mers, unordered_map<string, int>& k_minus_1_to_vertex, int& num_edges, const int k) {
    // Constructs the de Bruijn graph from a given set of k-mers and (k-1)-mer mappings

    // Determine the total number of vertices in the graph
    int n_vert = k_minus_1_to_vertex.size();

    // Initialize the adjacency list for the graph
    vector<vector<int>> adj(n_vert);

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
        adj[prefix_vertex].push_back(suffix_vertex);

        // Increment the total count of edge
        num_edges++;
    }

    // Return the constructed adjacency list
    return adj;
}

int main() {
    vector<string> k_mers;
    string k_mer;

    while (cin >> k_mer) {
        k_mers.push_back(k_mer);
    }

    int k = k_mers[0].size();

    // Create a mapping from unique (k-1)-mers to integer vertex IDs
    unordered_map<string, int> k_minus_1_to_vertex = get_k_minus_1_to_vertex(k_mers, k);

    // Create the inverse mapping (vertex ID back to (k-1)-mer string)
    unordered_map<int, string> vertex_to_k_minus_1 = get_vertex_to_k_minus_1(k_minus_1_to_vertex);

    // Construct the de Bruijn graph using the k-mers and (k-1)-mer mappings
    int num_edges;
    vector<vector<int>> graph = bruijn_graph(k_mers, k_minus_1_to_vertex, num_edges, k);

    // Find an Eulerian cycle in the constructed de Bruijn graph
    list<int> cycle = get_eulerian_cycle(graph, num_edges);

    // Reconstruct and print the assembled genome
    for (int vertex : cycle) {
        cout << vertex_to_k_minus_1[vertex][0];
    }

    return 0;
}