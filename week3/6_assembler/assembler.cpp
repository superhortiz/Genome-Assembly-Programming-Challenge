#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <limits>

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

vector<unordered_set<int>> bruijn_graph(vector<string>& k_mers, unordered_map<string, int>& k_minus_1_to_vertex, vector<unordered_map<int, int>>& edge_weights, vector<int>& incoming_edges, vector<int>& outgoing_edges, int& num_edges, const int k) {
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

        edge_weights[prefix_vertex][suffix_vertex]++;

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

void dfs_bubbles(const vector<unordered_set<int>>& graph, const vector<int>& incoming_edges, const vector<int>& outgoing_edges, vector<vector<int>>& paths, vector<int>& path, const int length, const int max_length, const int start_vertex, const int vertex) {
    // Perform DFS until the maximum length is reached

    // Add to the current path the current vertex
    path.push_back(vertex);

    // Add to the dictionary the path for the pair (start_vertex, vertex)
    if (incoming_edges[vertex] > 1 && vertex != start_vertex)
        paths.push_back(path);

    // If the maximum length is reached or there are no more vertices to visit, return
    if (length >= max_length || outgoing_edges[vertex] == 0) {
        path.pop_back();
        return;
    }

    // Recursively explore all the adjacent vertices of the current vertex
    for (const int next_vertex : graph[vertex]) {
        // Discard self loops
        if (next_vertex != vertex)
            dfs_bubbles(graph, incoming_edges, outgoing_edges, paths, path, length + 1, max_length, start_vertex, next_vertex);
    }
    path.pop_back();
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

int get_path_weight(vector<unordered_map<int, int>>& edge_weights, const vector<int>& path) {
    int weight = 0;
    for (int i = 0; i < path.size() - 1; ++i) {
        int v = path[i];
        int w = path[i + 1];
        int edge_weight = edge_weights[v][w];
        weight = max(weight, edge_weight);
    }
    return weight;
}

void delete_path(vector<unordered_set<int>>& graph, const vector<int>& path, vector<int>& incoming_edges, vector<int>& outgoing_edges) {
    for (int i = 0; i < path.size() - 1; ++i) {
        int v = path[i];
        int w = path[i + 1];

        if (graph[v].count(w) != 0)
            graph[v].erase(w);
            incoming_edges[w]--;
            outgoing_edges[v]--;
    }
}

void handle_bubbles(vector<unordered_set<int>>& graph, vector<unordered_map<int, int>>& edge_weights, vector<int>& incoming_edges, vector<int>& outgoing_edges, const int max_length, const int min_coverage) {
    int n = graph.size();

    // Iterate over all the vertices in the graph
    for (int start_vertex = 0; start_vertex < n; ++start_vertex) {
        int num_paths = outgoing_edges[start_vertex];

        // Discard the vertices that don't have branching
        if (num_paths <= 1) continue;

        // Perform DFS to find all the pairs (start_vertex, w),
        // where w is a vertex such incoming_edges[w] > 1
        vector<vector<int>> paths;
        vector<int> path;
        dfs_bubbles(graph, incoming_edges, outgoing_edges, paths, path, 0, max_length, start_vertex, start_vertex);

        for (int i = 0; i < paths.size(); ++i) {
            int path_weight = get_path_weight(edge_weights, paths[i]);

            if (path_weight <= min_coverage) {
                delete_path(graph, paths[i], incoming_edges, outgoing_edges);
            }
        }
    }
}

void dfs_get_contigs(const vector<unordered_set<int>>& graph, const vector<int>& incoming_edges, const vector<int>& outgoing_edges, vector<vector<int>>& paths, const int start_vertex) {
    vector<int> stack, path, levels;
    stack.push_back(start_vertex);
    levels.push_back(0);

    while (!stack.empty()) {
        int vertex = stack.back();
        int level = levels.back();
        stack.pop_back();
        levels.pop_back();
        path.push_back(vertex);

        if (level != 0 && (incoming_edges[vertex] > 1 || outgoing_edges[vertex] > 1)) {
            paths.push_back(path);

            if (!levels.empty()) {
                int prev_level = levels.back();
                for (int i = 0; i < level - prev_level + 1; ++i) {
                    path.pop_back();
                }
            }
            continue;
        }

        for (const int next_vertex : graph[vertex]) {
            stack.push_back(next_vertex);
            levels.push_back(level + 1);
        }
    }
}

vector<string> get_contigs(const vector<unordered_set<int>>& graph, const vector<int>& incoming_edges, const vector<int>& outgoing_edges, const unordered_map<int, string>& vertex_to_k_minus_1, const int k) {
    int n = graph.size();
    vector<vector<int>> contig_paths;

    // Iterate over all the vertices in the graph
    for (int start_vertex = 0; start_vertex < n; ++start_vertex) {
        int num_paths = graph[start_vertex].size();

        // Discard the vertices that don't have branching
        if (num_paths <= 1) continue;

        dfs_get_contigs(graph, incoming_edges, outgoing_edges, contig_paths, start_vertex);
    }

    vector<string> contigs;

    for (const vector<int> contig_path : contig_paths) {
        string contig = vertex_to_k_minus_1.at(contig_path[0]);
        for (int v = 1; v < contig_path.size(); ++v) {
            const string& k_mer = vertex_to_k_minus_1.at(contig_path[v]);
            contig.push_back(k_mer.back());
        }
        contigs.push_back(contig);
    }
    
    return contigs;
}

int find_optimal_k(vector<string>& reads) {
    // This function uses binary search to find the optimal value of k
    int lo = 0.2 * reads[0].size();
    int hi = 0.9 * reads[0].size();
    int optimal_k = -1;

    vector<string> contigs;

    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;

        // Break the reads down by the value k
        vector<string> k_mers = get_k_mers(reads, mid);

        // Create a mapping from unique (k-1)-mers to integer vertex IDs
        unordered_map<string, int> k_minus_1_to_vertex = get_k_minus_1_to_vertex(k_mers, mid);

        // Create the inverse mapping (vertex ID back to (k-1)-mer string)
        unordered_map<int, string> vertex_to_k_minus_1 = get_vertex_to_k_minus_1(k_minus_1_to_vertex);

        // Construct the de Bruijn graph using the k-mers and (k-1)-mer mappings
        int num_edges;
        int num_vertices = k_minus_1_to_vertex.size();
        vector<int> incoming_edges(num_vertices, 0);
        vector<int> outgoing_edges(num_vertices, 0);
        vector<unordered_map<int, int>> edge_weights(num_vertices);
        vector<unordered_set<int>> graph = bruijn_graph(k_mers, k_minus_1_to_vertex, edge_weights, incoming_edges, outgoing_edges, num_edges, mid);

        // Remove tips and handle bubbles
        int min_coverage = 3;
        remove_tips(graph, incoming_edges, outgoing_edges);
        handle_bubbles(graph, edge_weights, incoming_edges, outgoing_edges, mid, min_coverage);

        contigs.clear();
        contigs = get_contigs(graph, incoming_edges, outgoing_edges, vertex_to_k_minus_1, mid);

        if (contigs.size() > 0) {
            optimal_k = mid;
            lo = mid + 1;
        } else {
            hi = mid - 1;
        }
    }
    return optimal_k;
}


int main() {
    // Read data
    int n_reads;
    cin >> n_reads;

    vector<string> reads(n_reads);

    for (int i = 0; i < n_reads; ++i) {
        cin >> reads[i];
    }

    // Define k value
    int k = find_optimal_k(reads);

    // Break the reads down by the value k
    vector<string> k_mers = get_k_mers(reads, k);

    // Create a mapping from unique (k-1)-mers to integer vertex IDs
    unordered_map<string, int> k_minus_1_to_vertex = get_k_minus_1_to_vertex(k_mers, k);

    // Create the inverse mapping (vertex ID back to (k-1)-mer string)
    unordered_map<int, string> vertex_to_k_minus_1 = get_vertex_to_k_minus_1(k_minus_1_to_vertex);

    // Construct the de Bruijn graph using the k-mers and (k-1)-mer mappings
    int num_vertices = k_minus_1_to_vertex.size();
    int num_edges;
    vector<int> incoming_edges(num_vertices, 0);
    vector<int> outgoing_edges(num_vertices, 0);
    vector<unordered_map<int, int>> edge_weights(num_vertices);
    vector<unordered_set<int>> graph = bruijn_graph(k_mers, k_minus_1_to_vertex, edge_weights, incoming_edges, outgoing_edges, num_edges, k);

    // Remove tips and handle bubbles
    int min_coverage = 3;
    remove_tips(graph, incoming_edges, outgoing_edges);
    handle_bubbles(graph, edge_weights, incoming_edges, outgoing_edges, k, min_coverage);

    vector<string> contigs = get_contigs(graph, incoming_edges, outgoing_edges, vertex_to_k_minus_1, k);

    for (int i = 0; i < contigs.size(); ++i) {
        cout << ">CONTIG" << i + 1 << endl;
        cout << contigs[i] << endl;
    }

    return 0;
}
