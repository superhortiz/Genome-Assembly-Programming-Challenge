#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <limits>

using namespace std;

class FlowGraph {
// This class implements a scheme for storing edges of the graph
public:
    struct Edge {
        int from, to, capacity, flow;
    };

    int source, target;

private:
    int n, m, total_lower_bounds;

    // List of all - forward and backward - edges
    vector<Edge> edges;

    // These adjacency lists store only indices of edges in the edges list
    vector<vector<size_t>> graph;

    vector<int> in, out, lower_bounds;

public:
    explicit FlowGraph(size_t n, size_t m): n(n), m(m), source(n), target(n + 1), graph(n + 2), in(n, 0), out(n, 0), total_lower_bounds(0) {}

    void add_edge(const int from, const int to, const int lower_bound, const int capacity) {
        // Store lower bounds
        lower_bounds.push_back(lower_bound);
        total_lower_bounds += lower_bound;

        /* Note that we first append a forward edge and then a backward edge,
         * so all forward edges are stored at even indices (starting from 0),
         * whereas backward edges are stored at odd indices in the list edges */
        Edge forward_edge = {from, to, capacity - lower_bound, 0};

        // For backward edges, keep the capacity equals to zero, then
        // residual_capacity = capacity - flow in both cases.
        // Flow in backward edges is negative.
        Edge backward_edge = {to, from, 0, 0};

        // Add edges to the vectors graph and edges
        graph[from].push_back(edges.size());
        edges.push_back(forward_edge);
        graph[to].push_back(edges.size());
        edges.push_back(backward_edge);

        in[to] += lower_bound;
        out[from] += lower_bound;
    }

    void prepare_for_circulation() {
        // This function completes the graph by adding a "super source" and a "super sink"
        // to handle the lower bounds on edge capacities.
        for (int vertex = 0; vertex < n; ++vertex) {
            // Define edges to connect the source to the vertex
            Edge forward_s2v = {source, vertex, in[vertex], 0};
            Edge backward_s2v = {vertex, source, 0, 0};

            // Define edges to connect the vertex to the target
            Edge forward_v2t = {vertex, target, out[vertex], 0};
            Edge backward_t2v = {target, vertex, 0, 0};

            // Add edges to the vectors graph and edges
            graph[source].push_back(edges.size());
            edges.push_back(forward_s2v);
            graph[vertex].push_back(edges.size());
            edges.push_back(backward_s2v);
            graph[vertex].push_back(edges.size());
            edges.push_back(forward_v2t);
            graph[target].push_back(edges.size());
            edges.push_back(backward_t2v);
        }
    }

    size_t size() const {
        return graph.size();
    }

    const vector<size_t>& get_ids(int from) const {
        return graph[from];
    }

    const Edge& get_edge(size_t id) const {
        return edges[id];
    }

    void add_flow(const size_t id, const int flow) {
        /* To get a backward edge for a true forward edge (id is even), we should get id + 1
         * due to the described above scheme. On the other hand, when we have to get a "backward"
         * edge for a backward edge (get a forward edge for backward - id is odd), id - 1
         * should be taken.
         *
         * It turns out that id ^ 1 works for both cases.*/
        edges[id].flow += flow;
        edges[id ^ 1].flow -= flow;
    }

    vector<int> get_circulation_flows() {
        vector<int> circulation_flows(m);

        for (int edge = 0; edge < m; ++ edge) {
            circulation_flows[edge] = lower_bounds[edge] + edges[2 * edge].flow;
        }

        return circulation_flows;
    }

    int get_total_lower_bounds() {
        return total_lower_bounds;
    }
};

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

FlowGraph get_flow_graph(const vector<vector<int>>& graph, const int num_vertices, const int num_edges) {
    // Read and store the edges

    FlowGraph flow_graph(num_vertices, num_edges);
    for (int v = 0; v < num_vertices; ++v) {
        for (const int& w : graph[v]) {
            flow_graph.add_edge(v, w, 1, 1000);
        }
    }

    // Complete the graph by adding a "super source" and a "super sink"
    flow_graph.prepare_for_circulation();

    // Return the constructed graph
    return flow_graph;
}

void bfs_max_flow(FlowGraph& graph, vector<int>& edge_to, const int from) {
    // Perform BFS to find paths from source to target
    queue<int> queue;
    queue.push(from);
    edge_to[from] = 0;

    while (!queue.empty()) {
        int v = queue.front();
        queue.pop();
        const vector<size_t>& ids = graph.get_ids(v);
        for (const auto& id : ids) {
            const FlowGraph::Edge& edge = graph.get_edge(id);
            if (edge_to[edge.to] == -1 && edge.capacity - edge.flow > 0) {
                queue.push(edge.to);
                edge_to[edge.to] = id; // Get track of the edges
            }
        }
    }
}

int solve_max_flow(FlowGraph& graph) {
    // Initialize vector to store the path found by BFS for each iteration
    size_t n = graph.size();
    vector<int> edge_to(n);

    // Get source and target vertices
    int from = graph.source;
    int to = graph.target;

    // Initialize flow
    int flow = 0;

    while (true) {
        // Reset vector for each iteration
        fill(edge_to.begin(), edge_to.end(), -1);

        // Use BFS to find the shortest path from source (from) to target (to)
        bfs_max_flow(graph, edge_to, from);

        // If there is no path to target, return flow
        if (edge_to[to] == -1) return flow;

        // Find the minimum capacity in the found path
        int x = numeric_limits<int>::max();
        int curr = to;
        while (curr != from) {
            const FlowGraph::Edge& edge = graph.get_edge(edge_to[curr]);
            x = min(x, edge.capacity - edge.flow);
            curr = edge.from;
        }

        // Update the flow for edges in the path
        curr = to;
        while (curr != from) {
            const FlowGraph::Edge& edge = graph.get_edge(edge_to[curr]);
            graph.add_flow(edge_to[curr], x);
            curr = edge.from;
        }

        // Update total flow
        flow += x;
    }
    return flow;
}

vector<int> solve_circulation_flow(FlowGraph& graph) {
    int flow = solve_max_flow(graph);
    int total_lower_bounds = graph.get_total_lower_bounds();
    vector<int> flows = graph.get_circulation_flows();
    return flows;
}

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

list<int> get_eulerian_cycle(vector<vector<int>>& adj, vector<int>& incoming_edges, vector<int>& outgoing_edges, int& num_edges) {
    // Function to find the Eulerian cycle using the Hierholzer's algorithm

    // Find a strongly connected vertex to use as a starting vertex
    int start_vertex = 0;
    while (incoming_edges[start_vertex] == 0 || outgoing_edges[start_vertex] == 0) {
        start_vertex++;
    }

    // Get the first cycle for the starting vertex
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

int main() {
    // Read data
    vector<string> reads;
    string read;
    while (cin >> read) {
        reads.push_back(read);
    }

    // Define k value
    int k = 20;

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
    int min_coverage = 6;
    remove_tips(graph, incoming_edges, outgoing_edges);
    handle_bubbles(graph, edge_weights, incoming_edges, outgoing_edges, k, min_coverage);

    // Transform the graph representation to a vector to allow edge multiplicities.
    // Update the number of edges present in the graph
    num_edges = 0;
    vector<vector<int>> graph_vector(num_vertices);
    for (int v = 0; v < num_vertices; ++v) {
        for (const int& w : graph[v]) {
            graph_vector[v].push_back(w);
            num_edges++;
        }
    }

    // Solve circulation flow
    FlowGraph flow_graph = get_flow_graph(graph_vector, num_vertices, num_edges);
    vector<int> circulation_flow = solve_circulation_flow(flow_graph);
    
    // Based on the circulation flow, construct the complete graph
    vector<vector<int>> complete_graph(num_vertices);
    int i = 0;
    for (int v = 0; v < num_vertices; ++v) {
        for (const int& w : graph_vector[v]) {
            int flow = circulation_flow[i++];

            for (int f = 0; f < flow; ++f) {
                complete_graph[v].push_back(w);
            }
        }
    }

    // Find an Eulerian cycle in the complete graph
    list<int> cycle = get_eulerian_cycle(complete_graph, incoming_edges, outgoing_edges, num_edges);

    // Reconstruct and print the assembled genome
    for (int vertex : cycle) {
        cout << vertex_to_k_minus_1[vertex][0];
    }

    return 0;
}