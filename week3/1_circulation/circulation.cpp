#include <iostream>
#include <vector>
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

    void add_edge(int from, int to, int lower_bound, int capacity) {
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

    void add_flow(size_t id, int flow) {
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

FlowGraph read_data() {
    // Read the number of vertices and edges
    int vertex_count, edge_count;
    cin >> vertex_count >> edge_count;

    // Read and store the edges
    FlowGraph graph(vertex_count, edge_count);
    for (int i = 0; i < edge_count; ++i) {
        int u, v, lower_bound, capacity;
        cin >> u >> v >> lower_bound >> capacity;
        graph.add_edge(u - 1, v - 1, lower_bound, capacity);
    }

    // Complete the graph by adding a "super source" and a "super sink"
    graph.prepare_for_circulation();

    // Return the constructed graph
    return graph;
}

void bfs_max_flow(FlowGraph& graph, vector<int>& edge_to, int from) {
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

void solve_circulation_flow(FlowGraph& graph) {
    int flow = solve_max_flow(graph);
    int total_lower_bounds = graph.get_total_lower_bounds();

    if (flow != total_lower_bounds) {
        cout << "NO" << endl;
        return;
    }

    cout << "YES" << endl;
    vector<int> flows = graph.get_circulation_flows();

    for (const int& flow : flows) {
        cout << flow << endl;
    }
}

int main() {
    FlowGraph graph = read_data();
    solve_circulation_flow(graph);
    return 0;
}
