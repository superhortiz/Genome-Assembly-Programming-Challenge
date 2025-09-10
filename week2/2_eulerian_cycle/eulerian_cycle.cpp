#include <iostream>
#include <vector>
#include <list>

using namespace std;

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

list<int> get_cycle(vector<vector<int>>& adj, int& num_edges, int start_vertex) {
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

int main() {
    // Number of vertices (n) and the number of edges (m)
    int n, m;
    cin >> n >> m;

    // Read graph's edges
    vector<vector<int>> adj(n, vector<int>());
    vector<int> incoming_edges(n);
    vector<int> outgoing_edges(n);

    for (int i = 0; i < m; ++i) {
        int u, v;
        cin >> u >> v;
        adj[u - 1].push_back(v - 1);
        incoming_edges[v - 1]++;
        outgoing_edges[u - 1]++;
    }

    // Verify if the graph satisfies Euler's Theorem
    bool has_cycle = eulerian_cycle_exist(n, incoming_edges, outgoing_edges);
    cout << has_cycle << endl;

    // If it does not satisfy Euler's Theorem, return
    if (!has_cycle) return 0;

    // Get the Eulerian cycle and print the results
    list<int> cycle = get_eulerian_cycle(adj, m);
    for (int vertex : cycle) {
        cout << vertex + 1 << ' ';
    }

    return 0;
}