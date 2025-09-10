#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <utility>
#include <algorithm>

using namespace std;

int get_overlap(const string& a, const string& b) {
    // This function checks for an exact suffix-prefix overlap between two strings.
    // It determines if a suffix of string 'a' matches a prefix of string 'b'.
    int a_size = a.size();
    int b_size = b.size();

    // Iterate through all possible starting positions for a suffix in string 'a'.
    for (int i = 0; i < a_size; ++i) {
        // Compares suffix from string 'a' starting from i with prefix of suffix 'b'
        for (int j = 0; j < a_size - i; ++j) {
            // if there is a mismatch continue with a new position i
            if (a[i + j] != b[j]) {
                break;
            }

            // If we reach the end of 'a' then return the corresponding overlap
            if (j == a_size - i - 1) {
                return a_size - i;
            }
        }
    }
    // No overlapping was found, then return 0
    return 0;
}

unordered_map<string, set<int>> build_kmer_index(vector<string> &reads, const int k) {
    // Function to get the mapping from all k-mers present in reads
    // to its corresponding index in reads.

    // Map to store the indices of the corresponding k-mer
    unordered_map<string, set<int>> kmer_index;

    // Iterate over all the reads
    for (int read_id = 0; read_id < reads.size(); ++read_id) {
        // Get the current read
        const auto& curr_read = reads[read_id];

        // Get all the possible k-mers in the current read and add the read index to the mapping
        for (int i = 0; i <= curr_read.size() - k; ++i) {
            // Get the current k-mer
            string kmer = curr_read.substr(i, k);

            // If the kmer is new, a new vector is created; otherwise, the read_id is appended.
            kmer_index[kmer].emplace(read_id);
        }
    }
    return kmer_index;
}

set<pair<int, int>> get_candidate_pairs(vector<string> &reads, const int k_mer) {
    // To avoid computing overlaps between all pairs of reads (O(n^2)),
    // first we compute a list of all pairs of reads that share a k-mer.
    // They are the candidates to have an overlapping.

    // Get the mapping from all k-mers present in reads to its corresponding index in reads.
    unordered_map<string, set<int>> kmer_index = build_kmer_index(reads, k_mer);

    // Set to store unique pairs of possible candidates
    set<pair<int, int>> candidate_pairs;

    // Iterate over all the lists of indices obtained for each k-mer
    for (const auto& entry : kmer_index) {
        // Get the current vector that stores the indices
        const auto& curr_read_ids = entry.second;

        // If the vector has size 1, then discard it
        if (curr_read_ids.size() <= 1) continue;

        vector<int> vector_read_ids(curr_read_ids.begin(), curr_read_ids.end());

        // Iterate over all the indices present in the current vector
        // and obtain the candidate pair
        for (int i = 0; i < curr_read_ids.size(); ++i) {
            for (int j = i + 1; j < curr_read_ids.size(); ++j) {
                int read1_id = vector_read_ids[i];
                int read2_id = vector_read_ids[j];

                // Check that the first element is lower or equal than the second,
                // to avoid repetition.
                if (read1_id > read2_id) {
                    swap(read1_id, read2_id);
                }

                // Add the current pair to the set of candidates
                candidate_pairs.emplace(read1_id, read2_id);
            }
        }
    }
    return candidate_pairs;
}

vector<pair<int, int>> build_overlap_graph(vector<string> &reads, const int k_mer) {
    // Function to build the overlap graph. In this case we will use a greedy fashion,
    // only considering the edge with the maximum overlap

    // Get the candidate pairs to avoid computing overlaps between all pairs of reads (O(n^2))
    set<pair<int, int>> candidate_pairs = get_candidate_pairs(reads, k_mer);

    // Vector to store the resulting graph
    vector<pair<int, int>> overlap_graph(reads.size());

    // Iterate over all the candidate pairs
    for (const auto& pair : candidate_pairs) {
        // Extract the vertices
        int u = pair.first;
        int v = pair.second;

        // Get the overlap in both directions
        int overlap_uv = get_overlap(reads[u], reads[v]);
        int overlap_vu = get_overlap(reads[v], reads[u]);

        // Discard in case we dont get an overlap
        if ((overlap_uv <= 0 && overlap_vu <= 0)) continue;

        // Only consider the direction with the biggest overlap
        if (overlap_vu > overlap_uv) {
            swap(u, v);
            swap(overlap_uv, overlap_vu);
        }

        // Just keep the edge with the biggest overlap
        if (overlap_graph[u].second < overlap_uv) {
            overlap_graph[u] = make_pair(v, overlap_uv);
        }
    }

    // Return the resulting graph
    return overlap_graph;
}

string assemble(const vector<pair<int, int>>& graph, const vector<string> &reads, const int start_vertex) {
    // Function to assemble the genome given the overlap graph and corresponding reads.

    // Define starting variables
    string assemble = "";
    int curr_vertex = start_vertex;
    pair<int, int> next_edge = graph[start_vertex];

    // Iterate over the overlap graph until reaching the starting point again
    while (true) {
        // Get the current read
        string curr_read = reads[curr_vertex];

        // Add the corresponding additional characters to the string
        assemble += curr_read.substr(0, curr_read.size() - next_edge.second);
        
        // Update values for the next iteration
        curr_vertex = next_edge.first;
        next_edge = graph[curr_vertex];

        // Check if we have reach the starting vertex, update the string and break the while loop
        if (next_edge.first == start_vertex) {
            string curr_read = reads[curr_vertex];
            assemble += curr_read.substr(0, curr_read.size() - next_edge.second);
            break;
        }
    }

    // Return the resulting assembled genome
    return assemble;
}

int main() {
    // Read and store the input
    vector<string> reads;
    string curr_read;
    while (cin >> curr_read) {
        reads.push_back(curr_read);
    }

    // Get the overlap graph.
    // To avoid computing overlaps between all pairs of reads (O(n^2)),
    // we only consider all pairs of reads that share a k-mer.
    const int k_mer = 12;
    vector<pair<int, int>> graph = build_overlap_graph(reads, k_mer);

    // Given the resulting graph, get the assembled genome
    const int start_vertex = 0;
    cout << assemble(graph, reads, start_vertex) << endl;

    return 0;
}