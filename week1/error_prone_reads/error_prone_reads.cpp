#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <utility>
#include <algorithm>

using namespace std;

int get_overlap_with_mismatches(const string& a, const string& b, const int max_num_mismatches) {
    // This function checks for an exact suffix-prefix overlap between two strings.
    // It determines if a suffix of string 'a' matches a prefix of string 'b'.

    int a_size = a.size();
    int b_size = b.size();

    // Iterate through all possible starting positions for a suffix in string 'a'.
    for (int i = 0; i < a_size; ++i) {
        int num_mistaches = 0;
    
        // Compares suffix from string 'a' starting from i with prefix of suffix 'b'
        for (int j = 0; j < a_size - i; ++j) {
            // if there is a mismatch continue with a new position i
            if (a[i + j] != b[j]) {
                num_mistaches++;
            }

            if (num_mistaches > max_num_mismatches) {
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

vector<pair<int, int>> build_overlap_graph(vector<string> &reads, const int k_mer, const int max_num_mistaches) {
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
        int overlap_uv = get_overlap_with_mismatches(reads[u], reads[v], max_num_mistaches);
        int overlap_vu = get_overlap_with_mismatches(reads[v], reads[u], max_num_mistaches);

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

void get_reads_order(const vector<pair<int, int>>& graph, vector<int>& reads_order, vector<int>& overlaps) {
    // Function to get the indices of the reads in the correct graph order.
    // This make easier to iterate through the reads and their direct neighbors.
    // The vector overlaps gives the overlap between the current vertex and the next vertex.

    // Initialize variables and data structures
    int i = 0;
    int start_vertex = 0;
    int curr_vertex = start_vertex;
    reads_order[start_vertex] = start_vertex;

    // While loop to transverse the graph.
    // It will terminate as soon we get the starting vertex.
    while (true) {
        // Get the next vertex based in the current vertex
        pair<int, int> next_edge = graph[curr_vertex];
        int next_vertex = next_edge.first;
        int weigth = next_edge.second;

        // If we have reached the starting vertex, update the overlap of the starting
        // vertex and break the while loop.
        if (next_vertex == start_vertex) {
            overlaps[start_vertex] = weigth;
            break;
        }

        // Update the vectors with the retrieved information
        reads_order[++i] = next_vertex;
        overlaps[i] = weigth;

        // Update the current vertex for the next iteration
        curr_vertex = next_vertex;
    }
}

char boyer_moore_maj_vote(vector<char>& sequence) {
    // Boyer-Moore Majority Vote Algorithm implementation.
    // It works when it is guaranteed that a majority element (element
    // that appears more than half the time in the list) exists.
    int n = sequence.size();

    if (n == 0) return ' ';
    char candidate = sequence[0];
    int count = 1;

    for (int i = 0; i < n; ++i) {
        if (sequence[i] == candidate) {
            count++;
        } else {
            count--;
        }
        if (count == 0) {
            candidate = sequence[i];
            count = 1;
        }
    }
    return candidate;
}

string assemble(vector<int>& reads_order, vector<int>& overlaps, const vector<string>& reads, const int num_voters) {
    // Function to assemble the genome given the vectors: reads_order, overlaps corresponding reads.
    // For each character to add to the assemble, it will use the Boyer-Moore Majority Vote Algorithm.
    // Then, it will consider 'num_voters' reads, where each read is a voter, to determine which character
    // will be the next to add to the sequence.

    string assemble = "";  // Start with an empty string
    int n = reads_order.size();  

    // Iterate over all the reads in order
    for (int i = 0; i < n; ++i) {
        int pos = reads_order[i];
        string curr_read = reads[pos];
        int num_new_chars = curr_read.size() - overlaps[(i + 1) % n];
        vector<char> candidates;

        // Iterate over all the new characters we need to add
        // given the current read
        for (int j = 0; j < num_new_chars; ++j) {
            candidates.clear();
            int char_pos = j;

            // Iterate over all the read voters to get the majority character
            // which will be added to the string
            for (int k = 0; k < num_voters; ++k) {
                int voter_read_pos = (i - k + n) % n;
                string voter_read = reads[reads_order[voter_read_pos]];
                candidates.push_back(voter_read[char_pos]);
                char_pos += voter_read.size() - overlaps[voter_read_pos];

                // If there are no more possible voters, then break the for loop
                if (char_pos >= voter_read.size()) {
                    break;
                }
            }
            // Add the character obtained by Boyer-Moore Majority Vote Algorithm
            assemble.push_back(boyer_moore_maj_vote(candidates));
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
    const int max_num_mistaches = 2;
    vector<pair<int, int>> graph = build_overlap_graph(reads, k_mer, max_num_mistaches);

    // To make easier to iterate through the reads and their direct neighbors
    // we use the function get_reads_order.
    // The vector reads_order will have the indices of reads in order, according to the graph.
    // The vector overlaps gives the overlap between the current vertex and the next vertex.
    int n = graph.size();
    vector<int> reads_order(n);
    vector<int> overlaps(n);
    get_reads_order(graph, reads_order, overlaps);

    // Get the assembled genome. It will consider 'num_voters' different reads to
    // determine the characters present in the assemble.
    int num_voters = 6;
    cout << assemble(reads_order, overlaps, reads, num_voters) << endl;

    return 0;
}