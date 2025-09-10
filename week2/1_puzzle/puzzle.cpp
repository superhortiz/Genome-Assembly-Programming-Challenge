#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

vector<string> split_string(string input_string, char delimiter) {
    // Function to split a string given a delimiter
    vector<string> tokens;
    int start_index = 1;
    int n = input_string.size();

    for (int i = 1; i < n; ++i) {
        if (input_string[i] == delimiter || i == n - 1) {
            tokens.push_back(input_string.substr(start_index, i - start_index));
            start_index = i + 1;
        }
    }
    return tokens;
}

vector<vector<vector<string>>> get_solution(const vector<vector<string>>& input, const int n, const int len) {
    // This function uses brute force to find a solution

    // Vector to store the solution
    vector<vector<vector<string>>> solution(n, vector<vector<string>>(n, vector<string>(len)));

    // Vector to store groups of pieces, given their location
    vector<vector<string>> upper;
    vector<vector<string>> lower;
    vector<vector<string>> left;
    vector<vector<string>> right;
    vector<vector<string>> central;

    for (const auto& el : input) {
        // Store directly in the solution the upper left corner
        if (el[0] == "black" && el[1] == "black") {
            solution[0][0] = el;
        
        // Store directly in the solution the upper right corner
        } else if (el[0] == "black" && el[3] == "black") {
            solution[0][n - 1] = el;

        // Store directly in the solution the lower left corner
        } else if (el[1] == "black" && el[2] == "black") {
            solution[n - 1][0] = el;

        // Store directly in the solution the lower right corner
        } else if (el[2] == "black" && el[3] == "black") {
            solution[n - 1][n - 1] = el;

        // Store pieces of the upper border
        } else if (el[0] == "black") {
            upper.push_back(el);

        // Store pieces of the left border
        } else if (el[1] == "black") {
            left.push_back(el);

        // Store pieces of the lower border
        } else if (el[2] == "black") {
            lower.push_back(el);

        // Store pieces of the right border
        } else if (el[3] == "black") {
            right.push_back(el);

        // // Store pieces of the center
        } else {
            central.push_back(el);
        } 
    }

    // Get solution of the upper border
    vector<int> indices = {0, 1, 2};
    do {
        vector<string>& piece1 = upper[indices[0]];
        vector<string>& piece2 = upper[indices[1]];
        vector<string>& piece3 = upper[indices[2]];

        if (solution[0][0][3] == piece1[1] && piece1[3] == piece2[1] && piece2[3] == piece3[1] && piece3[3] == solution[0][n-1][1]) {
            solution[0][1] = piece1;
            solution[0][2] = piece2;
            solution[0][3] = piece3;
            break;
        }
    } while (next_permutation(indices.begin(), indices.end()));

    // Get solution of the lower border
    indices = {0, 1, 2};
    do {
        vector<string>& piece1 = lower[indices[0]];
        vector<string>& piece2 = lower[indices[1]];
        vector<string>& piece3 = lower[indices[2]];

        if (solution[n - 1][0][3] == piece1[1] && piece1[3] == piece2[1] && piece2[3] == piece3[1] && piece3[3] == solution[n-1][n-1][1]) {
            solution[n - 1][1] = piece1;
            solution[n - 1][2] = piece2;
            solution[n - 1][3] = piece3;
            break;
        }
    } while (next_permutation(indices.begin(), indices.end()));

    // Get solution of the left border
    indices = {0, 1, 2};
    do {
        vector<string>& piece1 = left[indices[0]];
        vector<string>& piece2 = left[indices[1]];
        vector<string>& piece3 = left[indices[2]];

        if (solution[0][0][2] == piece1[0] && piece1[2] == piece2[0] && piece2[2] == piece3[0] && piece3[2] == solution[n-1][0][0]) {
            solution[1][0] = piece1;
            solution[2][0] = piece2;
            solution[3][0] = piece3;
        }
    } while (next_permutation(indices.begin(), indices.end()));

    // Get solution of the right border
    indices = {0, 1, 2};
    do {
        vector<string>& piece1 = right[indices[0]];
        vector<string>& piece2 = right[indices[1]];
        vector<string>& piece3 = right[indices[2]];

        if (solution[0][n - 1][2] == piece1[0] && piece1[2] == piece2[0] && piece2[2] == piece3[0] && piece3[2] == solution[n-1][n-1][0]) {
            solution[1][n - 1] = piece1;
            solution[2][n - 1] = piece2;
            solution[3][n - 1] = piece3;
        }
    } while (next_permutation(indices.begin(), indices.end()));

    // Get solution of the central pieces
    indices = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    do {
        vector<string>& piece1 = central[indices[0]];
        vector<string>& piece2 = central[indices[1]];
        vector<string>& piece3 = central[indices[2]];
        vector<string>& piece4 = central[indices[3]];
        vector<string>& piece5 = central[indices[4]];
        vector<string>& piece6 = central[indices[5]];
        vector<string>& piece7 = central[indices[6]];
        vector<string>& piece8 = central[indices[7]];
        vector<string>& piece9 = central[indices[8]];

        if (solution[1][0][3] == piece1[1] && piece1[3] == piece2[1] && piece2[3] == piece3[1] && piece3[3] == solution[1][n-1][1] && 
            solution[2][0][3] == piece4[1] && piece4[3] == piece5[1] && piece5[3] == piece6[1] && piece6[3] == solution[2][n-1][1] && 
            solution[3][0][3] == piece7[1] && piece7[3] == piece8[1] && piece8[3] == piece9[1] && piece9[3] == solution[3][n-1][1] && 

            solution[0][1][2] == piece1[0] && piece1[2] == piece4[0] && piece4[2] == piece7[0] && piece7[2] == solution[n-1][1][0] && 
            solution[0][2][2] == piece2[0] && piece2[2] == piece5[0] && piece5[2] == piece8[0] && piece8[2] == solution[n-1][2][0] && 
            solution[0][3][2] == piece3[0] && piece3[2] == piece6[0] && piece6[2] == piece9[0] && piece9[2] == solution[n-1][3][0]) {

            solution[1][1] = piece1;
            solution[1][2] = piece2;
            solution[1][3] = piece3;
            solution[2][1] = piece4;
            solution[2][2] = piece5;
            solution[2][3] = piece6;
            solution[3][1] = piece7;
            solution[3][2] = piece8;
            solution[3][3] = piece9;
        }
    } while (next_permutation(indices.begin(), indices.end()));

    return solution;
}

int main() {
    const int n = 5;
    const int len = 4;
    vector<vector<string>> input;
    string piece;

    while (cin >> piece) {
        vector<string> colors = split_string(piece, ',');
        input.push_back(colors);
    }

    vector<vector<vector<string>>> solution = get_solution(input, n, len);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << '(';
            for (int k = 0; k < len; ++k) {
                cout << solution[i][j][k];
                if (k < len - 1) {
                    cout << ',';
                }
            }
            cout << ')';
            if (j < n - 1) {
                cout << ';';
            }
        }
        if (i < n - 1) {
            cout << endl;
        }
    }
}