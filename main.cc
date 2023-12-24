#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]) {
    int NoOfDimensions;
    std::vector<std::vector<int>> vectors;
    for (int i = 1; i < argc; ++i) {
        std::cout << "Argument " << i << ": " << argv[i] << std::endl;
        std::istringstream iss(argv[i]);
        std::vector<int> currentInputvalues;
        char x;
        iss >> x;

        while (iss >> x) {
            if (x == ',' || x == ']') {
                continue;
            }
            int value;
            if (iss >> value) {
                currentInputvalues.push_back(value);
            } else {
                std::cerr << "Error reading the input string." << std::endl;
                return 1;
            }
        }
    }
    std::cout << "Input strings have been extracted to vectors.";
    for (const std::vector<int>& array : vectors) {
        for (int value : array) {
            std::cout << value << ' ';
        }
        std::cout << std::endl;
    }

    return 0;
}