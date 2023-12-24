#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>

using namespace std;

int main(int argc, char *argv[]) {
    // Take in arguments and format the vectors as vectors for the program
    // Reads input variables, parses them and stores them in the vectors array
    std::vector<std::vector<int>> vectors;
    for (int i = 1; i < argc; ++i) {
        std::cout << "Argument " << i << ": " << argv[i] << std::endl;
        std::vector<int> currentInputvalues;
        std::string currentValue = argv[i];
        std::stringstream ss(currentValue);

        ss.ignore(1);

        int num;
        while (ss >> num) {
            currentInputvalues.push_back(num);
            ss.ignore(1);
        };
        vectors.push_back(currentInputvalues);
    }


    // Ouptuts the vectors array
    std::cout << "Input strings have been extracted to vectors." << std::endl;
    for (const auto &row : vectors) {
        for (const auto &elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}