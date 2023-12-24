#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char *argv[]) {
    int NoOfDimensions;

    for (int i = 0; i < argc; ++i) {
        std::cout << "Argument " << i << ": " << argv[i] << std::endl;

    }

    // Declares the number of dimensions passed 
    int NoOfVectors;
    NoOfVectors = argc - 1;

    int NoOfDimensions;


    return 0;
}