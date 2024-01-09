#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>

using namespace std;

// Construct matrix for problem
using Vector = vector<long double>;
using Matrix = vector<Vector>;

// This is baseline approach

void printMatrix(Matrix& x) {
    std::cout << "The Matrix:" << std::endl;
    for (const auto &row : x) {
        for (const auto &elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

void printVector(Vector& x) {
    cout << "The Vector:" << endl;
    for (const auto &row : x) {
        cout << row << endl;
    }
}

double dotProduct(const Vector& a, const Vector& b) {
    double product = 0.0;
    for (int i = 0; i < a.size(); i++)
        product += a[i] * b[i];
    return product;
}

// Subtracts two vectors from each other
Vector subtract(const Vector& a, const Vector& b) {
    Vector result(a.size());
    for (int i = 0; i < a.size(); i++)
        result[i] = a[i] - b[i];
    return result;
}

// Multiplies a vector and a scalar together
Vector multiply(const Vector v, double scalar) {
    Vector result(v.size());
    for (int i = 0; i < v.size(); i++)
        result[i] = v[i] * scalar;
    return result;
}

Vector multiplyVectors(Vector& v, Vector& x) {
    Vector product(v.size());

    for (size_t i = 0; i < v.size(); ++i) {
        product[i] = v[i] * x[i];
    }

    return product;
}

// Returns the euclidean norm of a vector
double eNorm(const Vector& v) {
    double sum = 0.0;
    for (double x : v) {
        sum += x * x;
    }
    return sqrt(sum);
};

Vector normalize(const Vector& v) {
    double magnitude = sqrt(dotProduct(v, v));
    return multiply(v, 1.0 / magnitude);
}

Vector combineVectors(const Vector& a, const Vector& b) {
    Vector result(a.size());
    for (int i = 0; i < a.size();i++) {
        result[i] = a[i] + b[i];
    }
    return result;
}

Matrix gramSchmidt(Matrix& basis) {
    int i,j,k;
    Matrix r(basis.size());
    Matrix v(basis.size());

    for (i = 0; i < basis.size(); i++) {
        for (j = 0; j < basis.size(); j++) {
            v[i][j] = basis[i][j];
        }
    }

    for (i = 0; i < basis.size(); i++) {
        r[i][i] = eNorm(v[i]);
        for (j = 0; j < basis.size(); j++) {
            basis[i][j] = v[i][j] / r[i][i];
        }
        for (k = i + 1; k < basis.size(); k++) {
            r[i][k] = dotProduct(basis[i], v[k]);
            for (j = 0; j < basis.size(); j++) {
                v[k][j] = v[k][j] - r[i][k] * basis[i][j];
            }
        }
    }
    return v;
}

int main(int argc, char *argv[]) {
    // Take in arguments and format the vectors as vectors for the program
    // Reads input variables, parses them and stores them in the vectors array

    // Sets a variable for dimension of vectors, to see if they differ as this will be an invalid input
    int dimension;
    dimension = 0;

    Matrix matrix;
    for (int i = 1; i < argc; ++i) {
        std::cout << "Argument " << i << ": " << argv[i] << std::endl;
        Vector currentInputvalues;
        std::string currentValue = argv[i];
        std::stringstream ss(currentValue);

        // ignores the [
        ss.ignore(1);

        int tempDimension;
        tempDimension = 0;
        int num;

        if (dimension == 0) {
            while (ss >> num) {
                dimension += 1;
                currentInputvalues.push_back(num);
                ss.ignore(1);
            };
        } else {
            while (ss >> num) {
                tempDimension += 1;
                currentInputvalues.push_back(num);
                ss.ignore(1);
            };
        };

        // Checks if the current vector matches the dimensions of others and if not it errors out
        if (tempDimension == dimension || tempDimension == 0) {
            matrix.push_back(currentInputvalues);
        } else {
            std::cerr << "The vectors given are an inconsistent number of dimensions" << std::endl;
        }        
    }
    std::cout << "Dimension of matrix: " << dimension << std::endl;

    Matrix x = gramSchmidt(matrix);
    printMatrix(x);
    // Creates output file and writes euclidean norm of shortest vector to it
    // ofstream myfile;
    // myfile.open("result.txt");
    // myfile << shortestNorm;
    // myfile.close();

    // Ends program
    return 0;
}