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

struct lllreturn {
    Matrix a;
    Matrix b;
};

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

// Adds two vectors together
Vector addVector(const Vector& a, const Vector& b) {
    Vector result(a.size());
    for (int i = 0; i < a.size(); i++)
        result[i] = a[i] + b[i];
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

double eNormSquared(const Vector& v) {
    double sum = 0.0;
    for (double x : v) {
        sum += x * x;
    }
    return sum;
}

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

// Generates all possible vectors based off given basis
Vector Enumeration(Matrix& basis) {
    Vector shortest = basis[0];

    for (int i = 0; i <basis.size();i++){
        for (int j = 0; j < basis.size(); j++) {
            Vector n = combineVectors(basis[i], basis[j]);
            if (eNorm(n) < eNorm(shortest)) {
                shortest = n;
            }
        }
    }
    return shortest;
}

Matrix enumumer(Matrix& basis, int limits) {
    Matrix lattice;

    for (int i = -limits; i <= limits; i++) {
        for (int j = -limits; j <= limits; j++) {
            Vector vec(basis.size(), 0.0);

            for (size_t k = 0; k < vec.size(); k++) {
                vec[k] = i * basis[0][k] + j * basis[1][k];
            }
            lattice.push_back(vec);
        }
    }
    return lattice;
}

Vector shortestVector(Matrix& lattice) {
    double minNorm = numeric_limits<double>::max();
    Vector shortestVec;

    for (const auto& vec : lattice) {
        double currentNorm = eNorm(vec);
        if (currentNorm < minNorm && currentNorm != 0) {
            minNorm = currentNorm;
            shortestVec = vec;
        }
    }
    return shortestVec;
}

// Vector projectionFunction(Vector& )

// Makes orthogonal basis
Matrix gs(Matrix& basis) {
    Matrix result = basis;

    Vector v1 = basis[0];
    // v1 = normalize(v1);
    result[0] = v1;

    for (int i = 1; i < basis.size(); i++) {
        Vector vn = basis[i];
        for (int j = 0; j < i; j++) {
            vn = subtract(vn, multiply(result[j], dotProduct(basis[i], result[j]) / (eNorm(result[j]) * eNorm(result[j]))));
        }
        // vn = normalize(vn);
        result[i] = vn;
    }
    return result;
}

// makes orthonormal basis
Matrix gs2(Matrix& basis) {
    Matrix result = basis;

    Vector v1 = basis[0];
    v1 = normalize(v1);
    result[0] = v1;

    for (int i = 1; i < basis.size(); i++) {
        Vector vn = basis[i];
        for (int j = 0; j < i; j++) {
            vn = subtract(vn, multiply(result[j], dotProduct(basis[i], result[j]) / (eNorm(result[j]) * eNorm(result[j]))));
        }
        vn = normalize(vn);
        result[i] = vn;
    }
    return result;
}

double mu(Matrix& basis, Matrix& oBasis, int i, int j) {
    return (dotProduct(basis[i], oBasis[j]) / dotProduct(oBasis[j], oBasis[j]));
}

Matrix LLL(Matrix basis, double delta) {
    Matrix oBasis = gs(basis);

    int n = basis.size();
    int k = 1;

    while (k < n) {
        for (int j = k - 1; j >= 0; --j) {
            if (abs(mu(basis, oBasis, k, j)) > 0.5) {
                basis[k] = subtract(basis[k], multiply(basis[j], round(mu(basis, oBasis,k,j))));
                oBasis = gs(basis);
            }
        }
        
        if (dotProduct(oBasis[k], oBasis[k]) >= ((delta - abs(mu(basis, oBasis, k, k-1))*(delta - abs(mu(basis, oBasis, k, k-1)))*(dotProduct(oBasis[k-1], oBasis[k-1]))))) {
            ++k;
        } else {
            basis[k], basis[k-1] = basis[k-1], basis[k];
            oBasis = gs(basis);
            k = max(k-1, 1);
        }
    }
    return basis;
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

    Matrix result = LLL(matrix, 0.75);

    printMatrix(result);

    Vector x = shortestVector(result);
    double shortestNorm = eNorm(x);
    cout << "The shortest vector is: " << endl;
    printVector(x);

    cout << "The norm of the shortest vector in the lattice is " << shortestNorm << endl;


    // Creates output file and writes euclidean norm of shortest vector to it
    ofstream myfile;
    myfile.open("result.txt");
    myfile << shortestNorm;
    myfile.close();

    // Ends program
    return 0;
}