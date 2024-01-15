#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <chrono>
#include <sys/time.h>
#include <sys/resource.h>

using namespace std;
using namespace std::chrono;

// Construct matrix for problem
using Vector = vector<double>;
using Matrix = vector<Vector>;

// Outputs matrix
void printMatrix(Matrix x) {
    std::cout << "The Matrix:" << std::endl;
    for (const auto &row : x) {
        for (const auto &elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

// Outputs a vector
void printVector(Vector x) {
    cout << "The Vector:" << endl;
    for (const auto &row : x) {
        cout << row << endl;
    }
}

// Returns the dot product of two given vectors
double dotProduct(const Vector& a, const Vector& b) {
    double product = 0.0;
    for (int i = 0; i < a.size(); i++)
        product += a[i] * b[i];
    return product;
}

// Subtracts two vectors from each other and returns it
Vector subtract(const Vector& a, const Vector& b) {
    Vector result(a.size());
    for (int i = 0; i < a.size(); i++)
        result[i] = a[i] - b[i];
    return result;
}

// Adds two vectors together and returns it
Vector addVector(const Vector& a, const Vector& b) {
    Vector result(a.size());
    for (int i = 0; i < a.size(); i++)
        result[i] = a[i] + b[i];
    return result;
}

// Multiplies a vector and a scalar together and returns it
Vector multiply(const Vector v, double scalar) {
    Vector result(v.size());
    for (int i = 0; i < v.size(); i++)
        result[i] = v[i] * scalar;
    return result;
}

// Multiplies two vectors together and returns the vector
Vector multiplyVectors(Vector v, Vector x) {
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
}

// Returns the euclidian norm squared of a given vector
double eNormSquared(const Vector& v) {
    double sum = 0.0;
    for (double x : v) {
        sum += x * x;
    }
    return sum;
}

// Normalises a given vector
Vector normalize(const Vector &v) {
    double magnitude = sqrt(dotProduct(v, v));
    return multiply(v, 1.0 / magnitude);
}

// Loops through given lattice vectors
// and returns the vector with the shortest euclidean norm
Vector shortestVector(Matrix lattice) {
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

// Makes orthogonal basis via
// the gram-schmdit process (used in LLL)
Matrix gs(Matrix &basis) {
    Matrix result = basis;

    Vector v1 = basis[0];
    result[0] = v1;

    for (int i = 1; i < basis.size(); ++i) {
        Vector vn = basis[i];
        for (int j = 0; j < i; ++j) {
            double temp = (eNorm(result[j]) * eNorm(result[j]));
            double temp2 = dotProduct(basis[i], result[j]);
            vn = subtract(vn, multiply(result[j], temp2 / temp));
        }
        result[i] = vn;
    }
    return result;
}

// Used in the LLL algorithm as part of _________
double mu(Matrix & basis, Matrix & oBasis, int i, int j) {
    float dp = dotProduct(oBasis[j], oBasis[j]);
    if (dp == 0) {
        cerr << "Division by zero error" << endl;
        exit(EXIT_FAILURE);
    } else {
        return (dotProduct(basis[i], oBasis[j]) / dp);
    }
}

bool condition(Matrix& basis, Matrix& oBasis, double delta, int k) {
    double temp;
    temp = abs(mu(basis, oBasis, k, k-1)) * abs(mu(basis, oBasis, k, k-1));
    double temp2;
    temp2 = (dotProduct(oBasis[k-1], oBasis[k-1]));
    double temp3;
    temp3 = dotProduct(oBasis[k], oBasis[k]);

    if (temp3 >= (delta - temp * temp2)) {
        return true;
    } else {
        return false;
    }
}

// LLL algorithm for lattice reduction
Matrix LLL(Matrix& basis, double delta) {
    Matrix oBasis = gs(basis);

    int n = basis.size();
    int k = 1;

    while (k < n) {
        for (int j = k - 1; j >= 0; --j) {
            if (abs(mu(basis, oBasis, k, j)) > 0.5) {
                double temp = round(mu(basis, oBasis, k, j));
                Vector temp2 = multiply(basis[j], temp);
                basis[k] = subtract(basis[k], temp2);
                oBasis = gs(basis);
            }
        }
        if (condition(basis, oBasis, delta, k)) {
            ++k;
        } else {
            basis[k], basis[k-1] = basis[k-1], basis[k];
            oBasis = gs(basis);
            k = max(k-1, 1);
        }
    }
    return basis;
}

// Simple enumeration algorithm that
// is bound by the current known shortest vector
Vector enumerate(const Matrix& basis, Vector shortest) {
    for (int i = 0; i < basis.size(); ++i) {
        double temp = abs(eNorm(addVector(shortest, basis[i])));
        if ((temp < eNorm(shortest)) && temp != 0) {
            shortest = addVector(shortest, basis[i]);
        }
        temp = abs(eNorm(subtract(shortest, basis[i])));
        if ((temp < eNorm(shortest)) && temp != 0) {
            shortest = subtract(shortest, basis[i]);
        }
    }
    return shortest;
}

// LU decomposition and then followed by the determinant
double Determinant(Matrix basis) {
    Matrix lower(basis.size(), Vector(basis.size(),0));
    Matrix upper(basis.size(), Vector(basis.size(), 0));
    for (int i = 0; i < basis.size(); i++) {
            lower[i][i] = 1;
    }
    for (int k = 0; k < basis.size(); ++k) {
        // Initialise lower
        upper[k][k] = basis[k][k];
        for (int i = k+1; i < basis.size(); i++) {
            lower[i][k] = basis[i][k] / upper[k][k];
            upper[k][i] = basis[k][i];
        }
        for (int i = k+1; i < basis.size(); i++) {
            for (int j = k+1; j < basis.size(); j++) {
                basis[i][j] -= lower[i][k] * upper[k][j];
            }
        }
    }
    double det = 1;
    for (int i = 0; i < basis.size(); i++) {
        cout << upper[i][i] << endl;
        det = det * upper[i][i];
    }
    return det;
}

// Main algorithm
int main(int argc, char* argv[]) {
    // For memory usage
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    int NoOfVectors = sqrt(argc - 1);

    if (NoOfVectors * NoOfVectors != argc-1) {
        cerr << "Incorrect number of arguments." << endl;
        exit(EXIT_FAILURE);
    }

    Matrix basis;

    for (int vec = 1; vec < argc; vec += NoOfVectors) {
        Vector currentVector;
        for (int x = 0; x < NoOfVectors; ++x) {
            string arg = argv[vec+x];
            arg.erase(remove(arg.begin(), arg.end(), '['), arg.end());
            arg.erase(remove(arg.begin(), arg.end(), ']'), arg.end());
            try {
                currentVector.push_back(stof(arg));
            } catch (exception& e) {
                cerr << "Incorrect vector numbers" << endl;
            }
        }
        basis.push_back(currentVector);
    }

    double det = 0.0;
    det = Determinant(basis);
    cout << "Determinant: " << det << endl;
    if (det == 0 || isnan(det)) {
        cerr << "Error - input vectors are linearly dependent." << endl;
        exit(EXIT_FAILURE);
    }

    // // Sets start time
    auto startTime = high_resolution_clock().now();

    // // Calls LLL on input basis
    Matrix result = LLL(basis, 0.75);

    // // Takes note of shortest vector from LLL reduced basis
    Vector x = shortestVector(result);
    double shortestNorm = eNorm(x);

    // // Enumerates lattice around LLL-reduced basis
    // to see if any other vectors are shorter than one in the basis
    x = enumerate(result, x);
    shortestNorm = eNorm(x);

    // // Records end time
    auto endTime = high_resolution_clock().now();

    // // Outputs running time of algorithm
    auto duration = duration_cast<microseconds>(endTime - startTime);
    cout << "Run-time for given basis: " << duration.count() << endl;

    // // Writes euclidean norm of shortest vector to text file
    std::ofstream myfile("result.txt");
    if (myfile.is_open()) {
        myfile << shortestNorm;
        myfile.close();
    } else {
        cerr << "Error creating file" << endl;
    }
    // Ends program
    return 0;
}
