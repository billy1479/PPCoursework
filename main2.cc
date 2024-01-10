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

double u(Matrix& basis, Matrix& oBasis,int i, int k) {
    return (dotProduct(basis[k], oBasis[i]) / dotProduct(basis[i], basis[i]));
}

Matrix step1ofLLL(Matrix& basis) {
    Matrix oBasis = gs(basis);
    for (int i = 0; i < basis.size(); ++i) {
        cout << "i: " << i << endl;
        for (int k = i - 1; k > -1; k--) {
            cout << "k: " << k << endl; 
            double m = round(u(basis, oBasis, i, k));
            basis[i] = subtract(basis[i], multiply(basis[k], m));
        }
    }
    return basis, oBasis;
}

Matrix LLL(Matrix& basis, float delta) {
    Matrix newBasis, oBasis = step1ofLLL(basis);
    for (int i = 0; i < newBasis.size() - 1; ++i) {
        if (eNormSquared(addVector(oBasis[i + 1], multiply(oBasis[i], u(newBasis, oBasis, i, i + 1)))) < (0.75 * eNormSquared(oBasis[i]))) {
            newBasis[i+1], newBasis[i] = newBasis[i], newBasis[i+1];
            newBasis, oBasis = step1ofLLL(newBasis);
            // cout << "This is running" << endl;
        }
        cout << "this part is running" << endl;
    }
    return basis;
}

Vector reduce(Vector& v, Matrix& basis) {
    for (auto& l : basis) {
        int counter = 0;
        while (bool state = true and counter <= basis.size()) {
            Vector difference = subtract(v, l);
            cout << "run" << endl;
            if (eNorm(difference) >= eNorm(v)) {
                state = false;
            }
            v = difference;
            counter++;
        }
    }
    return v;
}

// Maybe this is better
Matrix sieve(Matrix& basis, int maxIt) {
    Matrix tempMatrix;
    int iterations = 0;

    for (auto& x : basis) {
        if (iterations >= maxIt) {
            return tempMatrix;
        }
        x = reduce(x, basis);
        if (eNorm(x) == 0) {
        } else {
            iterations++;
            tempMatrix.push_back(x);
        }
    }

    return tempMatrix;
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

    // Generates all vectors in latice based off basis
    // Matrix newBasis = gs(matrix);

    // printMatrix(newBasis);

    // Matrix lattice = enumumer(newBasis, 10);

    // Matrix lattice = sieve(matrix,10);

    Matrix result = LLL(matrix, 0.75);

    printMatrix(result);

    // Vector x = shortestVector(lattice);
    // double shortestNorm = eNorm(x);
    // cout << "The shortest vector is: " << endl;
    // printVector(x);

    // cout << "The norm of the shortest vector in the lattice is " << shortestNorm << endl;


    // Creates output file and writes euclidean norm of shortest vector to it
    // ofstream myfile;
    // myfile.open("result.txt");
    // myfile << shortestNorm;
    // myfile.close();

    // Ends program
    return 0;
}