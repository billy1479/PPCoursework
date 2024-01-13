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

void printMatrix(Matrix x) {
    std::cout << "The Matrix:" << std::endl;
    for (const auto &row : x) {
        for (const auto &elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

void printVector(Vector x) {
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
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] + b[i];
    }
    return result;
}

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

// Vector projectionFunction(Vector& )

// Makes orthogonal basis
Matrix gs(Matrix basis) {
    Matrix result = basis;

    Vector v1 = basis[0];
    // v1 = normalize(v1);
    result[0] = v1;

    for (int i = 1; i < basis.size(); i++) {
        Vector vn = basis[i];
        for (int j = 0; j < i; j++) {
            double temp = (eNorm(result[j]) * eNorm(result[j]));
            double temp2 = dotProduct(basis[i], result[j]);
            vn = subtract(vn, multiply(result[j], temp2 / temp));
        }
        // vn = normalize(vn);
        result[i] = vn;
    }
    return result;
}

// makes orthonormal basis
Matrix gs2(Matrix basis) {
    Matrix result = basis;

    Vector v1 = basis[0];
    v1 = normalize(v1);
    result[0] = v1;

    for (int i = 1; i < basis.size(); i++) {
        Vector vn = basis[i];
        for (int j = 0; j < i; j++) {
            double temp = (eNorm(result[j]) * eNorm(result[j]));
            double temp2 = dotProduct(basis[i], result[j]);
            vn = subtract(vn, multiply(result[j], temp2 / temp));
        }
        vn = normalize(vn);
        result[i] = vn;
    }
    return result;
}

double mu(Matrix basis, Matrix oBasis, int i, int j) {
    return (dotProduct(basis[i], oBasis[j]) / dotProduct(oBasis[j], oBasis[j]));
}

bool condition(Matrix basis, Matrix oBasis, double delta, int k) {
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

Matrix LLL(Matrix basis, double delta) {
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

int main(int argc, char* argv[]) {
    
    int NoOfVectors = sqrt(argc - 1);

    if (NoOfVectors * NoOfVectors != argc-1) {
        cerr << "Incorrect number of arguments." << endl;
        exit(EXIT_FAILURE);
    }
    Matrix basis;

    for (int vec = 1; vec < argc; vec += NoOfVectors) {
        Vector nbasis;
        cout << "Vector:" << endl;
        for (int x = 0; x < NoOfVectors; ++x) {
            string arg = argv[vec+x];
            arg.erase(remove(arg.begin(), arg.end(), '['), arg.end());
            arg.erase(remove(arg.begin(), arg.end(), ']'), arg.end());
            cout << arg << endl;
            try {
                float a = stof(arg);
                nbasis.push_back(stof(arg));
            } catch (exception& e) {
                cerr << "Incorrect vector numbers" << endl;
            }
        }
        basis.push_back(nbasis);
    }
    printMatrix(basis);
    // cout << "Dimension of matrix: " << dimension << endl;

    ofstream myfile("./result.txt");

    Matrix result = LLL(basis, 0.75);

    printMatrix(result);

    Vector x = shortestVector(result);
    double shortestNorm = eNorm(x);
    printVector(x);

    // Creates output file and writes euclidean norm of shortest vector to it
    myfile << shortestNorm;
    myfile.close();

    // Ends program
    return 0;
}
