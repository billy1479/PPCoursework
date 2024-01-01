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

void fetchVectorLength(Matrix& x) {
    Vector shortestVector = x[0];
    int NoOfDimensions = x.size();

    float shortestLength = 0.0;

    for (int i = 0; i < NoOfDimensions; i++) {
        shortestLength += shortestVector[i] * shortestVector[i];
    };
    
    std::cout << "Shortest vector length: " << shortestLength << std::endl;
};

Matrix rref(Matrix& x) {
    // performs row operations to convert matrix to row echelon form

    // get rows and columns of matrix
    int m = x.size(); // rows
    int n = x[0].size(); // columns
    int lead = 0; // for indexing

    for (int r = 0; r < m; r++) {
        if (lead >= n) {
            return x;
        };
        int i = r;
        while (x[i][lead] == 0) {
            ++i;
            if (i == m) {
                i = r;
                if(++lead == n) {
                    return x;
                }
            }
        };
        std::swap(x[i],x[r]);
        double lv = x[r][lead];

        for (auto& mrx : x[r]) mrx /= lv;

        for (i = 0; i < m; ++i) {
            if (i != r) {
                lv = x[i][lead];
                for (int j = 0; j < n; ++j) {
                    x[i][j] -= lv * x[r][j];
                }
            }
        }
        ++lead;
    }
    return x;
};

void printMatrix(Matrix& x) {
    std::cout << "The Matrix:" << std::endl;
    for (const auto &row : x) {
        for (const auto &elem : row) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}

double dotProduct(const Vector& a, const Vector& b) {
    double product = 0.0;
    for (int i = 0; i < a.size(); i++)
        product += a[i] * b[i];
    return product;
}

Vector subtract(const Vector& a, const Vector& b) {
    Vector result(a.size());
    for (int i = 0; i < a.size(); i++)
        result[i] = a[i] - b[i];
    return result;
}

Vector multiply(const Vector v, double scalar) {
    Vector result(v.size());
    for (int i = 0; i < v.size(); i++)
        result[i] = v[i] * scalar;
    return result;
}

Vector normalize(const Vector& v) {
    double magnitude = sqrt(dotProduct(v, v));
    return multiply(v, 1.0 / magnitude);
}

Matrix gramSchmidt(const Matrix& vectors) {
    Matrix result;
    for (const auto& v : vectors) {
        Vector temp = v;
        for (const auto& u : result)
            temp = subtract(temp, multiply(u, dotProduct(v, u)));
        result.push_back(normalize(temp));
    }
    return result;
}

// This works but obvioulsy is bad in general
Matrix LLL(Matrix& basis, double delta) {
    size_t n = basis.size();
    Matrix mu(n, Vector(n, 0.0));

    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            mu[i][j] = round(basis[i][j] / basis[j][j]);
            for (size_t k = 0; k < j; ++k) {
                mu[i][j] -= round(mu[i][k] * basis[j][k] / basis[j][j]);
            }
        }

        if (std::abs(basis[i][i] - mu[i][i]) > 0.5) {
            for (size_t k = 0; k <= i - 1; ++k) {
                basis[i][k] -= round(mu[i][k] * basis[i][i] / basis[i][i]);
            }
        }
    }

    for (size_t i = n - 1; i > 0; --i) {
        for (size_t j = i - 1; j < i; --j) {
            double ratio = round(mu[i][j] * basis[i][j] / basis[j][j]);
            for (size_t k = 0; k <= j; ++k) {
                basis[i][k] -= ratio * basis[j][k];
            }
        }
    }
    return basis;
};

Matrix bkz(Matrix& basis) {
    // Number of basis vectors given
    int n = basis.size();

    for (int i = 0; i < n; i++) {
        Matrix basisPrime = LLL(basis, 0.75);
        
    }

    return basis;
};

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

    // Ensures matrix is linearly independent
    // Matrix newMatrix = makeLinearlyIndependent(matrix);

    // Applies gram-schmidt process to matrix
    // Matrix newMatrix2 = gs2(matrix);
    
    // Applies LLL algorithm to matrix
    double delta = 0.75;
    // LLL(newMatrix, delta);

    Matrix newMatrix = gramSchmidt(matrix);

    // Prints new matrix as a result
    printMatrix(newMatrix);

    // Fetches shortest vector
    fetchVectorLength(newMatrix);

    // Outputs it to txt file
    // to be done at the end
    return 0;
}

