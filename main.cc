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
    Matrix basis_prime = gramSchmidt(basis);
    int index = 1;

    while (index < basis_prime.size()) {
        for (int j = index - 1; j >= -1; j--) {
            double mu = round(dotProduct(basis[index], basis_prime[j]) / dotProduct(basis_prime[j], basis_prime[j]));
            if (mu != 0) {
                basis[index] = subtract(basis[index], multiply(basis[j], mu));
                basis_prime = gramSchmidt(basis);
            }
        }
        if (normalize(basis_prime[index]) >= (delta - multiply((dotProduct(basis[index], basis_prime[index - 1]) * dotProduct(basis[index], basis_prime[index - 1]), normalize(basis_prime[index - 1]))))) {
            index += 1;
        } else {
            basis[index], basis[index-1] = basis[index-1], basis[index];
            basis_prime = gramSchmidt(basis);
            index = std::max(index-1, 1);            
        }
    }

    return basis;
};

Matrix bkz(Matrix& basis) {
    // Number of basis vectors given
    int n = basis.size();
    Matrix basisPrime = LLL(basis, 0.75);

    for (int i = 0; i < n; i++) {
        
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

    // LLL test below - need to sort first 
    // "[1,0,0,1345]" "[0,1,0,35]" "[0,0,1,154]" should be "[0,9,-2,7]" "[1,1,-9,-6]" "[1,-3,-8,8]"
    Matrix newMatrix = LLL(matrix, delta);

    // Prints new matrix as a result
    printMatrix(newMatrix);

    // Fetches shortest vector
    // fetchVectorLength(newMatrix);

    // Outputs it to txt file
    // to be done at the end
    return 0;
}

