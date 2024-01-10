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

Matrix gramSchmidt(const Matrix& vectors) {
    Matrix result = vectors;
    int dimension = vectors.size();

    Matrix r(dimension, Vector(dimension));
    Matrix v = vectors;

    // for (int i = 0; i < dimension; ++i) {
    //     // v[i][i] = round(eNorm(v[i]));
    //     // for (int j = 0; j < dimension; ++j) {
    //     //     if (v[i][i] == 0) {
    //     //         cout << "Error - Division by zero (Linearly Dependent vectors)" << endl;
    //     //     } else {
    //     //         result[i][j] = v[i][j] / v[i][i];
    //     //     }
    //     // }
    //     for(int k = i + 1; k < dimension; ++k) {
    //         r[i][k] = dotProduct(result[i], v[k]);
    //         for (int j = 0;j < k;++j) {
    //             // v[k][j] = v[k][j] - (r[i][k] * result[i][j]);
    //             v[k] = subtract(v[k], multiply(v[k], dotProduct(v[k], v[j]) / dotProduct(v[j],v[j])));
    //         }
    //         double normOfColumn = eNorm(v[k]);
    //         for (int j = 0;j<dimension;++j) {
    //             v[k][j] = v[k][j] / normOfColumn;
    //         }
    //     }
    // }

    for (int i = 0;i<dimension;i++) {
        for (int j = 0; j < i; j++) {
            v[i] = subtract(v[i], multiply(v[i], dotProduct(v[i], v[j]) / dotProduct(v[j],v[j])));
        }
        // double normOfColumn = eNorm(v[i]);
        // for (int z = 0;z<dimension;++z) {
        //     v[i][z] = v[i][z] / normOfColumn;
        // }
    }
    return v;
}

Matrix gramSchmidt2(Matrix& vectors) {
    Matrix result = vectors;
    for (int i = 0;i < vectors.size(); i++) {
        Vector currentColumn = vectors[i];
        for (int j = 0; j < i; j++) {
            double x = dotProduct(vectors[j],vectors[i]);
            currentColumn = subtract(currentColumn, multiply(result[j], x));
        }
        double z = eNorm(vectors[i]);
        for (int y = 0; y < vectors.size();y++) {
            result[i][y] = result[i][y] / z;
        }
    }
    return result;
}

// Dont touch this
Matrix gs(Matrix& basis) {
    Matrix result = basis;

    Vector v1 = basis[0];
    v1 = normalize(v1);
    result[0] = v1;

    for (int i = 1; i < basis.size(); i++) {
        Vector vn = basis[i];
        for (int j = 0; j < i; j++) {
            vn = subtract(vn, multiply(result[j], dotProduct(basis[i], result[j])));
        }
        vn = normalize(vn);
        result[i] = vn;
    }

    return result;
}

// Dont touch this
Matrix LLL(Matrix& basis, double delta) {
    Matrix basis_prime = gramSchmidt(basis);
    int index = 1;

    while (index < basis_prime.size()) {
        cout << index << endl;
        for (int j = index - 1; j >= 0; --j) {
            double mu = dotProduct(basis[index], basis[j]);
            if (abs(mu) > 0.5) {
                basis[index] = subtract(basis[index], multiply(basis[j], mu));
                basis_prime = gramSchmidt(basis);
            }
        }
        if (eNorm(basis_prime[index]) >= (delta - pow((dotProduct(basis[index], basis_prime[index-1])),2) * eNorm(basis_prime[index-1]))) {
            index += 1;
        } else {
            basis[index], basis[index-1] = basis[index-1], basis[index];
            basis_prime = gramSchmidt(basis);
            index = std::max(index-1, 1);            
        }
    }
    return basis;
};

Matrix LLL2(Matrix& basis, double delta) {
    Matrix result;


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

    // SUBTRACT works
    // cout << "Test subtract: " << endl;
    // Vector v1 = {2,2,2};
    // Vector v2 = {1,1,1};
    // Vector v3 = subtract(v1, v2);
    // printVector(v3);

    // MULTIPLY works
    // cout << "Test multiply: " << endl;
    // Vector v1 = {1,1,1};
    // int scalar = 3;
    // Vector v2 = multiply(v1, scalar);
    // printVector(v2);

    // Dot product works
    // cout << "Test dot product: " << endl;
    // Vector v1 = {2,7,1};
    // Vector v2 = {8,2,8};
    // double dp = dotProduct(v1, v2);
    // cout << "DP = " << dp << endl;

    // eNorm works
    // cout << "ENORM test: " << endl;
    // Vector v1 = {1,1,1,1};
    // double n = eNorm(v1);
    // cout << n << endl;

    // normalize 
    // cout << "Normalise test: " << endl;
    // Vector v1 = {1,1,1};
    // Vector v2 = normalize(v1);
    // printVector(v2);

    // Applies gram-schmidt process to matrix
    // Matrix newMatrix2 = gramSchmidt(matrix);
    // Matrix gM = gramSchmidt2(matrix);
    // printMatrix(newMatrix2);
    // printMatrix(gM);
    // Matrix newMatrix = gs(matrix);
    // printMatrix(newMatrix);
    Matrix gm2 = gramSchmidt2(matrix);
    printMatrix(gm2);
    // Applies LLL algorithm to matrix
    // double delta = 0.5;
    // LLL(newMatrix, delta);

    // LLL test below - need to sort first 
    // "[1,0,0,1345]" "[0,1,0,35]" "[0,0,1,154]" should be "[0,9,-2,7]" "[1,1,-9,-6]" "[1,-3,-8,8]"
    // Matrix newMatrix = LLL(gM, delta);

    // Vector shortest = Enumeration(newMatrix);
    // printVector(shortest);
    // cout << eNorm(shortest) << endl;

    // Prints new matrix as a result
    // printMatrix(gM);

    // Outputs it to txt file
    // to be done at the end
    return 0;
}

