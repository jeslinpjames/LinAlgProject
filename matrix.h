// Assignment
// [Due Date: 04 September 2023 Midnight]
// Implement a C++ class named Matrix to represent a matrix. The file matrix.h shall contain the
// implementation. Following implementation should be attempted.
// Private members
// nb_row and nb_col store the number of rows and the number of columns of the matrix.
// mat_data is a std::vector<double>, it stores the matrix data, the ij values.
// Constructors
// The class has 5 constructors. A constructor to generate random matrices (uniformly or normaly),
// a constructor to load data from a CSV file and create a matrix, a constructor to create a matrix
// using an initializer list or a vector of vectors.
// Overloaded operators
// Multiplication *, *=
// Addition +, +=
// Subtraction -, -=
// Division \
// Equal ==
// Difference !=
// Member functions
// column and row to get a column or a row of a matrix.
// sub_matrix returns a submatrix of a matrix.
// shape print the dimension of a matrix.
// reshape reshape (change the number of rows and columns) of a matrix.
// add_row and add_column add a new column or a new row to a matrix.
// remove_column delete a column of a matrix.
// reorder_column sort a matrix column, or flat matrix.
// sort_matrix sort a matrix by column, using indexes
// T() and transpose (a friend function): return the transpose of the matrix.
// Id create a unitary matrix.
// sum returns the sum a flattened matrix.
// avg compute the average of a flattened matrix.
// head print first row of matrix
// print formatted print of a matrix, a value, a string.
// to_csv save a matrix in a csv file.
// Non Member Functions
// diag_matrix returns a diagonal matrix.
// is_triangular Check if a matrix is triangular
// LU decompostion of a matrix using the basic LU algorithm
// QR_factorization Matrix factorization using QR-factorization
// eigen_decomposition eigen decompostion of a matrix using the basic QR algorithm
// SVD Compute the singular value decomposition of a matrix using the QR algorithm

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <initializer_list>
#include <random> 
#include <fstream>   
#include <sstream> 
using namespace std;

class Matrix {
private:
    size_t r;
    size_t c;
    vector<double> mat_data;

public:
    // Constructors
    Matrix(size_t rows, size_t cols);
    Matrix(const vector<vector<double>>& data);
    Matrix(initializer_list<initializer_list<double>> data);
    Matrix(size_t rows, size_t cols, bool useUniform, double mean = 0.0, double stddev = 1.0);
    Matrix(const string& filename);
    // Overloaded operators
    Matrix operator*(const Matrix& other) const;
    Matrix& operator*=(const Matrix& other);
    Matrix operator+(const Matrix& other) const;
    Matrix& operator+=(const Matrix& other);
    Matrix operator-(const Matrix& other) const;
    Matrix& operator-=(const Matrix& other);
    Matrix operator/(double scalar) const;
    bool operator==(const Matrix& other) const;
    bool operator!=(const Matrix& other) const;

    // Member functions
    vector<double> column(size_t col) const;
    vector<double> row(size_t row) const;
    Matrix sub_matrix(size_t start_row, size_t start_col, size_t end_row, size_t end_col) const;
    pair<size_t, size_t> shape() const;
    void reshape(size_t new_rows, size_t new_cols);
    void add_row(const vector<double>& new_row);
    void add_column(const vector<double>& new_col);
    void remove_column(size_t x);
    void remove_row(size_t x);
    void reorder_column(size_t col, const vector<size_t>& order);
    Matrix T() const ;
    friend Matrix transpose(const Matrix& matrix);
    void print() const;
    void to_csv(const string& filename) const;

    // Non-member functions
    static Matrix diag_matrix(const vector<double>& diag);
    bool is_triangular() const;
    // Other matrix factorization functions

    // Other utility functions

};

Matrix::Matrix(size_t rows, size_t cols){
    r = rows;  
    c = cols;  
}
Matrix::Matrix(const vector<vector<double>>& data){
    r = data.size();  
    if(r>0){ 
        c = data[0].size(); 
        mat_data.reserve(r*c);
        for(size_t i =0;i<data.size();++i){
            const vector<double>& row = data[i];
            if(row.size()!=c){  
                throw invalid_argument("Input matrix rows have different sizes");
            }
            mat_data.insert(mat_data.end(),row.begin(),row.end());
        }
    }else
    c =0;  
}
Matrix:: Matrix(initializer_list<initializer_list<double>> data){
   r = data.size();  
    if(r>0){  
         c = data.begin()->size();  
         mat_data.reserve(r*c);  
         for(const auto& row : data){
              if(row.size()!=c){  
                throw invalid_argument("Input matrix rows have different sizes");
              }
              mat_data.insert(mat_data.end(),row.begin(),row.end());
         } 
    } else{
        c =0;  
    }
}
Matrix::Matrix(size_t rows, size_t cols, bool useUniform, double mean, double stddev) {
    r = rows;  
    c = cols;  

    random_device rd;
    mt19937 generator(rd());

    if (useUniform) {
        uniform_real_distribution<double> uniformDistribution(0.0, 1.0);

        mat_data.reserve(r * c);  
        for (size_t i = 0; i < r * c; ++i) {  
            mat_data.push_back(uniformDistribution(generator));
        }
    } else {
        normal_distribution<double> normalDistribution(mean, stddev);

        mat_data.reserve(r * c);  
        for (size_t i = 0; i < r * c; ++i) {  
            mat_data.push_back(normalDistribution(generator));
        }
    }
}
Matrix::Matrix(const string& filename) {
    ifstream file(filename);  // Open the file

    if (!file.is_open()) {
        throw runtime_error("Failed to open CSV file");
    }

    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        double value;
        vector<double> rowValues;

        while (iss >> value) {
            rowValues.push_back(value);
            if (iss.peek() == ',') {
                iss.ignore();
            }
        }

        if (!rowValues.empty()) {
            if (c == 0) {  
                c = rowValues.size();  
            } else if (rowValues.size() != c) {  
                throw runtime_error("CSV rows have different sizes");
            }
            mat_data.insert(mat_data.end(), rowValues.begin(), rowValues.end());
            r++;  
        }
    }
}

Matrix Matrix:: operator*(const Matrix& other) const{
    if(c!=other.r){  
        throw invalid_argument("Matrix dimensions are not compatible");
    }
    Matrix result(r,other.c);  
    for(size_t i =0;i<r;++i){
        for(size_t j =0;j<other.c;++j){
            double sum =0.0;
            for(size_t k =0;k<c;++k){
                sum+=mat_data[i*c+k]*other.mat_data[k*other.c+j];
            }
            result.mat_data[i*other.c+j] = sum;
        }
    }
    return result;
}
Matrix& Matrix:: operator*=(const Matrix& other){
    *this = *this * other;
    return *this;
}
Matrix Matrix:: operator+(const Matrix& other) const{
    if(r!=other.r || c!=other.c){  
        throw invalid_argument("Matrix dimensions are not compatible");
    }
    Matrix result(r,c);  
    for(size_t i =0;i<r;++i){
        for(size_t j =0;j<c;++j){
            result.mat_data[i*c+j] = mat_data[i*c+j] + other.mat_data[i*c+j];
        }
    }
    return result;
}
Matrix& Matrix:: operator+=(const Matrix& other){
    *this = *this + other;
    return *this;
}
Matrix Matrix:: operator-(const Matrix& other) const{
    if(r!=other.r || c!=other.c){  
        throw invalid_argument("Matrix dimensions are not compatible");
    }
    Matrix result(r,c);  
    for(size_t i =0;i<r;++i){
        for(size_t j =0;j<c;++j){
            result.mat_data[i*c+j] = mat_data[i*c+j] - other.mat_data[i*c+j];
        }
    }
    return result;
}
 Matrix& Matrix:: operator-=(const Matrix& other){
    *this = *this - other;
    return *this;
}
Matrix Matrix:: operator/(double scalar) const{
    if (scalar == 0.0) {
        throw runtime_error("Division by zero");
    }
    Matrix result(r,c);  
    for(size_t i =0;i<r;++i){
        for(size_t j =0;j<c;++j){
            result.mat_data[i*c+j] = mat_data[i*c+j] / scalar;
        }
    }
    return result;
}
bool Matrix:: operator==(const Matrix& other) const{
    if(r!=other.r || c!=other.c){  
        return false;
    }
    for(size_t i =0;i<r;++i){
        for(size_t j =0;j<c;++j){
            if(mat_data[i*c+j] != other.mat_data[i*c+j]){
                return false;
            }
        }
    }
    return true;
}
bool Matrix:: operator!=(const Matrix& other) const{
    return !(*this == other);
}
vector<double> Matrix:: column(size_t col) const{
    if (col >= c) {
        throw out_of_range("Column index out of range");
    }
    vector<double> columnData;
    columnData.reserve(r);
    for (size_t i = 0; i < r; ++i) {
        columnData.push_back(mat_data[i * c + col]);
    }
}
vector<double> Matrix :: row(size_t row) const{
    if (row >= r) {
        throw out_of_range("Row index out of range");
    }
    vector<double> rowData;
    rowData.reserve(c);
    for (size_t i = 0; i < c; ++i) {
        rowData.push_back(mat_data[row * c + i]);
    }
}

Matrix Matrix ::sub_matrix(size_t start_row, size_t start_col, size_t end_row, size_t end_col) const{
    if (start_row >= r || start_col >= c || end_row >= r || end_col >= c) {
        throw out_of_range("Submatrix indices are out of range");
    }
    
    size_t sub_rows = end_row - start_row + 1;
    size_t sub_cols = end_col - start_col + 1;
    Matrix submatrix(sub_rows, sub_cols);
    for (size_t i = start_row; i <= end_row; ++i) {
        for (size_t j = start_col; j <= end_col; ++j) {
            submatrix.mat_data.push_back(mat_data[i * c + j]);
        }
    }

    return submatrix;
}
pair<size_t, size_t> Matrix:: shape() const{
    return make_pair(r,c);
}
void Matrix:: reshape(size_t new_rows, size_t new_cols){
    if (r * c != new_rows * new_cols) {
        throw invalid_argument("New shape must be compatible with the original shape");
    }
    r = new_rows;
    c = new_cols;
}
void Matrix:: add_row(const vector<double>& new_row){
    if (new_row.size() != c) {
            throw invalid_argument("New row size does not match matrix's number of columns");
        }
    mat_data.insert(mat_data.end(), new_row.begin(), new_row.end());
    r++;    
}
void Matrix::add_column(const vector<double>& new_col){
    if (new_col.size() != r) {
            throw invalid_argument("New column size does not match matrix's number of rows");
        }
    for (size_t i = 0; i < r; ++i) {
            mat_data.insert(mat_data.begin() + (i * (c + 1)), new_col[i]);
        }
    c++;    
}
void Matrix:: remove_column(size_t x){
    if (x >= c) {
        throw out_of_range("Column index out of range");
    }
    for (size_t i = 0; i < r; ++i) {
            mat_data.erase(mat_data.begin() + (i * c) + x);
        }
    c--;    
}
void Matrix::remove_row(size_t x){
    if (x >= r) {
        throw out_of_range("Row index out of range");
    }
    mat_data.erase(mat_data.begin() + (x * c), mat_data.begin() + ((x + 1) * c));
    r--;    
}
void Matrix::reorder_column(size_t col, const vector<size_t>& order) {
    if (col >= c) {
        throw out_of_range("Column index out of range");
    }
    if (order.size() != r) {
        throw invalid_argument("Order vector size does not match matrix's number of rows");
    }
    vector<double> reorderedColumn;
    reorderedColumn.reserve(r);
    for (size_t i = 0; i < r; ++i) {
        reorderedColumn.push_back(mat_data[i * c + col]);
    }
    for (size_t i = 0; i < r; ++i) {
        mat_data[i * c + col] = reorderedColumn[order[i]];
    }
}
 Matrix Matrix::T() const {
    Matrix transposed(c, r);  // Create a new matrix with swapped dimensions
    for (size_t i = 0; i < r; ++i) {
        for (size_t j = 0; j < c; ++j) {
            transposed.mat_data[j * r + i] = mat_data[i * c + j]; // Transpose the elements
        }
    }
    return transposed;
}

Matrix transpose(const Matrix& matrix) {
    Matrix transposed(matrix.c, matrix.r); // Create a new matrix with swapped dimensions
    for (size_t i = 0; i < matrix.r; ++i) {
        for (size_t j = 0; j < matrix.c; ++j) {
            transposed.mat_data[j * matrix.r + i] = matrix.mat_data[i * matrix.c + j]; // Transpose the elements
        }
    }
    return transposed;
}

// static Matrix Matrix::Id(size_t size){
    
// }
#endif // MATRIX_H
