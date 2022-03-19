#ifndef MATRIX_H
#define MATRIX_H
#include "Vector_.hpp"
#include <cstddef>
#include <iostream>

// forward declaration for Vector class to multiply Matrix obj * Vector obj
class Vector;

class Matrix
{
public:
    Matrix();
    explicit Matrix(int nrow, int ncol);

    Matrix(Matrix const& m) noexcept;
    Matrix& operator=(const Matrix& m) noexcept;

    double** allocateMem_(size_t nrow, size_t ncol);

    Matrix(Matrix&& m) noexcept;
    Matrix& operator=(Matrix&& m) noexcept;

    Matrix(std::initializer_list<std::initializer_list<double> > init);

    virtual ~Matrix();

     Matrix(double** m, int nrow, int ncol) ;
     Matrix(double m[][1], int nrow, int ncol);

     static Matrix Tr_Matrix(Matrix& src);
     static Matrix pinv_Matrix(Matrix& src);

    void setmNumRows(int nrows) { mNumRows = nrows; }
    void setmNumCols(int ncols) { mNumCols = ncols; }

    const int getmNumRows() const {return mNumRows;}
    const int getmNumCols() const {return mNumCols;}

    void setmData(double** data) { mData = data; }
    double** getmData() const { return mData; }

    friend std::ostream& operator<<(std::ostream& os, const Matrix& obj);
    friend std::istream& operator>>(std::istream& is, Matrix& obj);

    double& operator()(int lhs, int rhs) const;
    double& operator()(int lhs, int rhs);

    bool operator[](const Matrix& rhs) const;
    bool operator[](const Matrix& rhs);

    // Note:  only in method
    //    1. return a constant can still change the value after that
    //    2. return constant is the same as return a class in method
    // Arithmetic Operator
    const Matrix operator+(const double& rhs) = delete;
    const Matrix operator-(const double& rhs) = delete;
    const Matrix operator*(const double& rhs);
    const Matrix operator/(const double& rhs);
    const Matrix operator%(const double& rhs) = delete;

    const Matrix operator+(const double& rhs) const = delete;
    const Matrix operator-(const double& rhs) const = delete;
    const Matrix operator*(const double& rhs) const;
    const Matrix operator/(const double& rhs) const = delete;
    const Matrix operator%(const double& rhs) const = delete;

    friend Matrix operator+(const double& rhs, const Matrix& lhs) = delete;
    friend Matrix operator-(const double& rhs, const Matrix& lhs) = delete;
    friend Matrix operator*(const double& rhs, const Matrix& lhs);
    friend Matrix operator/(const double& rhs, const Matrix& lhs) = delete;
    friend Matrix operator%(const double& rhs, const Matrix& lhs) = delete;

    // Matrix multiply Vector and return Vector type
    friend Vector operator*(const Matrix& mat, const Vector& vec);

    const Matrix operator+(const Matrix& rhs) const;
    const Matrix operator-(const Matrix& rhs) const;
    const Matrix operator*(const Matrix& rhs) const;
    const Matrix operator/(const Matrix& rhs) const = delete;
    const Matrix operator%(const Matrix& rhs) const = delete;

    const Matrix operator+(const Matrix& rhs);
    const Matrix operator-(const Matrix& rhs);
    Matrix operator*(const Matrix& rhs);
    const Matrix operator/(const Matrix& rhs) = delete;
    const Matrix operator%(const Matrix& rhs) = delete;

    // Assignment Operators
    const Matrix& operator+=(const Matrix& rhs);
    const Matrix& operator-=(const Matrix& rhs);
    const Matrix& operator*=(const Matrix& rhs);
    const Matrix& operator/=(const Matrix& rhs);
    const Matrix& operator%=(const Matrix& rhs);

    const Matrix& operator+=(const double& rhs) = delete;
    const Matrix& operator-=(const double& rhs) = delete;
    const Matrix operator*=(const double& rhs);
    const Matrix& operator/=(const double& rhs) = delete;
    const Matrix& operator%=(const double& rhs) = delete;

    const Matrix& operator+=(const double& rhs) const = delete;
    const Matrix& operator-=(const double& rhs) const = delete;
    const Matrix operator*=(const double& rhs) const;
    const Matrix& operator/=(const double& rhs) const = delete;
    const Matrix& operator%=(const double& rhs) const = delete;

    // Unary operator
    Matrix operator++() = delete; // ++c
    Matrix operator++(int) = delete; // c++
    Matrix operator--() = delete;
    Matrix operator--(int) = delete;

    const Matrix operator++() const = delete;
    const Matrix operator++(int) const = delete;
    const Matrix operator--() const = delete;
    const Matrix operator--(int) const = delete;

    // Binary operator:
    // Relational operator
    bool operator==(const Matrix& rhs);
    bool operator!=(const Matrix& rhs);
    bool operator<=(const Matrix& rhs);
    bool operator>=(const Matrix& rhs);
    bool operator<(const Matrix& rhs);
    bool operator>(const Matrix& rhs);

    bool operator==(const Matrix& rhs) const;
    bool operator!=(const Matrix& rhs) const;
    bool operator<=(const Matrix& rhs) const;
    bool operator>=(const Matrix& rhs) const;
    bool operator<(const Matrix& rhs) const;
    bool operator>(const Matrix& rhs) const;

    friend bool operator==(const double& lhs, const Matrix& rhs) = delete;
    friend bool operator!=(const double& lhs, const Matrix& rhs) = delete;
    friend bool operator<=(const double& lhs, const Matrix& rhs) = delete;
    friend bool operator>=(const double& lhs, const Matrix& rhs) = delete;
    friend bool operator<(const double& lhs, const Matrix& rhs) = delete;
    friend bool operator>(const double& lhs, const Matrix& rhs) = delete;

    bool operator==(const double& rhs) = delete;
    bool operator!=(const double& rhs) = delete;
    bool operator<=(const double& rhs) = delete;
    bool operator>=(const double& rhs) = delete;
    bool operator<(const double& rhs) = delete;
    bool operator>(const double& rhs) = delete;

    bool operator==(const double& rhs) const = delete;
    bool operator!=(const double& rhs) const = delete;
    bool operator<=(const double& rhs) const = delete;
    bool operator>=(const double& rhs) const = delete;
    bool operator<(const double& rhs) const = delete;
    bool operator>(const double& rhs) const = delete;

    // Logical Operators
    bool operator&&(const Matrix& rhs) const = delete;
    bool operator||(const Matrix& rhs) const = delete;
    bool operator!() const = delete;

    bool operator&&(const Matrix& rhs) = delete;
    bool operator||(const Matrix& rhs) = delete;
    bool operator!() = delete;

    // Also I know how to transform the XOR operator worked as exponential operator i.e, A^T or A^-1
    // But if I try to implement that operator then it leads to wrong purpose of that operator's functionality.
    // Bitwise Operators
    Matrix operator>>(const Matrix& rhs) const = delete;
    Matrix operator<<(const Matrix& rhs) const = delete;
    Matrix operator^(const Matrix& rhs) const = delete;
    Matrix operator&(const Matrix& rhs) const = delete;
    Matrix operator|(const Matrix& rhs) const = delete;
    Matrix operator~() const = delete;
    // But anyway I still implement it if user want to use it
    Matrix operator^(const int& _inverseValue);
    Matrix operator^(const char& T);

    Matrix operator>>(const Matrix& rhs) = delete;
    Matrix operator<<(const Matrix& rhs) = delete;
    Matrix operator^(const Matrix& rhs) = delete;
    Matrix operator&(const Matrix& rhs) = delete;
    Matrix operator|(const Matrix& rhs) = delete;
    Matrix operator~() = delete;

    Matrix operator>>(const int& rhs) const = delete;
    Matrix operator<<(const int& rhs) const = delete;
    //Matrix operator^(const int& rhs) const = delete;
    Matrix operator&(const int& rhs) const = delete;
    Matrix operator|(const int& rhs) const = delete;

    Matrix operator>>(const int& rhs) = delete;
    Matrix operator<<(const int& rhs) = delete;
    //Matrix operator^(const int& rhs) = delete;
    Matrix operator&(const int& rhs) = delete;
    Matrix operator|(const int& rhs) = delete;

    // Determinant of square matrix
    static void getCofactor(const Matrix& src, Matrix& dest);
    void getCofactor(const int& start, const int& idx, const int& n);
    void getCofactor(Matrix& temp, int start, int idx);
    void getCofactor(const Matrix& src, Matrix& temp, int start, int idx, int n);
    double determinantMatix(const Matrix& temp, int n);
    static double determinantMatix(const Matrix& src);
    double determinantMatix();

    // Adjoint of matrix
    void adjointMatrix(Matrix& adj);
    void adjointMatrix(const Matrix& src, Matrix& adj);
    void adjointMatrix();

    // Calculate inverse matrix
    Matrix matrixInverse();
    static void matrixInverse(Matrix& src_mat, Matrix& inv);

    // Check Invertible matrix
    bool isInvertible(Matrix& mat);
    bool isInvertible();

    // Some extra method to swap row<->row and column<->column
    void swapRowMatrix(int row1, int row2, int rank_col);
    static void swapRowMatrix(const Matrix& mat, int row1, int row2, int rank_col);
    static void swapColMatrix(const Matrix& mat, int col1, int col2);
    void swapColMatrix(int col1, int col2);
    void swapRowMatrix(int row1, int row2);

    // display or print the result
    virtual void display() const;
    static void display(const Matrix& src);

    // Calculate Pseudo inverse (Moore-Penrose inverse)
    static void transposeMatrix(const Matrix& src, Matrix& temp);
    Matrix transposeMatrix();
    static int rankMatrix(const Matrix& mat);
    int rankMatrix();
    static void matrixMultiply(const Matrix& mat1, const Matrix& mat2, Matrix& matrix_product);
    void matrixMultiply(const Matrix& mat);
    Matrix pseudoInverse(Matrix& src);
    Matrix pseudoInverse();

    // Pseudo-inverse Tikhonov
    Matrix pseudoInverseTikhonov(const double& lambda);

    // Sub-matrix
    static void setSubMatrix(const Matrix& srcMat_, Matrix& subMat_, std::size_t block_size, std::size_t row, std::size_t col);
    static void setSubMatrix(const Matrix& src, Matrix& subMat_, size_t row, size_t col, size_t mrows, size_t ncols);
    void setSubMatrix(const Matrix& srcMat_, size_t row, size_t col, size_t mrows, size_t ncols);
    void setSubMatrix(const Matrix& srcMat_, std::size_t block_size, std::size_t row, std::size_t col);

    // Sub-matrix for other method
    Matrix setSubMatrix(std::size_t row, std::size_t col, std::size_t mrows, std::size_t ncols);
    Matrix setSubMatrix(size_t block_size, size_t i, size_t j);

    // Print sub-matrix
    void printSubMatrix(std::size_t block_size, std::size_t row, std::size_t col);
    static void printSubMatrix(const Matrix& srcMat_, std::size_t block_size, std::size_t row, std::size_t col);
    static void printSubMatrix(const Matrix& srcMat_, std::size_t mrows, size_t ncols, std::size_t row, std::size_t col);
    void printSubMatrix(std::size_t row, std::size_t col, std::size_t mrows, std::size_t ncols);

    // Free memory
    void deallocateMem(double** mat, int n);

    // Gaussian Elimination with pivoting
    static void free2DArray(double** arr, int row);
    double** create2DArray(int row, int col);
    /* function to reduce matrix to ref, then returns a value to indicate whether matrix is singular or not */
    int forwardElimination(Matrix& mat);
    /* function to calculate the values of the unknowns */
    static void backSubstitution(const Matrix& mat);
    Matrix backSubstitution();
    /* function to get matrix content */
    static void gaussianElimination(Matrix& mat);
    Matrix gaussianElimination();
    /* function for elementary operation of swapping 2 rows */
    static void swapRowMatrix(const Matrix& mat, int row1, int row2);
    /* Identity matrix */
    static void Identity(const Matrix& mat);
    void Identity();
    /* Make augmented matrix */
    Matrix makeAugementedMatrix(const Matrix& mat1, const Matrix& mat2);

    // Conjugate Gradient Method
    bool isSymmetric(const Matrix& mat);
    bool isSymmetric();
    Matrix uminus(Matrix& mat); // unary minus method working similar in matlab
    Matrix uminus();

    // test pseudo inverse
    Matrix geninv();
    static Vector diag(const Matrix& mat, int k = 0);

protected:
    int mNumRows;
    int mNumCols;
    double** mData;

    // Now LinearSystem class can access the private data member of Matrix class because parent class can not access the child class (inheritance case)
    friend class LinearSystem;
    // Now PosSymLinSystem class can access the private data member of Matrix class because parent class can not access the child class (inheritance case)
    friend class PosSymLinSystem;
    // Now Dataset class can access the private data member of Matrix class because parent class can not access the child class (inheritance case)
    friend class Dataset;

};
#endif // MATRIX_H
