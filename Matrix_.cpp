#include "Matrix_.hpp"
#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <time.h>
#include <initializer_list>
#include <cmath>
#include <assert.h>
#include <utility>
#include <iomanip>

Matrix::Matrix()
    : mNumRows(0)
    , mNumCols(0)
    , mData(nullptr) {}

Matrix::Matrix(int nrow, int ncol)
{
    try
    {
        if(nrow < 0 || ncol < 0)
            throw std::logic_error("Error: Negative rows or columns of matrix is not allowed");
        else{
            mNumRows = nrow;
            mNumCols = ncol;
            srand(time(nullptr));
            mData = new double*[mNumRows];

            for(int i = 0; i < mNumRows; ++i)
                mData[i] = new double[mNumCols];

            for(int i = 0; i < mNumRows; ++i){
                for(int j = 0; j < mNumCols; ++j){
                    mData[i][j] = 0;
                    //mData[i][j] = rand() % 100;
                }
            }
        }
    }

    catch(std::logic_error& le)
    {
        std::cout << le.what() << '\n';
        exit(EXIT_FAILURE);
    }
}

Matrix::Matrix(double** m, int nrow, int ncol)
    : mNumRows(nrow), mNumCols(ncol), mData(nullptr)
{
    mData = allocateMem_(mNumRows, mNumCols);

    for(int i = 0; i < nrow; ++i)
        for(int j = 0; j < ncol; ++j)
            mData[i][j] = m[i][j];
}

Matrix::Matrix(double m[][1], int nrow, int ncol)
    : mNumRows(nrow), mNumCols(ncol), mData(nullptr)
{
    mData = allocateMem_(mNumRows, mNumCols);

    for(int i = 0; i < mNumRows; ++i)
        for(int j = 0; j < mNumCols; ++j)
            mData[i][j] = m[i][j];
}

Matrix Matrix::Tr_Matrix(Matrix& src)
{
    Matrix dummy(src.mData, src.mNumCols, src.mNumRows);

    for(int i = 0; i < dummy.mNumRows; ++i)
        for(int j = 0; j < dummy.mNumCols; ++j)
            dummy.mData[i][j] = src.mData[j][i];

    return dummy;
}

Matrix Matrix::pinv_Matrix(Matrix& src)
{
    Matrix dummy(src.mData, src.mNumCols, src.mNumRows);

    for(int i = 0; i < dummy.mNumRows; ++i)
        for(int j = 0; j < dummy.mNumCols; ++j)
            dummy.mData[i][j] = src.mData[j][i];

    return dummy;
}

double** Matrix::allocateMem_(size_t nrow, size_t ncol)
{
    double **ptr = nullptr;
    try
    {
        if(nrow == 0 || ncol == 0)
            throw std::range_error("Error: There is not existed no row or column in matrix");
        else if(nrow != 0 && ncol != 0)
        {
            ptr = new double*[nrow];

            for(size_t irow = 0; irow < nrow; ++irow)
                ptr[irow] = new double[ncol];

        }
    }
    catch(std::range_error& re)
    {
        std::cout << re.what() << '\n';
        exit(EXIT_FAILURE);
    }
    return ptr;
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double> > init)
    : mData(nullptr)
{
    try
    {
        // The iter pointer point to each comma in the initializer list(if 3x3 there are only 2 commas, and 3 child init_list in parent init_list)
        // This iter -> each row via comma
        for(auto iter = init.begin(); iter < init.end() - 1; ++iter)
        {
            const size_t mSize = iter->size();

            if(mSize == (iter+1)->size()) continue;
            else
            {
                throw std::length_error("Your input value is not a matrix!");
            }
        }
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        exit(EXIT_FAILURE);
    }

    mNumRows = init.size();
    mNumCols = init.begin()->size(); // use begin() as a pointer point to size() method

    this->mData = new double*[mNumRows];
    for(int i = 0; i < mNumRows; ++i)
        this->mData[i] = new double[mNumCols];

    int i = 0, j = 0;
    // Traverse through each row of the parent initializer list
    for(auto iter = init.begin(); iter < init.end(); ++iter, ++i)
    {
        // Traverse through each col like each member in the child initializer list
        for(auto _iter = iter->begin();_iter < iter->end(); ++_iter, ++j)
        {
            this->mData[i][j] = *_iter;
        }
        j = 0;
    }
}

Matrix::Matrix(const Matrix& m) noexcept
    : mNumRows( m.mNumRows )
    , mNumCols( m.mNumCols )
    , mData( allocateMem_(m.mNumRows, m.mNumCols) )
{
    m.display();
    if(mData)
        for(int i = 0; i < mNumRows; ++i)
            for(int j = 0; j < mNumCols; ++j)
                mData[i][j] = m.mData[i][j];
}

Matrix::~Matrix()
{
    // (new [] then use delete [] ) or (new then use delete) to avoid "undefined behavior"
    if(mData)
    {
        for(int irow = 0; irow < mNumRows; ++irow)
            delete [] mData[irow];
    delete [] mData;
    mData = nullptr;
    }
}

Matrix& Matrix::operator=(const Matrix& m) noexcept
{
    if(this != &m)
    {
        if(mNumRows != m.mNumRows || mNumCols != m.mNumCols)
        {
            for(int i = 0; i < mNumRows; ++i) delete [] mData[i];
            delete [] mData;
            mNumCols = m.mNumCols;
            mNumRows = m.mNumRows;
            mData = allocateMem_(mNumRows, mNumCols);
        }
    }
    if(m.mData)
        for(int i = 0; i < mNumRows; ++i)
            for(int j = 0; j < mNumCols; ++j)
                mData[i][j] = m.mData[i][j];

    return (*this);
}

Matrix::Matrix(Matrix&& m) noexcept
    : mNumRows(0)
    , mNumCols(0)
    , mData(nullptr)
{
    mNumCols = m.mNumCols;
    mNumRows = m.mNumRows;
    mData = m.mData;

    m.mNumCols = m.mNumRows = 0;

    m.mData = nullptr;
}

Matrix& Matrix::operator=(Matrix&& m) noexcept
{

    for(int i = 0; i < mNumRows; ++i){
        delete [] mData[i];
    }
    delete [] mData;
    mNumRows = m.mNumRows;
    mNumCols = m.mNumCols;
    mData = m.mData;

    m.mData = nullptr;
    m.mNumCols = 0;
    m.mNumRows = 0;

    return (*this);
}

std::ostream& operator<<(std::ostream& os, const Matrix& obj)
{
    for(int i = 0; i < obj.mNumRows; ++i)
    {
        for(int j = 0; j < obj.mNumCols; ++j)
        {
            os << obj.mData[i][j] << "\t";
        }
        os << '\n';
    }
    return os;
}

std::istream& operator>>(std::istream& is, Matrix& obj)
{
    for(int i = 0; i < obj.mNumRows; ++i)
    {
        for(int j = 0; j < obj.mNumCols; ++j)
        {
            is >> obj.mData[i][j];
        }
    }
    return is;
}

double& Matrix::operator()(int lhs, int rhs)
{
    try
    {
        if(lhs < 0 || rhs < 0 || lhs > mNumRows || rhs > mNumCols)
            throw std::out_of_range("Error: Out of range!");
    }
    catch(const std::out_of_range& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }
    return mData[lhs-1][rhs-1];
}

double& Matrix::operator()(int lhs, int rhs) const
{
    try
    {
        if(lhs > mNumRows || rhs > mNumCols)
            throw std::out_of_range("Error: Out of range!");
    }
    catch(std::out_of_range& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }
    return mData[lhs - 1][rhs - 1];
}

bool Matrix::operator[](const Matrix& rhs)
{
    try
    {
        if(rhs.mNumCols < 0 || rhs.mNumRows < 0 || rhs.mNumRows > mNumRows || rhs.mNumCols > mNumCols)
            throw std::logic_error("Error: There is not existed the matrix that appropriate the matrix you index!");
    }
    catch(std::logic_error& le)
    {
        std::cout << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

bool Matrix::operator[](const Matrix& rhs) const
{
    try
    {
        if(rhs.mNumCols < 0 || rhs.mNumRows < 0 || rhs.mNumRows > mNumRows || rhs.mNumCols > mNumCols)
            throw std::logic_error("Error: There is not existed the matrix that appropriate the matrix you index!");
    }
    catch(std::logic_error& le)
    {
        std::cout << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

const Matrix Matrix::operator/(const double& rhs)
{
    Matrix res(mNumRows, mNumCols);

    for(int i = 0; i < res.mNumRows; ++i){
        for(int j = 0; j < res.mNumCols; ++j){
            res.mData[i][j] = mData[i][j] / rhs;
        }
    }
    return res;
}

const Matrix Matrix::operator*(const double& rhs)
{
    Matrix res(mNumRows, mNumCols);

    for(int i = 0; i < res.mNumRows; ++i)
    {
        for(int j = 0; j < res.mNumCols; ++j)
        {
            res.mData[i][j] = mData[i][j] * rhs;
        }
    }
    return res;
}

const Matrix Matrix::operator*(const double& rhs) const
{
    Matrix res(mNumRows, mNumCols);

    for(int i = 0; i < res.mNumRows; ++i)
    {
        for(int j = 0; j < res.mNumCols; ++j)
        {
            res.mData[i][j] = mData[i][j] * rhs;
        }
    }
    return res;
}

Matrix operator*(const double& lhs, const Matrix& rhs)
{
    Matrix res(rhs.mNumRows, rhs.mNumCols);

    for(int i = 0; i < res.mNumRows; ++i)
    {
        for(int j = 0; j < res.mNumCols; ++j)
        {
            res.mData[i][j] = lhs * rhs.mData[i][j];
        }
    }
    return res;
}

Vector operator*(const Matrix& mat, const Vector& vec)
{
    //static_assert(mat.mNumCols != vec.mSize, "Matrix object and Vector object is not the same size!");
    assert(mat.mNumCols == static_cast<int>(vec.mSize));

    Vector sol(mat.mNumRows);   // since mat(row*col) * vec(col*1) => sol(row*1);

    for(int i = 0; i < mat.mNumRows; ++i)
        for(int j = 0; j < mat.mNumCols; ++j)
            sol.mData[i] += mat.mData[i][j] * vec.mData[j];

    return sol;
}

const Matrix Matrix::operator+(const Matrix& rhs) const
{
    if(mNumRows != rhs.mNumRows || mNumCols != rhs.mNumCols)
    {
        std::cout << "Improper dimensions" << '\n';
        exit(EXIT_FAILURE);
    }
    Matrix res(mNumRows, mNumCols);

    for(int i = 0; i < res.mNumRows; ++i)
    {
        for(int j = 0; j < res.mNumCols; ++j)
        {
            res.mData[i][j] = mData[i][j] + rhs.mData[i][j];
        }
    }
    return res;
}

const Matrix Matrix::operator+(const Matrix& rhs)
{
    if(mNumRows != rhs.mNumRows || mNumCols != rhs.mNumCols)
    {
        std::cout << "Improper dimensions" << '\n';
        exit(EXIT_FAILURE);
    }
    Matrix res(mNumRows, mNumCols);

    for(int i = 0; i < res.mNumRows; ++i)
    {
        for(int j = 0; j < res.mNumCols; ++j)
        {
            res.mData[i][j] = mData[i][j] + rhs.mData[i][j];
        }
    }
    return res;
}

const Matrix Matrix::operator-(const Matrix& rhs) const
{
    if(mNumRows != rhs.mNumRows || mNumCols != rhs.mNumCols)
    {
        std::cout << "Improper dimensions" << '\n';
        exit(EXIT_FAILURE);
    }
    Matrix res(mNumRows, mNumCols);

    for(int i = 0; i < res.mNumRows; ++i)
        for(int j = 0; j < res.mNumCols; ++j)
            res.mData[i][j] = mData[i][j] - rhs.mData[i][j];

    return res;
}

const Matrix Matrix::operator-(const Matrix& rhs)
{
    if(mNumRows != rhs.mNumRows || mNumCols != rhs.mNumCols)
    {
        std::cout << "Improper dimensions" << '\n';
        exit(EXIT_FAILURE);
    }
    Matrix res(mNumRows, mNumCols);

    for(int i = 0; i < res.mNumRows; ++i)
        for(int j = 0; j < res.mNumCols; ++j)
            res.mData[i][j] = mData[i][j] - rhs.mData[i][j];

    return res;
}

const Matrix Matrix::operator*(const Matrix& rhs) const
{
    if(rhs.mNumRows == 1 && rhs.mNumCols == 1)
    {
        Matrix r(mNumRows, mNumCols);

        for(int i = 0; i < mNumRows; ++i)
            for(int j = 0; j < mNumCols; ++j)
                r.mData[i][j] = mData[i][j] * rhs.mData[0][0];
        return r;
    }

    if(mNumRows == 1 && mNumCols == 1)
    {
        Matrix r(rhs.mNumRows, rhs.mNumCols);

        for(int i = 0; i < rhs.mNumRows; ++i)
            for(int j = 0; j < rhs.mNumCols; ++j)
                r.mData[i][j] = rhs.mData[i][j] * mData[0][0];
        return r;
    }

    try
    {
        if(this->mNumCols != rhs.mNumRows)
            throw std::logic_error("Error: Column of the first matrix must be equal to row of second matrix");
        else if(this->mNumCols == rhs.mNumRows)
        {
            Matrix res(rhs.mNumRows, rhs.mNumCols);
            for(int i = 0; i < this->mNumRows; ++i)
            {
                for(int j = 0; j < rhs.mNumCols; ++j)
                {
                    double sum = 0;
                    for(int k = 0; k < this->mNumCols; ++k)
                        sum += this->mData[i][k] * rhs.mData[k][j];
                    res.mData[i][j] = sum;
                }
            }
            return res;
        }
        else
        {
            throw std::domain_error("Error: Column of the first matrix must be equal to row of second matrix");
        }
    }
    catch(const std::domain_error& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }
    catch(const std::logic_error& le)
    {
        std::cout << le.what() << '\n';
        exit(EXIT_FAILURE);
    }

}

Matrix Matrix::operator*(const Matrix& rhs)
{
    if(rhs.mNumRows == 1 && rhs.mNumCols == 1)
    {
        Matrix r(mNumRows, mNumCols);

        for(int i = 0; i < mNumRows; ++i)
            for(int j = 0; j < mNumCols; ++j)
                r.mData[i][j] = mData[i][j] * rhs.mData[0][0];
        return r;
    }

    if(mNumRows == 1 && mNumCols == 1)
    {
        Matrix r(rhs.mNumRows, rhs.mNumCols);

        for(int i = 0; i < rhs.mNumRows; ++i)
            for(int j = 0; j < rhs.mNumCols; ++j)
                r.mData[i][j] = rhs.mData[i][j] * mData[0][0];
        return r;
    }

    try
    {
        if(this->mNumCols != rhs.mNumRows)
            throw std::logic_error("Error: Column of the first matrix must be equal to row of second matrix-1st argument");
        else if(this->mNumCols == rhs.mNumRows)
        {
            Matrix res(mNumRows, rhs.mNumCols);
            for(int i = 0; i < this->mNumRows; ++i)
            {
                for(int j = 0; j < rhs.mNumCols; ++j)
                {
                    double sum = 0.0;
                    for(int k = 0; k < this->mNumCols; ++k)
                        sum += this->mData[i][k] * rhs.mData[k][j];
                    res.mData[i][j] = sum;
                }
            }
            return res;
        }
        else
        {
            throw std::domain_error("Error: Column of the first matrix must be equal to row of second matrix-2nd argument");
        }
    }
    catch(const std::domain_error& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }
    catch(const std::logic_error& le)
    {
        std::cout << le.what() << '\n';
        exit(EXIT_FAILURE);
    }

}

const Matrix& Matrix::operator+=(const Matrix& rhs)
{
    try
    {
        if(!((mNumRows == rhs.mNumRows) || (mNumCols == rhs.mNumCols)))
            throw std::logic_error("Error: there is not available an addition between 2 different size matrices");
        else
        {
            for(int i = 0; i < mNumRows; ++i)
                for(int j = 0; j < mNumCols; ++j)
                    mData[i][j] += rhs.mData[i][j];
        }
    }
    catch(std::logic_error& le)
    {
        std::cout << le.what() << '\n';
        exit(EXIT_FAILURE);
    }
    return (*this);
}

const Matrix& Matrix::operator-=(const Matrix& rhs)
{
    try
    {
        if(!((mNumRows == rhs.mNumRows) || (mNumCols == rhs.mNumCols)))
            throw std::logic_error("Error: there is not available an addition between 2 different size matrices");
        else
        {
            for(int i = 0; i < mNumRows; ++i)
                for(int j = 0; j < mNumCols; ++j)
                    mData[i][j] -= rhs.mData[i][j];
        }
    }
    catch(std::logic_error& le)
    {
        std::cout << le.what() << '\n';
        exit(EXIT_FAILURE);
    }
    return (*this);
}

const Matrix& Matrix::operator*=(const Matrix& rhs)
{
    //Matrix res(mNumRows, rhs.mNumCols);

    try
    {
        if(this->mNumCols != rhs.mNumRows)
            throw std::logic_error("Error: Column of the first matrix must be equal to row of second matrix");
        else if(this->mNumCols == rhs.mNumRows)
        {
            for(int i = 0; i < this->mNumRows; ++i)
            {
                for(int j = 0; j < rhs.mNumCols; ++j)
                {
                    double sum = 0;
                    for(int k = 0; k < this->mNumCols; ++k)
                        sum += this->mData[i][k] * rhs.mData[k][j];
                    this->mData[i][j] = sum;
                }
            }
        }
        else
        {// Never happen because situations is the same with the logic_error, but using this to understand the order of try catch throw
            throw std::out_of_range("Error: Column of the first matrix must be equal to row of second matrix");
        }
    }
    catch(const std::out_of_range& er)
    {
        std::cerr << er.what() << '\n';
        exit(EXIT_FAILURE);
    }
    catch(const std::logic_error& le)
    {
        std::cerr << le.what() << '\n';
        exit(EXIT_FAILURE);
    }
    return (*this);
}

const Matrix Matrix::operator*=(const double& rhs)
{
    Matrix res(mNumRows, mNumCols);

    for(int i = 0; i < res.mNumRows; ++i)
    {
        for(int j = 0; j < res.mNumCols; ++j)
        {
            res.mData[i][j] = mData[i][j] * rhs;
        }
    }
    return res;
}

const Matrix Matrix::operator*=(const double& rhs) const
{
    Matrix res(mNumRows, mNumCols);

    for(int i = 0; i < res.mNumRows; ++i)
    {
        for(int j = 0; j < res.mNumCols; ++j)
        {
            res.mData[i][j] = mData[i][j] * rhs;
        }
    }
    return res;
}

bool Matrix::operator==(const Matrix& rhs)
{
    try
    {
        if(!((this->mNumCols == rhs.mNumCols) || (this->mNumRows == rhs.mNumRows)))
            throw std::length_error("Error: There is not suitable for comparing 2 different size between 2 matrices");
        else
        {
            for(int i = 0; i < this->mNumRows; ++i)
                for(int j = 0; j < this->mNumCols; ++j)
                {
                    if(this->mData[i][j] == rhs.mData[i][j]) continue;
                    else return false;
                }
        }
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

bool Matrix::operator==(const Matrix& rhs) const
{
    try
    {
        if(!((this->mNumCols == rhs.mNumCols) || (this->mNumRows == rhs.mNumRows)))
            throw std::length_error("Error: There is not suitable for comparing 2 different size between 2 matrices");
        else
        {
            for(int i = 0; i < this->mNumRows; ++i)
                for(int j = 0; j < this->mNumCols; ++j)
                {
                    if(this->mData[i][j] == rhs.mData[i][j]) continue;
                    else return false;
                }
        }
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

bool Matrix::operator!=(const Matrix& rhs)
{
    try
    {
        if(!((this->mNumCols == rhs.mNumCols) || (this->mNumRows == rhs.mNumRows)))
            throw std::length_error("Error: There is not suitable for comparing 2 different size between 2 matrices");
        else
        {
            for(int i = 0; i < this->mNumRows; ++i)
                for(int j = 0; j < this->mNumCols; ++j)
                {
                    if(this->mData[i][j] != rhs.mData[i][j]) continue;
                    else return false;
                }
        }
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

bool Matrix::operator!=(const Matrix& rhs) const
{
    try
    {
        if(!((this->mNumCols == rhs.mNumCols) || (this->mNumRows == rhs.mNumRows)))
            throw std::length_error("Error: There is not suitable for comparing 2 different size between 2 matrices");
        else
        {
            for(int i = 0; i < this->mNumRows; ++i)
                for(int j = 0; j < this->mNumCols; ++j)
                {
                    if(this->mData[i][j] != rhs.mData[i][j]) continue;
                    else return false;
                }
        }
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

bool Matrix::operator<=(const Matrix& rhs)
{
    try
    {
        if(!((this->mNumCols == rhs.mNumCols) || (this->mNumRows == rhs.mNumRows)))
            throw std::length_error("Error: There is not suitable for comparing 2 different size between 2 matrices");
        else
        {
            for(int i = 0; i < this->mNumRows; ++i)
                for(int j = 0; j < this->mNumCols; ++j)
                {
                    if(this->mData[i][j] <= rhs.mData[i][j]) continue;
                    else return false;
                }
        }
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

bool Matrix::operator<=(const Matrix& rhs) const
{
    try
    {
        if(!((this->mNumCols == rhs.mNumCols) || (this->mNumRows == rhs.mNumRows)))
            throw std::length_error("Error: There is not suitable for comparing 2 different size between 2 matrices");
        else
        {
            for(int i = 0; i < this->mNumRows; ++i)
                for(int j = 0; j < this->mNumCols; ++j)
                {
                    if(this->mData[i][j] <= rhs.mData[i][j]) continue;
                    else return false;
                }
        }
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

bool Matrix::operator>=(const Matrix& rhs)
{
    try
    {
        if(!((this->mNumCols == rhs.mNumCols) || (this->mNumRows == rhs.mNumRows)))
            throw std::length_error("Error: There is not suitable for comparing 2 different size between 2 matrices");
        else
        {
            for(int i = 0; i < this->mNumRows; ++i)
                for(int j = 0; j < this->mNumCols; ++j)
                {
                    if(this->mData[i][j] >= rhs.mData[i][j]) continue;
                    else return false;
                }
        }
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

bool Matrix::operator>(const Matrix& rhs)
{
    try
    {
        if(!((this->mNumCols == rhs.mNumCols) || (this->mNumRows == rhs.mNumRows)))
            throw std::length_error("Error: There is not suitable for comparing 2 different size between 2 matrices");
        else
        {
            for(int i = 0; i < this->mNumRows; ++i)
                for(int j = 0; j < this->mNumCols; ++j)
                {
                    if(this->mData[i][j] > rhs.mData[i][j]) continue;
                    else return false;
                }
        }
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

bool Matrix::operator>(const Matrix& rhs) const
{
    try
    {
        if(!((this->mNumCols == rhs.mNumCols) || (this->mNumRows == rhs.mNumRows)))
            throw std::length_error("Error: There is not suitable for comparing 2 different size between 2 matrices");
        else
        {
            for(int i = 0; i < this->mNumRows; ++i)
                for(int j = 0; j < this->mNumCols; ++j)
                {
                    if(this->mData[i][j] > rhs.mData[i][j]) continue;
                    else return false;
                }
        }
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

bool Matrix::operator<(const Matrix& rhs)
{
    try
    {
        if(!((this->mNumCols == rhs.mNumCols) || (this->mNumRows == rhs.mNumRows)))
            throw std::length_error("Error: There is not suitable for comparing 2 different size between 2 matrices");
        else
        {
            for(int i = 0; i < this->mNumRows; ++i)
                for(int j = 0; j < this->mNumCols; ++j)
                {
                    if(this->mData[i][j] < rhs.mData[i][j]) continue;
                    else return false;
                }
        }
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

bool Matrix::operator<(const Matrix& rhs) const
{
    try
    {
        if(!((this->mNumCols == rhs.mNumCols) || (this->mNumRows == rhs.mNumRows)))
            throw std::length_error("Error: There is not suitable for comparing 2 different size between 2 matrices");
        else
        {
            for(int i = 0; i < this->mNumRows; ++i)
                for(int j = 0; j < this->mNumCols; ++j)
                {
                    if(this->mData[i][j] < rhs.mData[i][j]) continue;
                    else return false;
                }
        }
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }
    return true;
}

Matrix Matrix::operator^(const char& T)
{
    try
    {
        if(T != 'T') throw std::logic_error("Wrong purpose of ^T operator!");
        else{
            Matrix mat(mNumCols, mNumRows);    // Size of transpose matrix is reversed
            for(int i = 0; i < mat.mNumRows; i++)
                for(int j = 0; j < mat.mNumCols; j++)
                    mat.mData[i][j] = mData[j][i];

        return mat;
        }
    }
    catch(const std::logic_error& le)
    {
        std::cerr << "Please check it again" << '\n';
        exit(EXIT_FAILURE);
    }
}

Matrix Matrix::operator^(const int& _inverseValue)
{
    Matrix inv(mNumRows, mNumCols);
    try
    {
        if(_inverseValue != -1) throw "Wrong purpose of operator ^-1!";
        else{
            matrixInverse(*this, inv);
        }
    }
    catch(const std::basic_string<char> msg)
    {
        std::cout << msg << '\n';
        exit(EXIT_FAILURE);
    }
    catch(const std::logic_error& le)
    {
        return Matrix(0,0);
    }
    return inv;
}

void Matrix::deallocateMem(double** mat, int n)
{
    for(int i = 0; i < n; ++i) delete [] mat[i];
    delete [] mat;
}

void Matrix::getCofactor(const Matrix& src, Matrix& dest)
{
    // i -> row, j -> col
    int colCount = 0, rowCount = 0;

    // Looping for each element of the matrix
    for(int i = 0; i < src.mNumRows; ++i)
    {
        for(int j = 0; j < src.mNumCols; ++j)
        {
            // Copying into temporary matrix only those element which are not in given
            // row and column
            // De Morgan Rule
            //if( ( (i != row) && (j != col)) )
            if( !((i == src.mNumRows) || (j == src.mNumCols)) )
            {
                dest.mData[rowCount][colCount++] = src.mData[i][j];

                // Row is filled, so increment row index and reset col index
                if(colCount == src.mNumRows - 1)
                {
                    colCount = 0;
                    rowCount++;
                }
            }
        }
    }
}

void Matrix::getCofactor(Matrix& temp, int start, int idx)
{
    // i -> row, j -> col
    int colCount = 0, rowCount = 0;

    // Looping for each element of the matrix
    for(int i = 0; i < this->mNumRows; ++i)
    {
        for(int j = 0; j < this->mNumCols; ++j)
        {
            // Copying into temporary matrix only those element which are not in given
            // row and column
            // De Morgan Rule
            //if( ( (i != row) && (j != col)) )
            if( !((i == start) || (j == idx)) )
            {
                temp.mData[rowCount][colCount++] = this->mData[i][j];

                // Row is filled, so increment row index and reset col index
                if(colCount == this->mNumRows - 1)
                {
                    colCount = 0;
                    rowCount++;
                }
            }
        }
    }
    std::cout << "Matrix cofactor: " << '\n'
              << temp << '\n';
}

/*
void Matrix::getCofactor(const int& start, const int& idx, const int& n)
{
    Matrix dest(mData, mNumRows, mNumCols);

    getCofactor((*this), dest, start, idx, n);
}
*/

void Matrix::getCofactor(const Matrix& src, Matrix& temp, int start, int idx, int n)
{
    int i = 0, j = 0;

    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element which are not in given row and column
            if (row != start && col != idx)
            {
                temp.mData[i][j++] = src.mData[row][col];

                // Row is filled, so increase row index and reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

// do not use this method
double Matrix::determinantMatix(const Matrix& mat, int n)
{
    double det = 0; // Initialize result
    if(mat.mNumCols == mat.mNumRows)
    {
        //  Base case : if matrix contains single element
        if (n == 1)
            return mat.mData[0][0];
        else if(n == 2)
        {
            det = mat.mData[0][0] * mat.mData[1][1] - mat.mData[1][0] * mat.mData[0][1];
        }
        else
        {
            Matrix temp(mat.mNumRows, mat.mNumCols);

            int sign = 1; // To store sign multiplier

            // Iterate for each element of first row
            for (int i = 0; i < n; i++)
            {
                // Getting Cofactor of mat[0][f]
                getCofactor(mat, temp, 0, i, n);

                det += sign * mat.mData[0][i] * determinantMatix(temp, n - 1);

                // terms are to be added with alternate sign
                sign *= -sign;
            }
        }
    }
    else
    {
        std::cout << "When matrix is determinant, matrix must be a square matrix" << '\n';
    }
    //std::cout << "Determinant of matrix: " << '\n' << det << '\n';
    return det;

}

double Matrix::determinantMatix()
{
    Matrix orgMatrix(this->mData, this->mNumRows, this->mNumCols);
    return determinantMatix(orgMatrix);
}

double Matrix::determinantMatix(const Matrix& org_matrix)
{
    if(org_matrix.mNumRows != org_matrix.mNumCols)
        throw std::logic_error("Non-square matrix!");

    // 1 x 1 case
    if(org_matrix.mNumRows == 1)
        return org_matrix.mData[0][0];

    // Variable to hold the result
    double det = 0;
    // Variable to hold the current sign of the entry
    int sign = 1;
    // Create sub_matrix to hold cofactor
    Matrix sub_matrix(org_matrix.mNumRows - 1, org_matrix.mNumCols - 1);

    // Calculate determinant by using recursive function
    for(int i = 0; i < org_matrix.mNumRows; i++)
    {
        // If current element is 0, just skip it
        if(org_matrix.mData[0][i] != 0)
        {
            /* Assign cofactor to sub_matrix */
            // Make a matrix which excludes the current row and column
            int a = 0, b = 0;
            for(int row = 0; row < org_matrix.mNumRows; row++)
            {
                for(int col = 0; col < org_matrix.mNumCols; col++)
                {
                    if(row != 0 && col != i)
                    {
                        sub_matrix.mData[a][b++] = org_matrix.mData[row][col];
                    }

                    if(b == sub_matrix.mNumCols)
                    {
                        b = 0;
                        a++;
                    }
                }
            }

            /* Determinant formula */
            // Recursively computation of determinant
            det += sign * org_matrix.mData[0][i] * determinantMatix(sub_matrix);
        }
        // Change sign for every entry calculated
        sign *= -1;
    }

    return det;
}

void Matrix::adjointMatrix(const Matrix& src, Matrix& adj)
{
    adj = Matrix(src.mNumRows, src.mNumCols);
    if (src.mNumCols == 1)
    {
        adj.mData[0][0] = 1;
        return;
    }
    //std::cout << "\nOutput for source matrix: " << '\n' << src << '\n';
    int sign = 1;
    // temp is used to store cofactors of A[][]

    for(int row = 0; row < src.mNumRows; ++row){
        for(int col = 0; col < src.mNumCols; ++col){
            //adj.getCofactor(src, temp, row, col, src.mNumCols);

            Matrix temp(src.mNumRows-1, src.mNumCols-1);
            int a = 0, b = 0;
            for(int i = 0; i < src.mNumRows; i++)
            {
                for(int j = 0; j < src.mNumCols; j++)
                {
                    if(i != row && j != col)
                    {
                        temp.mData[a][b] = src.mData[i][j];
                        b++;
                    }

                    if(b == temp.mNumCols)
                    {
                        b = 0;
                        a++;
                    }
                }
            }
            //std::cout << "Output for temp matrix: " << '\n' << temp << '\n';
            sign = ( (row + col)%2 == 0 ) ? 1 : -1;
            adj.mData[row][col] = sign * determinantMatix(temp);
        }
    }
    adj = adj.transposeMatrix();
    //std::cout << "Output for adjoint matrix: " << '\n' << adj << '\n';
}

void Matrix::adjointMatrix(Matrix& adj)
{

    if(this->mNumCols == 1)
    {
        adj.mData[0][0] = 1;
        return;
    }
    int sign = 1;

    Matrix temp(mNumRows, mNumCols);

    for(int i = 0; i < mNumRows; ++i)
    {
        for(int j = 0; j < mNumCols; ++j)
        {
            this->getCofactor((*this), temp, i, j, mNumCols);

            // sign of adjoint positive if sum of row and column indexes is even
            sign = ( (i + j)%2 == 0) ? 1 : -1;

            // Interchanging rows and columns to get the transpose of the cofactor matrix
            adj.mData[j][i] = (sign)*(determinantMatix(temp, this->mNumRows - 1));
        }
    }
}

void Matrix::adjointMatrix()
{
    Matrix adj(mNumRows, mNumCols);
    adjointMatrix((*this), adj);

    adj.display();
}

void Matrix::matrixInverse(Matrix& src_mat, Matrix& inv)
{
    double det = determinantMatix(src_mat);
    if(det == 0)
    {
        std::cerr << "Singular/Non-invertible matrix, can not find its inverse";
    }
    if(src_mat.mNumRows == src_mat.mNumCols && src_mat.mNumRows == 1)
    {
        inv.mData[0][0] = 1/src_mat.mData[0][0];
    }
    else if(src_mat.mNumRows == src_mat.mNumCols && src_mat.mNumRows == 2)
    {
        inv.mData[0][0] = src_mat.mData[1][1];
        inv.mData[0][1] = -src_mat.mData[0][1];
        inv.mData[1][0] = -src_mat.mData[1][0];
        inv.mData[1][1] = src_mat.mData[0][0];

        inv = inv * (1.0 / det);
    }
    else
    {
        Matrix adj(src_mat.mNumRows, src_mat.mNumCols);
        src_mat.adjointMatrix(src_mat, adj);
        //std::cout << "Output for adjoint matrix in inverse method: " << '\n' << adj << '\n';
        //std::cout << "Output for det: " << det << '\n';

        for(int i = 0; i < src_mat.mNumRows; ++i)
            for(int j = 0; j < src_mat.mNumCols; ++j)
                inv.mData[i][j] = adj.mData[i][j] / det;
    }
    //std::cout << "Output for inv matrix in inverse method: " << '\n' << inv << '\n';
    //std::cout << "Check inverse matrix: " << '\n' << inv * src_mat << '\n';
    //std::cout << "Here";
}

Matrix Matrix::matrixInverse()
{
    Matrix inv(mNumRows, mNumCols);
    matrixInverse(*this, inv);
    return inv;
}

void Matrix::transposeMatrix(const Matrix& src, Matrix& trans_T)
{
    for(int i = 0; i < trans_T.mNumRows; ++i)
        for(int j = 0; j < trans_T.mNumCols; ++j)
            trans_T.mData[i][j] = src.mData[j][i];
}

Matrix Matrix::transposeMatrix()
{
    Matrix T_trans(mNumCols, mNumRows);

    for(int i = 0; i < T_trans.mNumRows; ++i)
        for(int j = 0; j < T_trans.mNumCols; ++j)
            T_trans.mData[i][j] = mData[j][i];
    return T_trans;
}

bool Matrix::isInvertible(Matrix& mat)
{
    if( (determinantMatix(mat, mat.mNumCols) < 1e-9) && (determinantMatix(mat, mat.mNumCols) > -1e-9) ) return false;
    else return true;
}

bool Matrix::isInvertible()
{
    if( (determinantMatix((*this), this->mNumRows) < 1e-9) && (determinantMatix((*this), this->mNumCols) > -1e9) ) return false;
    else return true;
}

void Matrix::swapRowMatrix(const Matrix& mat, int row1, int row2, int rank_col)
{
    for(int i = 0; i < rank_col; ++i)
    {
        double temp = mat.mData[row1][i];
        mat.mData[row1][i] = mat.mData[row2][i];
        mat.mData[row2][i] = temp;
    }
}

void Matrix::swapRowMatrix(int row1, int row2)
{
    for(int i = 0; i < this->mNumCols; ++i)
    {
        double temp = this->mData[row1][i];
        this->mData[row1][i] = this->mData[row2][i];
        this->mData[row2][i] = temp;
    }
}

void Matrix::swapRowMatrix(int row1, int row2, int rank_col)
{
    for(int i = 0; i < rank_col; ++i)
    {
        double temp = this->mData[row1][i];
        this->mData[row1][i] = this->mData[row2][i];
        this->mData[row2][i] = temp;
    }
}

void Matrix::swapColMatrix(const Matrix& mat, int col1, int col2)
{
    for(int i = 0; i < mat.mNumRows; ++i)
    {
        double temp = mat.mData[i][col1];
        mat.mData[i][col1] = mat.mData[i][col2];
        mat.mData[i][col2] = temp;
    }
}

void Matrix::swapColMatrix(int col1, int col2)
{
    for(int i = 0; i < this->mNumRows; ++i)
    {
        double temp = this->mData[i][col1];
        this->mData[i][col1] = this->mData[i][col2];
        this->mData[i][col2] = temp;
    }
}

void Matrix::display() const
{
    for(int i = 0; i < mNumRows; ++i)
    {
        for(int j = 0; j < mNumCols; ++j)
        {
            std::cout << std::setw(10) << this->mData[i][j] << "\t";
        }
        std::cout << '\n';
    }
}

void Matrix::display(const Matrix& src)
{
    for(int i = 0; i < src.mNumRows; ++i)
    {
        for(int j = 0; j < src.mNumCols; ++j)
        {
            std::cout << src.mData[i][j] << "\t";
        }
        std::cout << '\n';
    }
}

int Matrix::rankMatrix(const Matrix& mat)
{
    int Rank = mat.mNumCols;

    for(int row = 0; row < Rank; ++row)
    {
        // Set all = 0 except diagonal
        if(mat.mData[row][row])
        {
            for(int col = 0; col < mat.mNumRows; ++col)
            {
                if(col != row)
                {
                    // this is algorithm to show how to do
                    double temp = std::fmod(mat.mData[col][row], mat.mData[row][row]);

                    for(int i = 0; i < Rank; ++i)
                        mat.mData[col][i] -= temp * mat.mData[row][i];
                }
            }
        }
        // there 2 cases:
        // If there is a row below diagonal with non-zero entry, then swap this row with that row and process that row
        // If all elements in current column below mData[i][row] = 0, then remove this column by swapping it with last col and reducing col index by 1
        else
        {
            bool reduce = true;

            // find non-zero element in current column
            for(int i = row + 1; i < mat.mNumRows; ++i)
            {
                //swap row contains non-zero element with this row
                if(mat.mData[i][row])
                {
                    swapRowMatrix(mat, row, i, Rank);
                    reduce = false;
                    break;
                }
            }

            // If not found row with non-zero element in current column, then all values in this column = 0
            if(reduce)
            {
                // Reduce number of columns
                Rank--;

                // Copy last column here
                for(int i = 0; i < mat.mNumRows; ++i)
                    mat.mData[i][row] = mat.mData[i][Rank];
            }
            row--;
        }
        mat.display(mat);
        std::cout << '\n';
    }
    return Rank;
}

int Matrix::rankMatrix()
{
    int Rank = this->mNumCols;

    for(int row = 0; row < Rank; ++row)
    {
        if(this->mData[row][row])
        {
            for(int col = 0; col < this->mNumRows; ++col)
            {
                if(!(col == row))
                {
                    double mult = std::fmod(mData[col][row], mData[row][row]);

                    for(int i = 0; i < Rank; ++i)
                        mData[col][i] -= mult * mData[row][i];
                }
            }
        }
        else
        {
            bool reduce = true; // turn flag on

            for(int i = row + 1; i < this->mNumRows; ++i)
            {
                if(mData[i][row])
                {
                    this->swapRowMatrix(row, i, Rank);
                    reduce = false; //turn flag off
                    break; // out of for loop
                }
            }
            if(reduce)
            {
                Rank--;

                for(int i = 0; i < this->mNumRows; ++i)
                    mData[i][row] = mData[i][Rank];
            }
            row--; // decrement to process again
        }
        this->display();
        std::cout << '\n';
    }
    return Rank;
}

void Matrix::matrixMultiply(const Matrix& mat1, const Matrix& mat2, Matrix& matrix_product)
{
    /*
    std::cout << "Output for columns of matrix 1 : " << '\n';
    std::cout << mat1.mNumCols << '\n';
        std::cout << "Output for rows of matrix 2 : " << '\n';
    std::cout << mat2.mNumRows << '\n'; */
    if(mat1.mNumCols != mat2.mNumRows)
    {
        fprintf(stderr, "Wrong dimension for multiplying 2 matrix");
        return;
    }
    for(int i = 0; i < mat1.mNumRows; ++i)
    {
        for(int j = 0; j < mat2.mNumCols; ++j)
        {
            double sum = 0.0;
            for(int k = 0; k < mat1.mNumCols; ++k)
            {
                sum += mat1.mData[i][k] * mat2.mData[k][j];
            }
            matrix_product.mData[i][j] = sum;
        }
    }
}

void Matrix::matrixMultiply(const Matrix& mat)
{
    Matrix matrix_product(mData, mNumRows, mNumCols);

    if(mNumCols != mat.mNumRows)
    {
        fprintf(stderr, "Wrong dimension for multiplying 2 matrix");
        return;
    }
    for(int i = 0; i < mNumRows; ++i)
    {
        for(int j = 0; j < mat.mNumCols; ++j)
        {
            double sum = 0.0;
            for(int k = 0; k < mNumCols; ++k)
            {
                sum += mData[i][k] * mat.mData[k][j];
            }
            matrix_product.mData[i][j] = sum;
        }
    }
}

// Pseudo inverse for A matrix
Matrix Matrix::pseudoInverse(Matrix& srcMat_)
{
    std::cout << "source Matrix: " << '\n'
              << srcMat_ << '\n';
    Matrix matrixT(this->mNumCols, this->mNumRows);
    transposeMatrix(srcMat_, matrixT);
    std::cout << "Transpose Matrix: " << '\n'
              << matrixT << '\n';

    Matrix pinv(this->mNumCols, this->mNumRows);
    //Matrix pinv = pinv_Matrix(srcMat_);
    // Solve for under-determined system
    // If rows of matrix are linearly independent, then right inverse of A
    // A+ = A^T*(A*A^T)^-1, here A+ is a right inverse of A, i.e, A*A+ = E, where E is an identity matrix
    if( (srcMat_.mNumRows < srcMat_.mNumCols) )
    {
        Matrix matrix_product(srcMat_.mData, srcMat_.mNumRows, srcMat_.mNumRows);
        matrix_product = srcMat_ * matrixT;
        //matrixMultiply(srcMat_, matrixT, matrix_product);
        std::cout << "matrix_product: " << '\n'
                  << matrix_product << '\n';

        if(isInvertible(matrix_product))
        {
            Matrix inv(srcMat_.mData, srcMat_.mNumRows, srcMat_.mNumRows);
            matrixInverse(matrix_product, inv);

            pinv = matrixT * inv;
            //matrixMultiply(matrixT, inv, pinv);

            std::cout << "Pseudo inverse of matrix: " << '\n';
            //display(pinv); // with size NxM
        }
    }
    // Solve for over-determined system
    // If columns of matrix A are linearly independent, so A^T*A is invertible then left inverse of A
    // A+ = (A^T*A)^-1*A^T, here A+ is a left inverse of A, i.e, A+*A = E, where E is an identity matrix
    if((srcMat_.mNumRows > srcMat_.mNumCols) && (rankMatrix(srcMat_) == srcMat_.mNumCols))
    {
        Matrix matrix_product(srcMat_.mData, srcMat_.mNumCols, srcMat_.mNumCols);
        matrix_product = matrixT * srcMat_;
        //matrixMultiply(matrixT, srcMat_, matrix_product);

        std::cout << "matrix_product: " << '\n'
                  << matrix_product << '\n';

        std::cout << "invertible matrix: " << '\n'
                  << isInvertible(matrix_product) << '\n';
        if(isInvertible(matrix_product))
        {
            Matrix inv(srcMat_.mData, srcMat_.mNumCols, srcMat_.mNumCols);
            matrixInverse(matrix_product, inv);
            std::cout << "Matrix inverse: " << '\n'
                      << inv << '\n';

            pinv = inv * matrixT;
            //matrixMultiply(inv, matrixT, pinv);

            std::cout << "Pseudo inverse of matrix: " << '\n';
            //display(pinv); // with size NxM
        }
    }
    // If both the columns and the rows of the matrix are linearly independent, then the matrix is invertible and
    // the pseudo inverse is equal to the inverse of the matrix.
    if(srcMat_.mNumCols == srcMat_.mNumRows)
    {

        if(isInvertible(srcMat_))
        {
            matrixInverse(srcMat_, pinv);

            std::cout << "Pseudo inverse of matrix: " << '\n';
            //display(pinv);
        }
    }
    return pinv;
}

Vector Matrix::diag(const Matrix& mat, int k)
{
    Vector A;
    auto min = [](const auto& a, const auto& b)
    {
        return (a < b) ? a : b;
        //return std::min(a, b);
    };
    if(k > 0){
        if(mat.mNumCols - k >= 0){
            A = Vector(min(mat.mNumCols - k, mat.mNumRows));
            for(size_t i = 0; i + k < A.mSize; ++i)
            {
                A.mData[i] = mat.mData[i][i + k];
            }
        }
    }
    if(k < 0){
        int k_abs = std::abs(k);
        if(mat.mNumRows - k_abs >= 0){
            A = Vector(min(mat.mNumRows - k_abs, mat.mNumCols));
            for(size_t i = 0; i + k_abs <  A.mSize; ++i)
            {
                A.mData[i] = mat.mData[i + k_abs][i];
            }
        }
    }
    if(k == 0){
        A = Vector(min(mat.mNumRows, mat.mNumCols));
        for(size_t i = 0; i < A.mSize; ++i)
        {
            A.mData[i] = mat.mData[i][i];
        }
    }

    return A;
}

Matrix Matrix::geninv()
{
    try
    {
        if(!mData){
            throw std::logic_error("Empty matrix!");
        }
    }
    catch(const std::logic_error& le)
    {
        std::cout << "logic error: " << le.what() << '\n';
        return Matrix(0,0);
    }
    //std::cout << "Output for this matrix: " << '\n' << (*this) << '\n';
    Matrix A(0, 0);
    int n = mNumCols, m = mNumRows; // temp object
    // Returns the Moore-Penrose inverse of the argument
    bool transpose = false;
    if(mNumRows < mNumCols){
        transpose = true;
        A = (*this)*(this->transposeMatrix());
        n = m;
    }
    else{
        A = ((*this)^'T')*(*this);
    }

    // Full rank Cholesky factorization of A
    Vector dA = Matrix::diag(A);

    auto minValue = [](const auto& a, const auto& b)
    {
        return (a < b) ? a : b;
        //return std::min(a, b);
    };
    // Capture all the above object and lambda expression, which is called lambda capture
    // Choose elements in the Vector object then return the min element
    auto min = [&](const Vector& vec) -> double
    {
        double minIdx = 0;

        for(std::size_t i = 0; i < vec.mSize - 1; ++i)
        {
            if( minValue(vec.mData[i], vec.mData[i+1]) > 0){
                minIdx = minValue(vec.mData[i], vec.mData[i+1]);
            }
        }
        return minIdx;
    };
    double tol = min(dA)*1e-9;

/*
    double tol = std::abs(dA.mData[0]);
    for(std::size_t i = 1; i < dA.mSize; i++)
    {
        if(tol > std::abs(dA.mData[i]) && dA.mData[i] != 0)
        {
            tol = dA.mData[i];
        }
    }
    tol *= 1e-9;
*/
    auto zeros = [](std::size_t row, std::size_t col)
    {
        Matrix mat(row, col);
        for(int i = 0; i < mat.mNumRows; ++i){
            for(int j = 0; j < mat.mNumCols; ++j){
                mat.mData[i][j] = 0;
            }
        }
        return mat;
    };
    Matrix L = zeros(A.mNumRows, A.mNumCols);

    //auto sqrt = [](const Matrix& mat, int irow, int icol) { return std::sqrt(mat.mData[irow][icol]); };
    auto subMatrix = [&](const Matrix& mat, const std::pair<int, int>& row, const std::pair<int, int>& col)
    {
        Matrix temp(row.second - row.first + 1, col.second - col.first + 1);
        for(int i = row.first; i <= row.second; ++i){
            for(int j = col.first; j <= col.second; ++j){
                temp.mData[i - row.first][j - col.first] = mat.mData[i][j];
            }
        }
        return temp;
    }; // final error ??
    auto iRowMatrix = [&](const Matrix& mat, int irow, const std::pair<int, int>& col)
    {
        Matrix temp(1, col.second - col.first + 1);
        for(int i = col.first; i <= col.second; ++i){
            temp.mData[0][i - col.first] = mat.mData[irow][i];
        }
        return temp;
    }; // main error, step into -> tranposeMatrix() method
    auto iColMatrix = [&](const Matrix& mat, const std::pair<int, int>& row, int icol)
    {
        Matrix temp(row.second - row.first + 1, 1);
        for(int i = row.first; i <= row.second; ++i){
            temp.mData[i - row.first][0] = mat.mData[i][icol];
        }
        return temp;
    };
    int r = 0;
    // Check it out ???
    for(int k = 0; k < n; ++k)
    {
        //std::cout << "#" << k + 1 << " loop(s)" << "\n";
        Matrix vecContainer(0, 0);
        std::pair<int, int> row(k, n - 1);
        if(r == 0){
            vecContainer = iColMatrix(A, row, k);
        }
        else
        {
            std::pair<int, int> col(0, r - 1);
            //iColMatrix(A, row, k).display();
            //subMatrix(L, row, col).display();
            //auto xyz = iRowMatrix(L, k, col);
            //xyz.display();
            //xyz.transposeMatrix().display();
            vecContainer = iColMatrix(A, row, k) - subMatrix(L, row, col) * (iRowMatrix(L, k, col).transposeMatrix());
        }
        //vecContainer.display();
        for(int i = k; i < n; ++i)
        {
            L.mData[i][r] = vecContainer.mData[i - k][0];
        }

        if(L.mData[k][r] > tol){
            L.mData[k][r] = std::sqrt(L.mData[k][r]);

            if(k < n - 1){
                std::pair<int, int> rowtemp(k + 1, n - 1);
                Matrix vecContainer1 = iColMatrix(L, rowtemp, r) / L.mData[k][r];
                for(int i = k + 1; i < n; ++i){
                    L.mData[i][r] = vecContainer1.mData[i - (k + 1)][0]; // error number(garbage L.mData[i][r]);
                }
            }
        }
        else{
            r--;
            if(k == (n - 1)){
                break;
            }
        }
        //std::cout << "Output for L matrix: " << '\n' << L << '\n';
        r++;
    }
    if(r > n - 1){
        r = n - 1;
    }
    //std::cout << "Output for L before submatrix: " << '\n';
    //L.display();
    L = subMatrix(L, std::pair(0, L.mNumRows - 1), std::pair(0, r));
    //std::cout << "Output for L: " << '\n' << L << '\n';
    //std::cout << "Output for L transpose: " << '\n' << (L^'T') << '\n';
    // Computation of the generalized inverse of G
    //std::cout << "L transpose * L: " << '\n';
    //((L^'T') * L).display();
    Matrix M = ( (L^'T') * L)^-1;
    //std::cout << "Output for M matrix with inverse: " << '\n';
    //M.display();

    /*Matrix M = ( ( L^'T' )*L );
    M.matrixInverse();*/
    Matrix Y(0, 0); // geninv function
    if(transpose){
        Y = ((*this)^'T') * L * M * M * (L^'T');
    }
    else{
        Y = L * M * M * (L^'T') * ((*this)^'T');
    }
    return Y;
}

Matrix Matrix::pseudoInverse()
{
    std::cout << "source Matrix: " << '\n'
              << (*this) << '\n';
    Matrix matrixT(this->mNumCols, this->mNumRows);
    transposeMatrix((*this), matrixT);
    std::cout << "Transpose Matrix: " << '\n'
              << matrixT << '\n';

    Matrix pinv(this->mNumCols, this->mNumRows);
    // Solve for under-determined system
    // If rows of matrix are linearly independent, then right inverse of A
    // A+ = A^T*(A*A^T)^-1, here A+ is a right inverse of A, i.e, A*A+ = E, where E is an identity matrix
    if( (mNumRows < mNumCols) )
    {

        // Note: M -> row size of source Matrix , N -> column size of source Matrix

        // (A * A transpose) means new matrix with size MxM
        // source matrix with MxN size and transpose matrix has NxM size
        Matrix matrix_product(mData, mNumRows, mNumRows);
        //matrixMultiply((*this), matrixT, matrix_product);
        matrix_product = (*this) * matrixT;
        std::cout << "matrix_product: " << '\n'
                  << matrix_product << '\n';

        if(isInvertible(matrix_product))
        {
            Matrix inv(this->mNumRows, this->mNumRows);
            matrixInverse(matrix_product, inv);

            std::cout << "matrix inverse: " << '\n'
                      << inv << '\n';

            // And multiply between A transpose with size NxM and new inverse matrix size MxM
            // pinv matrix with size NxM
            pinv = matrixT * inv;
            //matrixMultiply(matrixT, inv, pinv);

            std::cout << "Pseudo inverse of matrix: " << '\n';
            //display(pinv);
            //return pinv;
        }
    }
    // Solve for over-determined system
    // If columns of matrix A are linearly independent, so A^T*A is invertible then left inverse of A
    // A+ = (A^T*A)^-1*A^T, here A+ is a left inverse of A, i.e, A+*A = E, where E is an identity matrix
    if( (mNumRows > mNumCols) && (rankMatrix(*this) == mNumCols) )
    {

        Matrix matrix_product(mData, mNumCols, mNumCols);
        matrix_product = matrixT * (*this);
        //matrixMultiply(matrixT, (*this), matrix_product);

        std::cout << "matrix_product: " << '\n'
                  << matrix_product << '\n';

        std::cout << "invertible matrix: " << '\n'
                  << isInvertible(matrix_product) << '\n';


        if(isInvertible(matrix_product))
        {
            Matrix inv(mData, mNumCols, mNumCols);
            //inv = matrix_product^-1;
            matrixInverse(matrix_product, inv);
            std::cout << "Matrix inverse: " << '\n'
                      << inv << '\n';

            pinv = inv * matrixT;
            //matrixMultiply(inv, matrixT, pinv);

            std::cout << "Pseudo inverse of matrix: " << '\n';
            //display(pinv);
            //return pinv;
        }
    }
    // If both the columns and the rows of the matrix are linearly independent, then the matrix is invertible and
    // the pseudo inverse is equal to the inverse of the matrix.
    if( (mNumCols == mNumRows)  )
    {
        //Matrix pinv(this->mNumRows, this->mNumCols);
        if(isInvertible(*this))
        {
            matrixInverse((*this), pinv);

            std::cout << "Pseudo inverse of matrix: " << '\n';
            //display(pinv);
            //return pinv;
        }
    }
    //std::cout << "Pseudo inverse of matrix: " << '\n';
    return pinv;
}

Matrix Matrix::pseudoInverseTikhonov(const double& lambda)
{
    Matrix identity(mNumCols, mNumCols);
    identity.Identity();

    //  identity.display();
    Matrix Tmatrix(mNumCols, mNumRows);
    Tmatrix.transposeMatrix((*this), Tmatrix);

    Matrix inv(identity.mNumCols, identity.mNumCols);
    Matrix product(mNumCols, mNumCols);
    product = Tmatrix*(*this) + lambda*identity;
    inv.matrixInverse(product, inv);
    Matrix c = inv * Tmatrix;

    return c;
}

/* i-> index for row, j -> index for column */
// Square subMatrix
void Matrix::printSubMatrix(const Matrix& srcMat_, std::size_t block_size, std::size_t i, std::size_t j)
{
    for(std::size_t m = i * block_size; m < block_size * (i + 1); ++m)
    {
        for(std::size_t n = j * block_size; n < block_size * (j + 1); ++n)
        {
            std::cout << srcMat_.mData[m][n] << ' ';
        }
        std::cout << '\n';
    }
}

void Matrix::setSubMatrix(const Matrix& subMat_, std::size_t block_size, std::size_t i, std::size_t j)
{
    for(std::size_t m = i * block_size; m < block_size * (i + 1); ++m)
    {
        for(std::size_t n = j * block_size; n < block_size * (j + 1); ++n)
        {
            mData[m][n] = subMat_.mData[m][n];
        }
    }
}

void Matrix::printSubMatrix(std::size_t block_size, std::size_t i, std::size_t j)
{
    for(std::size_t m = i * block_size; m < block_size * (i + 1); ++m)
    {
        for(std::size_t n = j * block_size; n < block_size * (j + 1); ++n)
        {
            std::cout << mData[m][n] << ' ';
        }
        std::cout << '\n';
    }
}

void Matrix::setSubMatrix(const Matrix& src, Matrix& subMat_, std::size_t block_size, std::size_t i, std::size_t j)
{
    for(std::size_t m = i * block_size; m < block_size * (i + 1); ++m)
    {
        for(std::size_t n = j * block_size; n < block_size * (j + 1); ++n)
        {
            subMat_.mData[m][n] = src.mData[m][n];
        }
    }
}

Matrix Matrix::setSubMatrix(size_t block_size, size_t i, size_t j)
{
    Matrix subMat_(block_size, block_size);
    for(std::size_t m = i * block_size; m < block_size * (i + 1); ++m)
    {
        for(std::size_t n = j * block_size; n < block_size * (j + 1); ++n)
        {
            subMat_.mData[m][n] = mData[m][n];
        }
    }
    return subMat_;
}

// Rectangle subMatrix
void Matrix::setSubMatrix(const Matrix& src, Matrix& subMat_, std::size_t i, std::size_t j, std::size_t mrows, std::size_t ncols)
{
    for(std::size_t m = i * mrows; m < mrows * (i + 1); ++m)
    {
        for(std::size_t n = j * ncols; n < ncols * (j + 1); ++n)
        {
            subMat_.mData[m][n] = src.mData[m][n];
        }
    }
}

void Matrix::setSubMatrix(const Matrix& srcMat_, std::size_t i, std::size_t j, std::size_t mrows, std::size_t ncols)
{
    for(std::size_t m = i * mrows; m < mrows * (i + 1); ++m)
    {
        for(std::size_t n = j * ncols; n < ncols * (j + 1); ++n)
        {
            mData[m][n] = srcMat_.mData[m][n];
        }
    }
}

Matrix Matrix::setSubMatrix(std::size_t i, std::size_t j, std::size_t mrows, std::size_t ncols)
{
    Matrix subMat_(mrows, ncols);
    for(std::size_t m = i * mrows; m < mrows * (i + 1); ++m)
    {
        for(std::size_t n = j * ncols; n < ncols * (j + 1); ++n)
        {
            subMat_.mData[m][n] = mData[m][n];
        }
    }
    return subMat_;
}

void Matrix::printSubMatrix(const Matrix& srcMat_, std::size_t i, std::size_t j, std::size_t mrows, size_t ncols )
{
    for(std::size_t m = i * mrows; m < mrows * (i + 1); ++m)
    {
        for(std::size_t n = j * ncols; n < ncols * (j + 1); ++n)
        {
            std::cout << srcMat_.mData[m][n] << ' ';
        }
        std::cout << '\n';
    }
}

void Matrix::printSubMatrix(std::size_t i, std::size_t j, std::size_t mrows, std::size_t ncols)
{
    for(std::size_t m = i * mrows; m < mrows * (i + 1); ++m)
    {
        for(std::size_t n = j * ncols; n < ncols * (j + 1); ++n)
        {
            std::cout << mData[m][n] << ' ';
        }
        std::cout << '\n';
    }
}

void Matrix::Identity()
{
    for(int row = 0; row < mNumRows; ++row)
    {
        for(int col = 0; col < mNumCols; ++col)
        {
            if(row == col) mData[row][col] = 1;
            else mData[row][col] = 0;
        }
        std::cout << '\n';
    }
}

void Matrix::Identity(const Matrix& mat)
{
    for(int row = 0; row < mat.mNumRows; ++row)
    {
        for(int col = 0; col < mat.mNumCols; ++col)
        {
            if(row == col) mat.mData[row][col] = 1;
            else mat.mData[row][col] = 0;
        }
        std::cout << '\n';
    }
}

void Matrix::free2DArray(double** arr, int row)
{
    for(int i = 0; i < row; i++)
    {
        delete[] arr[i];
    }
    delete[] arr;
}

double** Matrix::create2DArray(int row, int column)
{
    double** arr = new double*[row];
    if(!arr)
    {
        std::cout << "Unable to allocate memory" << std::endl;
        return nullptr;
    }

    for(int i = 0; i < row; i++)
    {
        arr[i] = new double[column];
        if(!arr[i])
        {
            free2DArray(arr, i);
            return nullptr;
        }
    }

    return arr;
}

Matrix Matrix::makeAugementedMatrix(const Matrix& mat1, const Matrix& mat2)
{
    Matrix augmentedMatrix(mat1.mNumRows, mat1.mNumCols + mat2.mNumCols);
    if(!augmentedMatrix.mData)
    {
        augmentedMatrix.mData = nullptr;
    }

    // Fill in the corresponding value to the augmented matrix
    for(int i = 0; i < augmentedMatrix.mNumRows; i++)
    {
        for(int j = 0; j < augmentedMatrix.mNumCols; j++)
        {
            if(j < mat1.mNumCols)
            {
                augmentedMatrix.mData[i][j] = mat1.mData[i][j];
            }
            else
            {
                augmentedMatrix.mData[i][j] = mat2.mData[j - mat1.mNumCols][i];
            }
        }
    }

    return augmentedMatrix;
}

void Matrix::swapRowMatrix(const Matrix& mat, int row1, int row2)
{
    for(int i = 0; i < mat.mNumCols; ++i)
    {
        double temp = mat.mData[row1][i];
        mat.mData[row1][i] = mat.mData[row2][i];
        mat.mData[row2][i] = temp;
    }
}

/* Reduce matrix to ref */
int Matrix::forwardElimination(Matrix& mat)
{
    for(int k = 0; k < mat.mNumRows; ++k)
    {
        // Initialize maximum value and index for pivot
		int i_max = k;
		int v_max = mat.mData[i_max][k];

		/* find greater amplitude for pivot if any */
		for (int i = k + 1; i < mat.mNumRows; ++i)
			if (abs(mat.mData[i][k]) > v_max)
				v_max = mat.mData[i][k], i_max = i;

		/* if a principal diagonal element is zero, it denotes that matrix is singular, and will lead to a division-by-zero later. */
		if (!mat.mData[k][i_max])
			return k; // Matrix is singular

		/* Swap the greatest value row with current row */
		if (i_max != k)
			swapRowMatrix(mat, k, i_max);


		for (int i = k + 1; i < mat.mNumRows; ++i)
		{
			/* factor f to set current row kth element to 0, and subsequently remaining k-th column to 0 */
			double f = mat.mData[i][k]/mat.mData[k][k];

			/* subtract fth multiple of corresponding k-th row element*/
			for (int j = k + 1; j <= mat.mNumRows; ++j)
				mat.mData[i][j] -= mat.mData[k][j]*f;

			/* filling lower triangular matrix with zeros*/
			mat.mData[i][k] = 0;
		}

		mat.display();	 //for matrix state
	}
	mat.display();		 //for matrix state
	return -1;

}

void Matrix::backSubstitution(const Matrix& mat)
{
    Matrix temp(mat.mNumRows, mat.mNumCols);

    /* Starting calculating from last equation up to the first */
    for(int i =  mat.mNumRows - 1; i >= 0; --i)
    {
        /* Start with the RHS of equation */
        temp.mData[i][0] = mat.mData[i][mat.mNumRows];

        /* Initialize j to i+1 since matrix is upper triangular */
        for(int j = i + 1; j < mat.mNumRows; ++j)
        {
            /* Subtract all LHS values except the coefficient of the variable whose value is being calculated */
            temp.mData[i][0] -= mat.mData[i][j]*temp.mData[j][0];
        }
        /* divide the RHS by the coefficient of the unknown being calculated */
		temp.mData[i][0] = temp.mData[i][0]/mat.mData[i][i];
    }
    /* Print the output of the solution (just the solution so this under Vector form) */
    std::cout << "This is the output of the solution (just the solution so this under Vector form): " << '\n';
    for(int i = 0; i < mat.mNumRows; ++i)
        std::cout << temp.mData[i][0] << '\n';
}

Matrix Matrix::backSubstitution()
{
    Matrix temp(mNumRows, mNumCols);

    /* Starting calculating from last equation up to the first */
    for(int i =  mNumRows - 1; i >= 0; --i)
    {
        /* Start with the RHS of equation */
        temp.mData[i][0] = mData[i][mNumRows];
        std::cout << "here1: " << '\n'
                  << temp.mData[i][0] << '\n';
        /* Initialize j to i+1 since matrix is upper triangular */
        for(int j = i + 1; j < mNumRows; ++j)
        {
            /* Subtract all LHS values except the coefficient of the variable whose value is being calculated */
            temp.mData[i][0] -= mData[i][j]*temp.mData[j][0];
            std::cout << "here2: " << '\n';
            std::cout << temp.mData[i][0] << '\n';
        }
        /* divide the RHS by the coefficient of the unknown being calculated */
		temp.mData[i][0] = temp.mData[i][0]/mData[i][i];
		std::cout << "here3: " << '\n'
                  << temp.mData[i][0] << '\n';
    }
    //Print the output of the solution (just the solution so this under Vector form)
    std::cout << "This is the output of the solution (just the solution so this under Vector form): " << '\n';
    for(int i = 0; i < mNumRows; ++i)
        std::cout << temp.mData[i][0] << '\n';
    return temp;
}

void Matrix::gaussianElimination(Matrix& mat)
{
    /* reduction into ref */
    int singular_flag = mat.forwardElimination(mat);

    /* if matrix is singular */
    if(singular_flag != -1)
    {
        std::cout << "Singular matrix/ Non-invertible matrix" << '\n';

        /* if RHS of equation corresponding to zero row is 0, linear system has infinitely many solutions, else inconsistent */
        if(mat.mData[singular_flag][mat.mNumRows])
            std::cout << "Inconsistent system." << '\n';
        else
            std::cout << "May have infinitely many solutions" << '\n';

        return;
    }
    /* Get solution to system and output the solution using backward substitution */
    backSubstitution(mat);
}

Matrix Matrix::gaussianElimination()
{
    /* reduction into ref */
    int singular_flag = forwardElimination(*this);

    /* if matrix is singular */
    if(singular_flag != -1)
    {
        std::cout << "Singular matrix/ Non-invertible matrix" << '\n';

        /* if RHS of equation corresponding to zero row is 0, linear system has infinitely many solutions, else inconsistent */
        if(this->mData[singular_flag][this->mNumRows])
            std::cout << "Inconsistent system." << '\n';
        else
            std::cout << "May have infinitely many solutions" << '\n';

        exit(EXIT_FAILURE);
    }
    /* Get solution to system and output the solution using backward substitution */
    std::cout <<  "This is the output of the solution (just the solution so this under Vector form): " << '\n';
    return this->backSubstitution();
}

// Conjugate gradient Method
bool Matrix::isSymmetric(const Matrix& mat)
{
    for(int i = 0; i < mat.mNumRows; ++i)
    {
        for(int j = 0; j < mat.mNumCols; ++j)
        {
            if(mat.mData[i][j] != mat.mData[j][i])
                return false;
        }
    }
    return true;
}

bool Matrix::isSymmetric()
{
    for(int i = 0; i < mNumRows; ++i)
    {
        for(int j = 0; j < mNumCols; ++j)
        {
            if(mData[i][j] != mData[j][i])
                return false;
        }
    }
    return true;
}

Matrix Matrix::uminus(Matrix& mat)
{
    Matrix minusMatrix(mat.mNumRows, mat.mNumCols);
    for(int i = 0; i < mat.mNumRows; ++i)
        for(int j = 0; j < mat.mNumCols; ++j)
            minusMatrix.mData[i][j] = -mat.mData[i][j];
    return minusMatrix;
}

Matrix Matrix::uminus()
{
    Matrix minusMatrix(mNumRows, mNumCols);
    for(int i = 0; i < mNumRows; ++i)
        for(int j = 0; j < mNumCols; ++j)
            minusMatrix.mData[i][j] = -mData[i][j];
    return minusMatrix;
}
