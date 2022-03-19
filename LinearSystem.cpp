#include "LinearSystem.hpp"
#include "Matrix_.hpp"
#include "Vector_.hpp"
#include <iostream>
#include <initializer_list>
#include <stdexcept>
#include <cmath>
#include <cassert>
#include <iomanip>

LinearSystem::LinearSystem() :  mSize(0), mpA(nullptr), mpb(nullptr) {}

LinearSystem::~LinearSystem()
{
    delete mpA;
    delete mpb;
}

LinearSystem::LinearSystem(int row, int col, double** pM, unsigned int Size, long double* pV, int linSize, Matrix* var1, Vector* var2)
    :  Matrix(pM, row, col), Vector(pV, Size), mSize(linSize), mpA(nullptr), mpb(nullptr)
{
    Matrix mat(pM, row, col);
    mpA = new Matrix(mpA->mNumRows, mpA->mNumCols);
    Vector vec(pV, Size);
    mpb = new Vector(mpb->mSize);
    mpA = var1;
    mpb = var2;
}

// Copy matrix and vector so that original matrix and vector specified are unchanged by Gaussian Elimination
LinearSystem::LinearSystem(const Matrix& mat, const Vector& vec)
    : Matrix(mat),
    Vector(vec)
{
    //std::cout << "Test for copy constructor of LinearSystem class object: " << '\n';
    //mat.display();
    // Check matrix and vector are of compatible sizes
    int local_size =  mat.getmNumRows();
    try
    {
        /*
        if(mat.getmNumCols() != local_size)
            throw std::range_error("Error! This is not a square matrix!"); */

        if(vec.getSize() != local_size)
            throw std::logic_error("Error! Can not form a correct size of linear system!");

        mSize = local_size;
        mpA = new Matrix(mat);
        mpb = new Vector(vec);
    }
    catch(const std::logic_error& le)
    {
        std::cerr << "You should look for the size we warning!" << '\n';
        exit(EXIT_FAILURE);
    }
    /* Since we do not need the matrix is to be square since we just update the SolveAll() method for non-linear matrix
    catch(const std::range_error& re)
    {
        std::cerr << "You should look for the size we warning!" << '\n';
        exit(EXIT_FAILURE);
    }*/
}

void LinearSystem::display() const
{
    Matrix mat(this->mpA->mData, this->mpA->mNumRows, this->mpA->mNumCols);
    Vector vec(this->mpb->mData, this->mpb->mSize);
    for(int i = 0; i < mat.mNumRows; ++i){
        for(int j = 0; j < mat.mNumCols; ++j){
            std::cout << std::setw(10) << mat.mData[i][j] ;
        }
        std::cout << '\n';
    }
    std::cout << "Vector in LinearSystem method: " << '\n';
    for(std::size_t i = 0; i < vec.mSize; ++i){
        std::cout << std::setw(10) << vec.mData[i] ;
    }
    std::cout << '\n';
}

// Not done this inheritance initializer list yet
LinearSystem::LinearSystem(std::initializer_list<std::initializer_list<double> > initMat, std::initializer_list<long double> initVec, std::initializer_list<std::initializer_list<double> > initLin)
    : Matrix(initMat), Vector(initVec), mSize(0), mpA(nullptr), mpb(nullptr)
{
    mSize = initLin.size();

}

LinearSystem::LinearSystem(std::initializer_list<std::initializer_list<double> > init_mat, std::initializer_list<long double> init_vec)
    : Matrix(init_mat), Vector(init_vec)
{
    try
    {
        // The iter pointer point to each comma in the initializer list(if 3x3 there are only 2 commas, and 3 child init_mat in parent init_mat)
        // This iter -> each row via comma
        for(auto iter = init_mat.begin(); iter < init_mat.end() - 1; ++iter)
        {
            const size_t _mSize = iter->size();

            if(_mSize == (iter+1)->size()) continue;
            else
            {
                throw "Your input value is not a matrix!";
            }
        }
    }
    catch(const char* msg)
    {
        std::cout << msg << '\n';
        exit(EXIT_FAILURE);
    }

    // Memory allocate on the heap
    mpA = new Matrix(init_mat.size(), init_mat.begin()->size());

    int i = 0, j = 0;
    // Traverse through each row of the parent initializer list
    for(auto iter = init_mat.begin(); iter < init_mat.end(); ++iter, ++i)
    {
        // Traverse through each col like each member in the child initializer list
        for(auto _iter = iter->begin(); _iter < iter->end(); ++_iter, ++j)
        {
            mpA->mData[i][j] = *_iter;
        }
        j = 0;
    }
    try
    {
        if(init_mat.begin()->size() != init_vec.size())
            throw std::length_error("The size of matrix and the size of vector is not equal. Please input again!");
    }
    catch(const std::length_error& le)
    {
        std::cerr << le.what() << '\n';
        exit(EXIT_FAILURE);
    }
    // Allocate memory on the heap
    mpb = new Vector(init_vec.size());

    int k = 0;
    for(auto iter = init_vec.begin(); iter < init_vec.end(); ++iter, ++k)
    {
        mpb->mData[k] = *iter;
    }
}

void LinearSystem::printVec_linear(const LinearSystem& ls)
{
    if(!mpb)
    {
        std::cerr << "The pointer is now NULL!";
        mpb->mData = nullptr;
    }
    for(std::size_t i = 0; i < ls.mpb->mSize; ++i)
        std::cout << std::setw(10) << ls.mpb->mData[i] << '\n';
}

/*
void LinearSystem::print_initList(std::initializer_list<LinearSystem> init)
{
    if(!mpb)
    {
        std::cerr << "The pointer is now Null!";
        mpb->mData = nullptr;
    }
    int i = 0;
    for(auto iter = init.begin(); iter < init.end(); ++iter, ++i)
        std::cout << *iter << '\n';
}
*/

Matrix LinearSystem::AugmentedMatrixVector(const Matrix& mat, const Vector& vec)
{
    Matrix _solution(this->mpA->mNumRows, this->mpA->mNumCols + 1);
    if(!_solution.mData)
    {
        _solution.mData = nullptr;
    }
    for(int i = 0; i < _solution.mNumRows; ++i)
    {
        for(int j = 0; j < _solution.mNumCols; ++j)
        {
            if(j < mat.mNumCols)
                _solution.mData[i][j] = mat.mData[i][j];
            else
                _solution.mData[i][j] = vec.mData[i];
        }
    }

    return _solution;
}

Vector LinearSystem::SolveAll()
{
    Vector sol(this->mpb->mData, this->mpA->mNumRows);
    Matrix solution(this->mpA->mData, this->mpA->mNumRows, this->mpA->mNumCols);
    Matrix solution_(this->mpA->mNumRows, this->mpA->mNumCols + 1);

    if(solution.mNumRows == solution.mNumCols)
    {
        solution_ = AugmentedMatrixVector(solution, sol);
        /*
        std::cout << "Output of columns of solution " << '\n'
                  << solution_.mNumCols << '\n'
                  << "Output of row of solution " << '\n'
                  << solution_.mNumRows << '\n';
        */

        // The solution_ will return a vector and then the solution Matrix have the vector at the first column of matrix
        //solution = solution_.gaussianElimination();
        solution_.gaussianElimination();

        /*
        std::cout << "Solution of gauss method: " << '\n'
                  << solution_ << '\n';
        */

        for(int i = 0; i < solution.mNumRows; ++i)
            sol.mData[i] = solution_.mData[i][0];
    }
    else{
        //solution_ = AugmentedMatrixVector(solution, sol);
        sol = solution.geninv() * sol;
    }
    return sol;
}


// Solve linear system using Gaussian elimination
// This method changes the content of the matrix mpA
Vector LinearSystem::Solve()
{
    Vector m(mSize); //See description in Appendix A
    Vector solution(mSize);

    // We introduce references to make the syntax readable
    Matrix rA = *mpA;
    Vector rb = *mpb;

    // Forward sweep of Gaussian elimination
    for (int k = 0; k < mSize - 1; k++)
    {
        // See if pivoting is necessary
        double max = 0.0;
        int row = -1;
        for (int i = k; i < mSize; i++)

        {
            if (fabs(rA(i + 1, k + 1)) > max)
            {
                row = i;
            }
        }
        assert(row > 0);

        // Pivot if necessary
        if (row != k)
        {
            // Swap matrix rows k+1 with row+1
            for (int i = 0; i < mSize; i++)
            {
                double temp = rA(k + 1, i + 1);
                rA(k + 1, i + 1) = rA(row + 1, i + 1);
                rA(row + 1, i + 1) = temp;
            }
            // Swap vector entries k+1 with row+1
            double temp = rb(k + 1);
            rb(k + 1) = rb(row + 1);
            rb(row + 1) = temp;
        }

        // Create zeros in lower part of column k
        for (int i = k + 1; i < mSize; i++)
        {
            m(i + 1) = rA(i + 1, k + 1) / rA(k + 1, k + 1);
            for (int j = k; j < mSize; j++)
            {
                rA(i + 1, j + 1) -= rA(k + 1, j + 1) * m(i + 1);
            }
            rb(i + 1) -= rb(k + 1) * m(i + 1);
        }
    }

    // Back substitution
    for (int i = mSize - 1; i > -1; i--)
    {
        solution(i + 1) = rb(i + 1);
        for (int j = i + 1; j < mSize; j++)
        {
            solution(i + 1) -= rA(i + 1, j + 1) * solution(j + 1);
        }
        solution(i + 1) /= rA(i + 1, i + 1);
    }

    return solution;
}

