#ifndef POSSYMLINSYSTEM_H
#define POSSYMLINSYSTEM_H
#include "LinearSystem.hpp"
#include "Vector_.hpp"
#include "Matrix_.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

class PosSymLinSystem : public LinearSystem
{
public:
    PosSymLinSystem();
    virtual ~PosSymLinSystem();
    PosSymLinSystem(const Matrix& mat, const Vector& vec);

    //  positive definite symmetric linear systems
    bool CholeskyDecomposition();
    bool CholeskyDecomposition(const Matrix& mat);
    void dpotrf_(const Matrix& mat, int* info);

    // Conjugate Gradient Method
    double norm(Matrix& mat);
    double norm();
    Vector ConjugateGradientMethodSquareMatrix(Matrix& _mat, Vector& _b, Vector& _x, double _tol, unsigned _maxIter);
    Matrix ConjugateGradientMethodSquareMatrix(Matrix& _mat, Matrix& _b, Matrix& _x, double _tol, unsigned _maxIter);
    Vector ConjugateGradientMethodNonLinearMatrix(Matrix& _mat, Vector& _b, Vector& _x, double _tol, unsigned _maxIter);
    Matrix ConjugateGradientMethodNonLinearMatrix(Matrix& _mat, Matrix& _b, Matrix& _x, double tol, unsigned _maxIter);
    Matrix ConjugateGradientMethodLinearMatrix(Matrix& _A, Matrix& _b, Matrix& _x, double tol, unsigned _maxIter);
    Vector Solve() override;

    // Solve linear system when over-determined or over-determined or square matrix
    Vector SolveAll() override;
};
#endif // POSSYMLINSYSTEM_H



