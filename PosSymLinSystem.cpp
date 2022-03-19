#include "PosSymLinSystem.hpp"
#include "Matrix_.hpp"
#include "Vector_.hpp"
#include "LinearSystem.hpp"
#include <iostream>
#include <cctype>
#include <iomanip>

PosSymLinSystem::PosSymLinSystem()
    : LinearSystem() { }

PosSymLinSystem::~PosSymLinSystem() {}

double PosSymLinSystem::norm(Matrix& mat)
{
    double sum = 0.0;

    for(int i = 0; i < mat.mNumRows; ++i)
        sum += mat.mData[i][0]*mat.mData[i][0];

    return std::sqrt(sum);
}

PosSymLinSystem::PosSymLinSystem(const Matrix& mat, const Vector& vec)
 : LinearSystem(mat, vec)
{
    // Check matrix and vector are of compatible sizes
    int local_size =  mat.getmNumRows();
    try
    {
        if(mat.getmNumCols() != local_size)
            throw std::range_error("Error! This is not a square matrix!");

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
    catch(const std::range_error& re)
    {
        std::cerr << "You should look for the size we warning!" << '\n';
        exit(EXIT_FAILURE);
    }
}

double PosSymLinSystem::norm()
{
    double sum = 0.0;

    for(int i = 0; i < this->mpA->mNumRows; ++i)
        sum += this->mpA->mData[i][0]*this->mpA->mData[i][0];

    return std::sqrt(sum);
}

// New Level-3 BLAS Kernels for Cholesky Factorization
/*
 DPOTRF computes the Cholesky factorization of a real symmetric
 DPOTF2 is applied in this method
 positive definite matrix A.

 The factorization has the form
    A = U**T * U,  if UPLO = 'U', or
    A = L  * L**T,  if UPLO = 'L',
 where U is an upper triangular matrix and L is lower triangular.

 This is the block version of the algorithm, calling Level 3 BLAS.

 However, in this snippet code only applying sparse upper triangular matrix( U )

 Note:
 ajj is variable alike to index a[j][j]
 ddot is external variable, something like double precision
 NM1 is variable like the second last dimension
 JM1 is the dimension for column - 1
 dum1 is worked like container variable
 */
// Before use this make sure MATRIX IS A SQUARE MATRIX
void PosSymLinSystem::dpotrf_(const Matrix& A_Matrix, int* info)
{
    // For this program, it is assumed that UPLO = 1, so it will be hard-coded to use the upper
	// triangular part of the matrix. And the code will be unblocked: use DPOTF2.

	// On output:
	// info = 0:		successful exit
	// info = k > 0:	the leading minor of order k is not positive definite,
	//					and the factorization could not be completed.

	int NM1 = A_Matrix.mNumRows - 1,  i = NM1, j, JM1, jj;
	double ajj = A_Matrix.mData[0][0], ddot, dum1;

	// BEGIN DPOTF2

	// Deal with j = 0 case outside of main loop

	if (ajj <= 0){
		*info = 1;
		return;
	}

	A_Matrix.mData[0][0] = ajj = sqrt(ajj);

	while (i > 0) A_Matrix.mData[0][i--] /= ajj;

	// Now deal with the rest of the j cases

	for (j = 1; j < A_Matrix.mNumRows; ++j){

		ajj = A_Matrix.mData[j][j];
		ddot = 0.0;

		// DDOT
		i = j;
		do {
			--i;
			dum1 = A_Matrix.mData[i][j];
			ddot += dum1 * dum1;
		} while (i > 0);
		// End DDOT

		ajj -= ddot;

		if (ajj <= 0){
			A_Matrix.mData[j][j] = ajj;
			*info = j + 1;
			return;
		}

		A_Matrix.mData[j][j] = ajj = std::sqrt(ajj);

		// Compute elements J+1 to N of row J

		if (j < NM1){

			// BEGIN DGEMV

			// Start the operations. In this version, the elements of A are accessed sequentially
			// with one pass through A.

			// Since this program assumes we are working with the upper diagonal matrix of A,
			// form y = alpha*A^T * x + y

			JM1 = j - 1;
			dum1 = A_Matrix.mData[JM1][j];

			for (jj = NM1; jj > j; --jj){

				i = JM1;
				ddot = dum1 * A_Matrix.mData[i][jj];

				while (i > 0){
					--i;
					ddot += A_Matrix.mData[i][jj] * A_Matrix.mData[i][j];  // A * X
				}

				A_Matrix.mData[j][jj] = (A_Matrix.mData[j][jj] - ddot) / ajj;

			} // End for jj

			// END DGEMV

		} // End if (j < NM1)

	} // End for j

	// END DPOTF2

	return;
}

bool PosSymLinSystem::CholeskyDecomposition(const Matrix& MatrixA)
{
    int matrixDim = MatrixA.mNumRows;

    int decimal;
    std::cout << "The results are calculated to double precision-- decimal places: ";
    std::cin >> decimal;

    if(MatrixA.mNumRows < 1)
    {
        std::cout << "\nInvalid dimension entered. Program terminated." << '\n';
        return false;
        exit(EXIT_FAILURE);
    }

    int info = 0;
    bool erflag = false;


    if(erflag == this->mpA->isSymmetric(MatrixA))
    {
        std::cout << "Non-symmetry in matrix detected. No further action." << '\n';
        return false;
    }

    if((matrixDim-1 > 0) && (erflag != isSymmetric(MatrixA)))
    {
        Matrix MatrixB(MatrixA.mNumRows, MatrixA.mNumCols);

        for(int i = 0; i < MatrixA.mNumCols; ++i)
            for(int j = i; j < MatrixA.mNumCols; ++j)
                MatrixB.mData[i][j] = MatrixA.mData[i][j];

        dpotrf_(MatrixB, &info);

        if(info != 0)
        {
            std::cout << "The leading is not minor of order " << info << "is not positive definite, "
                      << "and the factorization could not be completed." << '\n';
            std::cout << "the " << info << "th argument had an illegal value" << '\n';
            return false;
        }
        else // if(info == 0)
        {
            std::cout << "The factorization was completed successfully." << '\n';
            std::cout << std::ios_base::fixed << std::ios_base::showpoint;
            std::cout << std::setprecision(decimal);

            std::cout << "The factored matrix follows: " << '\n';

            for(int i = 0; i < MatrixB.mNumRows; ++i)
                for(int j = 0; j < MatrixB.mNumCols; ++j)
                    std::cout << MatrixB.mData[i][j] << '\n';
            return true;
        }
    }
    return true;
}

/*
Matrix& PosSymLinSystem::ConjugateGradientMethodNonLinearMatrix(Matrix& _A, Matrix& _b, Matrix& _x, double _tol, unsigned _maxIter)
{
    // Cholesky decomposition seems not correct
    if(_A.mNumRows > _A.mNumCols)
    {

    }
}
*/

/*
Note:
    _A is source matrix
    _b is vector of solution
    _x is x1,x2,...,.. is vector of x variable
    f stands for gradient direction or function  phi: -g(k)
    l stands for conjugate direction or s(k)
    d0 stands for -(Ax(k) - b) also negative of gradient direction abbreviate for g(k)
    pom stands for x^T * b
    matrix_product stands for g(k)^T * A * g(k)
    l0 = g(k) * matrix_product^T
    d = d0 use for not change the d variable
    inv_mat = matrix_product^-1
    matrix_product2 = T_d * _A * d;
    inv_mat2 = matrix_product2^-1
    mi = ( f^T * _A * d ) * inv_mat2
    d = mi * d + g(k)
    pom1 = (-g(k))^T * d
    matrix_product3 = T_d * _A * d;
    inv_mat3 = matrix_product^-1
    l = -( (-g(k))^T * d ) * inv_mat3
    _x += l * g(k)

Constraints:
    max iteration always smaller than 10
    norm(f) or norm of function phi < tolerance(_tol)
*/
Matrix PosSymLinSystem::ConjugateGradientMethodLinearMatrix(Matrix& _A, Matrix& _b, Matrix& _x, double _tol, unsigned _maxIter)
{
        /*
        std::cout << "Output for matrix A inside PosSymLinSystem: " << '\n'
                  << _A << '\n'
                  << "Output for matrix b inside PosSymLinSystem: " << '\n'
                  << _b << '\n';*/
        if(!isSymmetric(_A))
        {
            std::cout << "Matrix is not symmetric" << '\n';
            return Matrix(0, 0);
        }
        else{
            /*
            //Matrix product1(_A.mNumRows, _b.mNumCols);
            //Matrix f0(_A.mNumRows, _b.mNumCols);

            //product1.matrixMultiply(_A, _x , product1);*/
            Matrix f0 = _A * _x - _b;
            /*
            //f0 = product1 - _b;

            //Matrix d0(_A.mNumRows, _b.mNumCols); */
            Matrix d0 = uminus(f0);
            /*
            //Matrix T_f0(_b.mNumCols, _A.mNumRows);
            //this->mpA->transposeMatrix(f0, T_f0);

            //Matrix pom(_b.mNumCols, _b.mNumCols);
            //pom.matrixMultiply(T_f0, d0, pom); // correct */
            Matrix pom = (f0^'T') * d0;
            /*
            Matrix T_d0(_b.mNumCols, _A.mNumRows);
            d0.transposeMatrix(d0, T_d0); // correct

            Matrix matrix_product(_b.mNumCols, _b.mNumCols);
            Matrix temp(_b.mNumCols, _A.mNumRows);

            //matrix_product = T_d0 * _A * d0;
            //temp.matrixMultiply(T_d0, _A, temp); // correct
            //matrix_product.matrixMultiply(temp, d0, matrix_product); // correct
            matrix_product = T_d0 * _A * d0;

            Matrix inv_mat(_b.mNumCols, _b.mNumCols);
            matrixInverse(matrix_product, inv_mat); // correct

            Matrix l0(_b.mNumCols, _b.mNumCols);
            //l0 = uminus(pom) * inv_mat;

            l0.matrixMultiply(pom.uminus(), inv_mat, l0); // correct */
            Matrix l0 = uminus(pom) * ( (d0^'T') * _A * d0)^-1;
            _x = _x + d0 * l0;
            /*
            Matrix temp2(_A.mNumRows, _b.mNumCols);
            temp2.matrixMultiply(d0, l0, temp2);
            _x = _x + temp2; // correct
            Matrix d(_A.mNumRows, _b.mNumCols);*/
            Matrix d = d0;

            // MaxIter should be in the range from _A.mNumRows to _A.mNumRows+1
            for(unsigned i = 1; i <= _maxIter; ++i)
            {
                Matrix f = _A *_x - _b;
                /*
                Matrix f(_A.mNumRows, _b.mNumCols);
                f.matrixMultiply(_A, _x, f);
                f = f - _b; // correct

                if(norm(f) < _tol)
                    return _x; // correct*/
                Matrix mi = ( (f^'T') * _A * d) * ( ( (d^'T') * _A * d)^-1 );
                /*
                Matrix T_f(_b.mNumCols, _A.mNumRows);
                f.transposeMatrix(f, T_f); // correct
                Matrix T_d(_b.mNumCols, _A.mNumRows);
                d.transposeMatrix(d, T_d); // correct
                Matrix matrix_product2(_b.mNumCols, _b.mNumCols); //correct
                Matrix temp3(_b.mNumCols, _A.mNumRows); //correct

                //matrix_product2 = T_d * _A * d;
                temp3.matrixMultiply(T_d, _A, temp3);

                matrix_product2.matrixMultiply(temp3, d, matrix_product2);

                Matrix inv_mat2(_b.mNumCols, _b.mNumCols);
                matrixInverse(matrix_product2, inv_mat2);

                //mi = (T_f * _A * d) * inv_mat2;
                Matrix mi(_b.mNumCols, _b.mNumCols);
                Matrix temp4(_b.mNumCols, _A.mNumRows);
                temp4.matrixMultiply(T_f, _A, temp4);
                Matrix temp5(_b.mNumCols, _b.mNumCols);
                temp5.matrixMultiply(temp4, d, temp5);
                mi.matrixMultiply(temp5, inv_mat2, mi); // correct*/
                d = mi*d - f;
                /*
                d.matrixMultiply(d, mi, d);
                d = d - f; // correct*/
                Matrix pom1 = f.transposeMatrix() * d;
                /*
                Matrix pom1(_b.mNumCols, _b.mNumCols);
                pom1.matrixMultiply(T_f, d, pom1); // correct*/
                Matrix l = (uminus(pom1)) * ( (d.transposeMatrix() * _A * d)^-1 );
                /*
                //matrix_product3 = T_d * _A * d;
                Matrix matrix_product3(_b.mNumCols, _b.mNumCols);
                Matrix temp6(_b.mNumCols, _A.mNumRows);

                d.transposeMatrix(d, T_d);
                temp6.matrixMultiply(T_d, _A, temp6);
                matrix_product3.matrixMultiply(temp6, d, matrix_product3);

                Matrix inv_mat3(_b.mNumCols, _b.mNumCols);
                matrixInverse(matrix_product3, inv_mat3);*/
                _x = _x + l * d;
                /*
                //l = (uminus(pom1)) * inv_mat3;
                Matrix l(_b.mNumCols, _b.mNumCols);

                l.matrixMultiply(pom1.uminus(), inv_mat3, l); // correct

                Matrix temp7(_A.mNumRows, _b.mNumCols);
                temp7.matrixMultiply(d, l, temp7);

                _x = _x + temp7; // correct*/
                std::cout << "Output for _x: " << '\n'
                        << _x << '\n';
            }
            return _x; // correct
        }
}

Vector PosSymLinSystem::SolveAll()
{
    Vector _b(this->mpb->mData, this->mpA->mNumRows);
    Matrix _A(this->mpA->mData, this->mpA->mNumRows, this->mpA->mNumCols);
    Matrix _x(this->mpA->mData, this->mpA->mNumRows, 1);
    Vector _sol(this->mpb->mData, this->mpA->mNumRows);
    Matrix _bMatrix(this->mpA->mData, this->mpA->mNumRows, 1);

    for(int i = 0; i < _bMatrix.mNumRows; ++i)
        _bMatrix.mData[i][0] = _b.mData[i];

    for(int i = 0; i < _x.mNumRows; ++i)
        _x.mData[i][0] = -_b.mData[i];

    unsigned maxIter;
    std::cout << "Enter the max number of iterations" << '\n';
    std::cin >> maxIter;

    double tol;
    std::cout << "Enter tolerance: " << '\n';
    std::cin >> tol;

    ConjugateGradientMethodLinearMatrix(_A, _bMatrix, _x, tol, maxIter);

    for(int i = 0; i < _A.mNumRows ; ++i)
        _sol.mData[i] = _x.mData[i][0];

    std::cout << "Solution of matrix: " << '\n' << _sol << '\n';

    return _sol;
}

Vector PosSymLinSystem::Solve()
{
    // Pseudo-Penrose_Moore method
    Vector _b(this->mpb->mData, this->mpb->mSize);
    Matrix _A(this->mpA->mData, this->mpA->mNumRows, this->mpA->mNumCols);
    Matrix _AMatrix(this->mpA->mNumRows, this->mpA->mNumCols);
    Matrix _solution(this->mpA->mNumRows, this->mpA->mNumCols + 1);

    Matrix _bMatrix(this->mpA->mNumRows,1);

    Vector _sol(this->mpb->mData, this->mpb->mSize);

    // over-determined linear system
    if(_A.mNumRows > _A.mNumCols)
    {
        _solution = AugmentedMatrixVector(_A, _b);
        std::cout << "Output for the whole linear system : " << '\n'
                  << _solution << '\n';

        std::cout << "Solution of over-determined linear system: " << '\n';
        _AMatrix = _A.pseudoInverse();
        std::cout << _AMatrix << '\n';

        for(unsigned i = 0; i < _b.mSize; ++i)
            _bMatrix.mData[i][0] = _b.mData[i];

        _solution = _AMatrix * _bMatrix;

        for(int i = 0; i < _A.mNumRows; ++i)
            _sol.mData[i] = _solution.mData[i][0];

        std::cout << "Solution of vector _sol: " << '\n'
                  << _sol << '\n';
    }
    // under-determined linear system
    if(_A.mNumRows < _A.mNumCols)
    {
        _solution = AugmentedMatrixVector(_A, _b);

        std::cout << "Output for the whole linear system : " << '\n'
                  << _solution << '\n';

        std::cout << "Solution of under-determined just for A matrix: " << '\n';
        _AMatrix = _A.pseudoInverse();
        std::cout << _AMatrix << '\n';

        for(unsigned i = 0; i < _b.mSize; ++i)
            _bMatrix.mData[i][0] = _b.mData[i];

        _solution = _AMatrix * _bMatrix;

        for(int i = 0; i < _A.mNumRows; ++i)
            _sol.mData[i] = _solution.mData[i][0];

        std::cout << "Solution of vector _sol also the solution for linear system : " << '\n'
                  << _sol << '\n';
    }
    // square linear system
    if(_A.mNumRows == _A.mNumCols)
    {
        _solution = AugmentedMatrixVector(_A, _b);

        std::cout << "Output for the whole linear system : " << '\n'
                  << _solution << '\n';

        std::cout << "Solution of square matrix _A after inverse: " << '\n';
        _AMatrix = _A.pseudoInverse();
        std::cout << _AMatrix << '\n';

        for(unsigned i = 0; i < _b.mSize; ++i)
            _bMatrix.mData[i][0] = _b.mData[i];

        _solution = _AMatrix * _bMatrix;

        for(int i = 0; i < _A.mNumRows; ++i)
            _sol.mData[i] = _solution.mData[i][0];

        std::cout << "Solution of vector _sol also the solution for linear system : " << '\n'
                  << _sol << '\n';
    }
    return _sol;
}

/* Before using this method make sure that the matrix is non-linear */
/*
void PosSymLinSystem::ConjugateGradientMethodNonLinearMatrix(Matrix& _mat, Vector& _b, Vector& _x, double _tol, unsigned _maxIter)
{

}
*/
