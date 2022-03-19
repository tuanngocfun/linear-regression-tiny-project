#ifndef LINEAR_SYSTEM
#define LINEAR_SYSTEM
#include "Matrix_.hpp"
#include "Vector_.hpp"
#include <initializer_list>

// This class is somewhat alike to singleton design
// Assuming the system is nonsingular
class LinearSystem : protected Matrix, protected Vector
{
public:
    LinearSystem();
    LinearSystem(int row, int col, double** pM, unsigned int Size, long double* pV, int linSize, Matrix* var1, Vector* var2);
    LinearSystem(std::initializer_list<std::initializer_list<double> > init_mat, std::initializer_list<long double> init_vec);
    LinearSystem(std::initializer_list<std::initializer_list<double> > initMat, std::initializer_list<long double> initVec, std::initializer_list<std::initializer_list<double> > initLin);
    LinearSystem(const Matrix& mat, const Vector& vec);
    virtual ~LinearSystem();

    void printVec_linear(const LinearSystem& ls);
    void print_initList(std::initializer_list<LinearSystem> init);

    virtual void display() const override;

    // Merge matrix and vector to form a new matrix
    Matrix AugmentedMatrixVector(const Matrix& mat, const Vector& vec);

    //friend std::ostream& operator<<(std::ostream& out, const LinearSystem& obj);

    // Solve linear algebra
    virtual Vector SolveAll();
    virtual Vector Solve();

protected:
    int mSize;
    Matrix* mpA;
    Vector* mpb;

    friend class PosSymLinSystem;
private:
    // Only allow constructor that specifies matrix and vector to be used. Copy constructor is private
    LinearSystem(const LinearSystem& other);

};
#endif // LINEAR_SYSTEM

