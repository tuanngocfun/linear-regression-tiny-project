#include "readFile.hpp"


using namespace std;

bool doubleComp(long double firstNum, long double secNum)
{
    long double err = std::fabs(firstNum - secNum);

    return (err < 1e-9) ? true : false;
}

int main(int argc, char** argv)
{
    /*++++++++++++++++Welcome everyone by choosing this source code for executing or implement some math function++++++++++++*/
    // Note: some function are supported to solve linear system, or some sub-functions to solve problems related to linear algebra
    // Please note that if you want to see how to solve the problems by type a syntax alike to my snippet code.
    // Here is some tips to do/implement it:
    // 1. Read carefully the syntax to implement each function
    // 2. Some functions are upgraded to make user-friendly interface (approx. line 300~490)
    // 3. There maybe some bugs, since I have not tested with a wide range of cases, please tell me if you found it by this mail below:
    // ## 15809@student.vgu.edu.vn
    // Introduction for testing the tiny project from scracth
    /*
    Vector myVector(5);
    Vector p1(myVector);

    cout << fixed << showpoint << setprecision(2);

    cout << myVector.getData() << '\n';
    cout << myVector.getSize() << '\n';

    Vector v1(5);
    cout << "output for v1 constructor:" << '\n';
    cout << v1 << endl;
    cout << "\n\n";

    long double d[3]{1.1, 1.2, 3.3};
    Vector v2(d, 3);
    long double d1[3]{2.0, 3.0, 4.0};
    Vector v6(d1, 3);
    cout << "output for v2 constructor:" << '\n';
    cout << v2 << endl;
    cout << "Output for vector multiply with a scalar" << '\n';
    cout << (v2*2) << '\n';
    cout << "Output for v6 constructor:" << '\n';
    cout << v6 << '\n';
    cout << "Output for vector multiply with a vector" << '\n';
    cout << (v2 * v6) << '\n';
    cout << "\n\n";

    v1(3) = 8.8;
    cout << "output for v1 after v(3) = 8.8 constructor:" << '\n';
    cout << v1(3) << endl;
    cout << "\n\n";

    Vector v3(v2);
    v3 = v1;
    cout << "Output for v1 after v3 = v1: " << '\n';
    cout << v1 << endl;
    cout << "output for v3 = v1 constructor:" << '\n';
    cout << v3 << endl;
    cout << "\n\n";


    unsigned int SIZE = 3;
    Vector v4(SIZE);
    cout << "Input the number into v4: " << '\n';
    for(unsigned int i = 1; i <= SIZE; ++i)
        cin >> v4(i);
    if(cin.good())
        {
            cout << "output for v4 with size = 3 constructor:" << '\n';
            cout << v4 << endl;
        }
    else
        cerr << "Invalid input" << endl;


    cout << "Output v1 before comparison: " << '\n';
    cout << v1 << '\n';
    Vector v5(5);
    for(unsigned int i = 1; i <= v5.Size(); ++i)
        cin >> v5(i);
    if(cin.good())
        {
            cout << "output for v5 with size = 5 constructor:" << '\n';
            cout << v5 << endl;
        }
    else
        cerr << "Invalid input" << endl;

    Vector u(1);
    u = v1;

    cout << "Output for comparison between v5(i) and v1(i): " << '\n';
    for(unsigned int i = 1; i <= v5.Size(); ++i)
    {
        cout << doubleComp(v5(i), v1(i)) << "\n";
    }
    cout << endl;
    */
    /*
    long double d[3]{1.1, 2.2, 3.3};
    Vector v1(d, 3);
    cout << "Output of v1 constructor: " << '\n';
    long double d1[4]{1.1, 2.2, 3.3, 4.4};
    Vector v2(d1, 4);
    cout << "Output of v2 constructor: " << '\n';
    cout << "Output for addition of v1 and v2: " << '\n';
    cout << (v2 + v1) << '\n';
    */
    /*
    long double d[3]{1.1, 2.2, 3.3};
    Vector v1(d, 3);

    Vector v2 = v1;
    cout << "Output for += btw 2 vector: " << '\n';
    cout << (v2 += v1) << '\n';
    cout << "Output fot v2 after += " << '\n';
    cout << v2 << '\n';
    cout << "Output fot v1 after += " << '\n';
    cout << v1 << '\n';

    Vector v3(3);
    v3(3) = 1.2;
    cout << "Output for *= btw 2 vector: " << '\n';
    cout << (v3 *= v1) << '\n';
    cout << "Output for v3 when v3 = v1 - v2" << '\n';
    cout << (v3 = v1 - v2) << '\n';
    cout << "Output for v3 when v3 = v1 + v2" << '\n';
    cout << (v3 = v1 + v2) << '\n';
    cout << "Output for v1 * v2" << '\n';
    cout << (v1 * v2) << '\n';
    cout << "Output for v1 <= v2 : " << '\n';
    cout << (v1 <= v2) << '\n';
    cout << "Output for v1 >= v2 : " << '\n';
    cout << (v1 >= v2) << '\n';
    */
    /*
    long double d[3]{1.1, 2.2, 3.3};
    Vector v1(d, 3);
    cout << "Output v1: " << '\n';
    cout << v1 << '\n';
    long double d1[3]{2.2, 1.1, 3.3};
    Vector v2(d1, 3);
    cout << "Output v2: " << '\n';
    cout << v2 << '\n';
    for(unsigned int i = 1; i <= v2.Size(); ++i)
        cout << (v1(i) >= v2(i)) << '\n';
    */
    /*
    long double d[3]{1.1, 2.2, 3.3};
    Vector v(d, 3), v1(d, 3);
    Vector v2(3);
    cout << "Output for v2 = v + v1: " << '\n';
    cout << (v2 = v + v1) << '\n';
    */
    /*
    double** d = new double*[3];
    for(int i = 0; i < 3; ++i) d[i] = new double[3];

    std::cout << "Input your number: " << '\n';
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
            cin >> d[i][j];
    Matrix m1(d, 3, 3);
    cout << "Output m1: " << '\n';
    cout << m1 << '\n';
    cout << "Output matrix multiply scalar:" << '\n';
    cout << (m1 * 2) << '\n';
    */
    /*
    long double d[3] {1.1, 2.2, 3.3};
    Vector v1(d, 3);
    cout << "Output of scalar multiply a vector: " << '\n';
    cout << (2 * v1) << '\n';
    */
    /*
    long double myArr[4];
    std::initializer_list<long double> initList = {1.1, 2.2, 3.3, 4.4};

    int i = 0;
    for(auto iter = initList.begin(); iter < initList.end(); ++iter, ++i)
        myArr[i] = *iter;

    for(auto& iter : initList)
        myArr[i] = iter;

    std::cout << "Output for index i: " << '\n';
    for(auto& i : myArr)
        std::cout << i << '\n';

    std::cout << "Output for vector v: " << '\n';
    Vector v(myArr, 4);
    std::cout << v << '\n';

    double** myMatrix = new double*[3];
    for(int c = 0; c < 3; ++c)
        myMatrix[c] = new double[3];

    std::initializer_list<std::initializer_list<double> > initlist = { {10.0, 20.0 , 10.0},{4.0, 5.0, 6.0},{2.0, 3.0, 5.0} };

    int j = 0, k = 0;
    for(auto iter = initlist.begin(); iter < initlist.end(); ++iter)
    {
        for(auto _iter = iter->begin(); _iter < iter->end(); ++_iter, ++j)
        {
            myMatrix[k][j] = *_iter;
        }
        j = 0;
        ++k;
    }
    std::cout << "Output for matrix m1: " << '\n';
    Matrix m1(myMatrix, 3, 3);
    cout << m1 << '\n';

    double** myMatrix1 = new double*[3];
    for(int c = 0; c < 3; ++c)
        myMatrix1[c] = new double[3];

    std::initializer_list<std::initializer_list<double> > initList_ = { {3.0, 2.0, 4.0},{3.0, 3.0, 9.0},{4.0, 4.0, 2.0} };

    int m = 0, n = 0;
    for(auto iter = initList_.begin(); iter < initList_.end(); ++iter)
    {
        for(auto _iter = iter->begin(); _iter < iter->end(); ++_iter, ++n)
        {
            myMatrix1[m][n] = *_iter;
        }
        n = 0;
        ++m;
    }
    std::cout << "Output for matrix m2: " << '\n';
    Matrix m2(myMatrix1, 3, 3);
    cout << m2 << '\n';

    cout << "Matrix multiplication between m1 and m2: " << '\n';
    cout << (m1 * m2) << '\n';
    cout << "Matrix multiplication between m1 and m2 via m1 *= operator: " << '\n';
    cout << (m1 *= m2) << '\n';

    cout << "Addition between 2 matrix m1 and m2: " << '\n';
    cout << (m1 + m2) << '\n';

    cout << "Addition between 2 matrix m1 and m2 via m1 += operator: " << '\n';
    cout << (m1 += m2) << '\n';
    cout << "Subtract 2 matrix m1 and m2 via m1 -= operator: " << '\n';
    cout << (m1 -= m2) << '\n';

    cout << "Matrix m1 multiply with a scalar: " << '\n';
    cout << (m1 * 2) << '\n';
    cout << "Matrix m1 multiply with a scalar via m1 *= operator: " << '\n';
    cout << (m1 *= 2) << '\n';

    cout << "Comparison with itself: " << '\n';
    cout << (m1 == m1) << '\n';
    cout << "Comparison between m1 and m2: " << '\n';
    cout << (m1 == m2) << '\n';

    std::initializer_list<std::initializer_list<double> > _initlist = { {10.0, 20.0 , 10.0},{4.0, 5.0, 6.0},{2.0, 3.0, 5.0} };

    Matrix m3(_initlist);
    cout << "Initializer list constructor for Matrix class m3: " << '\n';
    cout << m3 << '\n';

    Matrix m4{ {1.1,2.2,3.3},{4.4,5.5,6.6},{7.7,8.8,9.9} };
    cout << "Initializer list constructor for Matrix class m4: " << '\n';
    cout << m4 << '\n';
    */
    /*
    double** myMatrix = new double*[3];
    for(int c = 0; c < 3; ++c)
        myMatrix[c] = new double[3];

    std::initializer_list<std::initializer_list<double> > initlist = { {10.0, 20.0 , 10.0},{4.0, 5.0, 6.0},{2.0, 3.0, 5.0} };

    int j = 0, k = 0;
    for(auto iter = initlist.begin(); iter < initlist.end(); ++iter)
    {
        for(auto _iter = iter->begin(); _iter < iter->end(); ++_iter, ++j)
        {
            myMatrix[k][j] = *_iter;
        }
        j = 0;
        ++k;
    }
    Matrix m1(myMatrix, 3, 3);
    cout << m1 << '\n';

    m1.pseudoInverse();
    */
    /*
    double** myMatrix_ = new double*[3];
    for(int c = 0; c < 3; ++c)
        myMatrix_[c] = new double[2];

    std::initializer_list<std::initializer_list<double> > initList ( {{1,-1},{0,1},{1,0}} );

    int a = 0, b = 0;
    for(auto iter = initList.begin(); iter < initList.end(); ++iter, ++a)
    {
        for(auto _iter = iter->begin(); _iter < iter->end(); ++_iter, ++b)
        {
            myMatrix_[a][b] = *_iter;
        }
        b = 0;
    }

    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 2; j++)
        {
            cout << myMatrix_[i][j] << "  ";
        }
        cout << "\n";
    }
    cout << endl;


    Matrix m2(myMatrix_, 3, 2);

    std::cout << m2 << '\n';

    m2.pseudoInverse();
    m2.pseudoInverse(m2);
    */
    /*
    double** myMatrix_ = new double*[3];
    for(int c = 0; c < 3; ++c)
        myMatrix_[c] = new double[3];

    std::initializer_list<std::initializer_list<double> > initList { {10.0, 20.0 , 10.0},{4.0, 5.0, 6.0},{2.0, 3.0, 5.0} };

    int a = 0, b = 0;
    for(auto iter = initList.begin(); iter < initList.end(); ++iter, ++a)
    {
        for(auto _iter = iter->begin(); _iter < iter->end(); ++_iter, ++b)
        {
            myMatrix_[a][b] = *_iter;
        }
        b = 0;
    }
    Matrix m1(myMatrix_, 3, 3);
    m1.display();
    Matrix M(3,3);
    m1.transposeMatrix(m1, M);
    M.display();
    m1.pseudoInverse(m1).display();
    m1.setSubMatrix(3, 0, 0).pseudoInverse().display();

    std::cout << '\n'
              << "Rectangular matrix: " << '\n'
              << "sub matrix is : " << '\n';
    m1.printSubMatrix(m1, 0, 0, 1, 3);
    std::cout << "sub matrix is : " << '\n';
    m1.setSubMatrix(0, 0, 1, 2).display();
    std::cout << "sub matrix is : " << '\n';
    m1.setSubMatrix(0, 0, 3, 3).display();
    std::cout << "sub matrix is : " << '\n';
    m1.setSubMatrix(0, 0, 2, 3).display();
    std::cout << "sub matrix is : " << '\n';
    m1.setSubMatrix(0, 0, 1, 3).display();

    std::cout << '\n'
              << "Square matrix: " << '\n'
              << "sub matrix is : " << '\n';
    std::cout << "sub matrix is : " << '\n';
    m1.setSubMatrix(1, 0, 0).display();
    std::cout << "sub matrix is : " << '\n';
    m1.setSubMatrix(2, 0, 0).display();
    std::cout << "sub matrix is : " << '\n';
    m1.setSubMatrix(3, 0, 0).display();
    */
    /*
    Matrix M = {{1, 4, 8}, {-2, 9, 0}, {2, -13, 5.5}, {6.5, -1, 55}};
    Matrix M_T(3, 4);
    Matrix::transposeMatrix(M, M_T);

    cout << "M:" << "\n"
         << M << "\n" << "\n"
         << "M_T:" << "\n"
         << M_T << "\n" << endl;
    */
    /*
    Matrix M = { {5, -2, 2, 7},
                    {1, 0, 0, 3},
                    {-3, 1, 5, 0},
                    {3, -1, -9, 4}};

    Matrix A(4,4);


    std::cout << "Determinant of Matrix: "
              <<  M.determinantMatix() << '\n';

    M.adjointMatrix(A);
    std::cout << "Adjoint of matrix: " << '\n'
              << A << '\n';

    M.adjointMatrix();

    M.matrixInverse(M, A);
    std::cout << "Inverse of Matrix: " << '\n'
              << A << '\n';

    M.matrixInverse();

    std::cout << "Invertible or not: " << '\n'
              << M.isInvertible() << '\n';

    std::cout << "Matrix rank is: " << '\n'
              << M.rankMatrix() << '\n';
    */
    /*
    double** myMatrix_ = new double*[2];
    for(int c = 0; c < 2; ++c)
        myMatrix_[c] = new double[3];

    std::initializer_list<std::initializer_list<double> > initList { {10.0, 20.0 , 10.0},{4.0, 5.0, 6.0} };

    int a = 0, b = 0;
    for(auto iter = initList.begin(); iter < initList.end(); ++iter, ++a)
    {
        for(auto _iter = iter->begin(); _iter < iter->end(); ++_iter, ++b)
        {
            myMatrix_[a][b] = *_iter;
        }
        b = 0;
    }
    Matrix m1(myMatrix_, 2, 3);
    m1.display();
    m1.pseudoInverse(m1).display();
    */
    /*
    double** myMatrix_ = new double*[3];
    for(int c = 0; c < 3; ++c)
        myMatrix_[c] = new double[4];

    std::initializer_list<std::initializer_list<double> > initList {{12000.0, 23999.0, 34.0, 33.0},{12.0, 45.0, 54.0, 15.0},{9.0, 88, 10.0, 14.0}};

    int a = 0, b = 0;
    for(auto iter = initList.begin(); iter < initList.end(); ++iter, ++a)
    {
        for(auto _iter = iter->begin(); _iter < iter->end(); ++_iter, ++b)
        {
            myMatrix_[a][b] = *_iter;
        }
        b = 0;
    }
    Matrix m1(myMatrix_, 3, 4);
    Vector v1(3);

    std::cout << "Output method of back-substitution have passing argument is m1 matrix class: " << '\n';
    m1.backSubstitution(m1);
    std::cout << "Output method of back-substitution have no passing argument is m1 matrix class: " << '\n';
    m1.backSubstitution().display();
    std::cout << "Output of linear system with no passing argument: " << '\n';
    m1.gaussianElimination().display();
    std::cout << "Output of linear system with passing argument: " << '\n';
    m1.gaussianElimination(m1);
    */
    /*
    Matrix m2{{1,2,3,4},{5,6,7,8},{9,10,11,12}};
    Matrix m3(3, 8);
    m3 = m1.makeAugementedMatrix(m1, 3, 4, m2, 4);
    std::cout << "This is augmented matrix m3: " << '\n'
              << m3 << '\n';
    */
    /*
    Matrix mat4{{3,12,32,4,111},
                {2,32,12,5,11},
                {9,8,7,6,5},
                {11,12,13,14,15},
                {12,32,24,43,56}};
    Vector v4{23,13,4,19,76};
    */
    /*+++++++++++++++++++++++++++++++ Test LinearSystem Class +++++++++++++++++++++++++++++++++++*/
    //LinearSystem lin(mat4, v4);
    //Vector v7 = lin.Solve();
    /*
    std::cout << "v7 is \n"
              << v7; */
    /*
   // Not able to differentiate the comma and bracket !!!
    LinearSystem linear {{{3,12,32,4,111},
                        {2,32,12,5,11},
                        {9,8,7,6,5},
                        {11,12,13,14,15},
                        {12,32,24,43,56}},
                        {23,13,4,19,76}};

    std::cout << "Output for constructor object linear by initializer list: " << '\n';
    //linear.printVec_linear(linear); // Need testing
    auto v2 = linear.SolveAll();
    std::cout << "v2 is " << '\n'
              << v2 << '\n';
    */
    /*
    LinearSystem linear { {{13, 15, 16},
                          {1, 3, 5},
                          {21, 12, 24}},
                          {12, 13, 14}  };


    auto v2 = linear.Solve();
    std::cout << "v2: " << '\n'
              << v2 << '\n';
    */
/*
    Matrix mat4{{3,12,32,4,12},
                {12,32,12,12,11},
                {32,12,7,13,24},
                {4,12,13,14,43},
                {12,11,24,43,56}};

    std::cout << "Output for matrix A: " << '\n'
              << mat4 << '\n';

    Vector v4{23,13,4,19,76}; */
/*
    Matrix temp(5,5);
    mat4.matrixInverse(mat4,temp);
    std::cout << "Output matrix inverse of A: " << '\n'
              << temp << '\n';

    Matrix mat1{{1,0,2,-1},{3,0,0,5},{2,1,4,-3},{1,0,5,0}};
    Vector v1{1,2,3,4};

    Matrix mat2{{1,2},{3,4}};
    Vector v2{5,6};

    std::cout << "Output for determinant of matrix mat1: " << '\n'
              << mat1.determinantMatix(mat4, 5) << '\n';


    double d[5][1]{23,13,4,19,76};
    std::cout << "Vector b but construct by Matrix class : " << '\n';
    Matrix b(d,5,1);
*/
    /*
    for(int i = 1; i <= 5; ++i)
        std::cin >> b(i,1);
    */
    /*
    std::cout << "Output for matrix b: " << '\n'
              << b << '\n';
    std::cout << "Enter the starting point" << '\n';
    double d1[5][1]{-23,-13,-4,-19,-76};
    Matrix x(d1, 5,1);*/
    /*
    Matrix x(5,1);
    for(int i = 1; i <= 5; ++i)
        cin >> x(i, 1);
    std::cout << "Output for starting point: " << '\n'
              << x << '\n';*/
    /*
    unsigned maxIter;
    std::cout << "Enter the max number of iterations" << '\n';
    std::cin >> maxIter;

    double tol;
    std::cout << "Enter tolerance: " << '\n';
    std::cin >> tol;


    // the matrix should be a symmetric matrix
    PosSymLinSystem p1;
    p1.ConjugateGradientMethodLinearMatrix(mat4, b, x, tol, maxIter);

    PosSymLinSystem p2(mat4, v4);
    p2.Solve();
    */
/*
    PosSymLinSystem p3(mat4,v4);
    p3.SolveAll();
*/
    /*
    // Not seem to be correct so do not use this method
    PosSymLinSystem p3;
    Matrix mat3 {{14,1,-13},{1,22,6},{-13,6,17}};
    std::cout << "Output for cholesky: " << '\n' << p3.CholeskyDecomposition(mat3) << '\n';*/
    /*
    // Some new implementation about overloading operation
    Matrix mat3 {{1,2,3},{4,22,35},{-1,6,17}};
    Matrix tMatrix(3,3);
    tMatrix = mat3^'T';
    tMatrix.display();
    Matrix inv(3,3);
    inv = mat3^-1;
    inv.display();
    */
    /*
    Matrix mat5{{3,12,32,4,12,   12,4,32,12},
                {12,32,12,12,11, 11,12,12,32},
                {32,12,7,13,24,  24,13,7,12},
                {4,12,13,14,43,  43,14,13,12},
                {12,11,24,43,56, 56,43,24,11},

                {1,2,3,4,5,      5,4,3,2},
                {23,24,25,26,27, 27,26,25,24},
                {32,43,45,65,45, 45,65,45,43},
                {9,8,7,6,5,      5,6,7,8},
                {21,32,43,54,65, 65,54,43,32}};

    Matrix mat(10,9);*/
    /*Matrix mat4{{3,12,32,4,12},
                {12,32,12,12,11},
                {32,12,7,13,24},
                {4,12,13,14,43},
                {12,11,24,43,56}};
    mat4.geninv().display();*/
    //mat = mat5.pseudoInverse();
    /*
    Matrix mat2{{21,311,41,25,63,7},
                {12,133,14,15,16,17},
                {22,23,24,25,261,27},
                {32,2,34,35,36,37},
                {423,43,44,452,46,4},
                {521,53,54,515,56,57}};

    std::cout << "Determinant of matrix: " << mat2.determinantMatix(mat2) << '\n';
    */
    //auto x = mat2.determinantMatix();
    //Matrix m;
    //std::cout << "\nOutput for determinant of matrix: " << '\n' << x;
    //mat2.adjointMatrix(mat2, m);

    //std::cout << "\nOutput for det of matrix: " << '\n' << mat2.determinantMatix(mat2);
    /*
    mat2.matrixInverse();

    mat5.geninv().display();
    LinearSystem ls;*/
    /*
    std::cout << "Output for pseudo-inverse matrix: " << '\n';
    mat.display();

    std::cout << "Output for matrix A: " << '\n'
              << mat5 << '\n';

    Vector v5{23,13,4,19,76, 76,19,4,13,23};
    LinearSystem lin(mat5, v5);
    */

    /*Vector sol(10);
    sol = lin.SolveAll();*/

    Dataset data;
    data.excecuteFile();
    //data.ExecuteFile(); still I have not come with this method yet....
    return 0;
}
