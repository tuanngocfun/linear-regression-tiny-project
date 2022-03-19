#ifndef VECTOR_CLASS
#define VECTOR_CLASS
#include <iostream>
#include <initializer_list>

// forward declaration for class Matrix to multiply Matrix obj * Vector obj
class Matrix;

class Vector
{
public:
    Vector() noexcept;
    explicit Vector(unsigned int n);
    Vector(Vector const& v) ;

    explicit Vector(long double* m, unsigned int x) noexcept;

    Vector(std::initializer_list<long double> init);
    // Assignment copy constructor
    Vector& operator=(Vector const& v) noexcept;

    virtual ~Vector();

    // Move constructor
    Vector(Vector&& v) noexcept;

    // Assignment move constructor
    Vector& operator=(Vector&& v) noexcept;

    // Getter
    int getSize() const {return mSize;}
    long double* getData() const {return mData;}

    // Setter
    void setSize(unsigned int x) {mSize = x;}
    void setData(long double *y) {mData = y;}

    // Size of an array
    const unsigned int Size();
    int length(const Vector &v);

    bool operator[](unsigned int idx);
    bool operator[](unsigned int idx) const;

    long double& operator()(unsigned int rhs);
    long double& operator()(unsigned int rhs) const;

    friend Vector operator+(const long double& lhs, const Vector& rhs) = delete;
    friend Vector operator-(const long double& lhs, const Vector& rhs) = delete;
    friend Vector operator*(const long double& lhs, const Vector& rhs);
    friend Vector operator/(const long double& lhs, const Vector& rhs) = delete;
    friend Vector operator%(const long double& lhs, const Vector& rhs) = delete;

    const Vector operator+(const long double& rhs) const = delete;
    const Vector operator-(const long double& rhs) const = delete;
    const Vector operator*(const long double& rhs) const;
    const Vector operator/(const long double& rhs) const = delete;
    const Vector operator%(const long double& rhs) const = delete;

    const Vector operator+(const long double& rhs) = delete;
    const Vector operator-(const long double& rhs) = delete;
    const Vector operator*(const long double& rhs);
    const Vector operator/(const long double& rhs) = delete;
    const Vector operator%(const long double& rhs) = delete;

    const Vector operator+(const Vector& rhs) const;
    const Vector operator-(const Vector& rhs) const;
    long double operator*(const Vector& rhs) const;
    const Vector operator/(const Vector& rhs) const = delete;
    const Vector operator%(const Vector& rhs) const = delete;

    const Vector operator+(const Vector& rhs);
    const Vector operator-(const Vector& rhs);
    long double operator*(const Vector& rhs);
    const Vector operator/(const Vector& rhs) = delete;
    const Vector operator%(const Vector& rhs) = delete;

    // Matrix multiply Vector and return Vector type
    friend Vector operator*(const Matrix& mat, const Vector& vec);

    friend std::ostream& operator<<(std::ostream& os, Vector const& obj);
    friend std::istream& operator>>(std::istream& is, Vector& obj);

    Vector& operator+=(const long double& rhs) = delete;
    Vector& operator-=(const long double& rhs) = delete;
    Vector& operator*=(const long double& rhs);
    Vector& operator/=(const long double& rhs) = delete;
    Vector& operator%=(const long double& rhs) = delete;

    Vector& operator+=(const Vector& rhs);
    Vector& operator-=(const Vector& rhs);
    long double operator*=(const Vector& rhs);
    Vector& operator/=(const Vector& rhs) = delete;
    Vector& operator%=(const Vector& rhs) = delete;

    /* Unary operator */
    Vector operator++() = delete;   // ++c
    Vector operator++(int) = delete; // c++
    Vector operator--() = delete;
    Vector operator--(int) = delete;

    const Vector operator++() const = delete;
    const Vector operator--() const = delete;
    Vector operator++(int) const = delete;
    Vector operator--(int) const = delete;
    /* end unary operator */

    bool operator==(const Vector& rhs);
    bool operator!=(const Vector& rhs);
    bool operator>(const Vector& rhs);
    bool operator<(const Vector& rhs);
    bool operator>=(const Vector& rhs);
    bool operator<=(const Vector& rhs);

    bool operator==(const Vector& rhs) const;
    bool operator!=(const Vector& rhs) const;
    bool operator>(const Vector& rhs) const;
    bool operator<(const Vector& rhs) const;
    bool operator>=(const Vector& rhs) const;
    bool operator<=(const Vector& rhs) const;

    friend bool operator==(const long double& lhs, const Vector& rhs) = delete;
    friend bool operator!=(const long double& lhs, const Vector& rhs) = delete;
    friend bool operator>(const long double& lhs, const Vector& rhs) = delete;
    friend bool operator<(const long double& lhs, const Vector& rhs) = delete;
    friend bool operator>=(const long double& lhs, const Vector& rhs) = delete;
    friend bool operator<=(const long double& lhs, const Vector& rhs) = delete;

    bool operator==(const long double& rhs) = delete;
    bool operator!=(const long double& rhs) = delete;
    bool operator>(const long double& rhs) = delete;
    bool operator<(const long double& rhs) = delete;
    bool operator>=(const long double& rhs) = delete;
    bool operator<=(const long double& rhs) = delete;

    bool operator==(const long double& rhs) const = delete;
    bool operator!=(const long double& rhs) const = delete;
    bool operator>(const long double& rhs) const = delete;
    bool operator<(const long double& rhs) const = delete;
    bool operator>=(const long double& rhs) const = delete;
    bool operator<=(const long double& rhs) const = delete;

    // if need any of them then omit delete and define the method later
    bool operator&&(const Vector& rhs) = delete;
    bool operator||(const Vector& rhs) = delete;
    bool operator!() = delete;

    bool operator&&(const Vector& rhs) const = delete;
    bool operator||(const Vector& rhs) const = delete;
    bool operator!() const = delete;

    // It will be lost the precision when try to convert double to int
    Vector operator>>(const Vector& rhs)  = delete;
    Vector operator<<(const Vector& rhs) = delete;
    Vector operator^(const Vector& rhs) = delete;
    Vector operator|(const Vector& rhs) = delete;
    Vector operator&(const Vector& rhs) = delete;
    Vector operator~();

    Vector operator>>(const Vector& rhs) const = delete;
    Vector operator<<(const Vector& rhs) const = delete;
    Vector operator^(const Vector& rhs) const = delete;
    Vector operator|(const Vector& rhs) const = delete;
    Vector operator&(const Vector& rhs) const = delete;
    Vector operator~() const;

    Vector operator>>(const int& rhs) = delete;
    Vector operator<<(const int& rhs) = delete;
    Vector operator^(const int& rhs) = delete;
    Vector operator|(const int& rhs) = delete;
    Vector operator&(const int& rhs) = delete;

    Vector operator>>(const int& rhs) const = delete;
    Vector operator<<(const int& rhs) const = delete;
    Vector operator^(const int& rhs) const = delete;
    Vector operator|(const int& rhs) const = delete;
    Vector operator&(const int& rhs) const = delete;

    virtual void display() const;

protected:
    unsigned int mSize;
    long double *mData;

    // Now LinearSystem class can access the private data member of Vector class because parent class can not access the child class (inheritance case)
    friend class LinearSystem;
    // Now Matrix class can access the private data member of Vector class because parent class can not access the child class (inheritance case)
    friend class Matrix;
    // Now PosSymLinSystem class can access the private data member of Vector class because parent class can not access the child class (inheritance case)
    friend class PosSymLinSystem;
    // Now Dataset class can access the private data member of Vector class because parent class can not access the child class (inheritance case)
    friend class Dataset;
};
#endif // VECTOR_CLASS
