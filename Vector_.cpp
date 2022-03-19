#include "Vector_.hpp"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <initializer_list>
#include <new>
#include <iomanip>

Vector::Vector() noexcept
    : mSize(0)
    , mData(new long double[mSize])
{
    mData = nullptr;
    //std::cout << "Default constructor Vector class.." << '\n';
}

Vector::Vector(unsigned int n)
{
    try
    {
        if(n == 0)
            throw std::out_of_range("error: Value is greater than 0");
    }
    catch(std::out_of_range& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }

    mSize = n;
    mData = new long double[mSize];
    for(unsigned int i = 0; i < mSize; ++i)
        mData[i] = 0.0;
}

Vector::Vector(const Vector& v)
    : mSize(v.mSize)
    , mData(nullptr)
{
    //std::cout << "The size of member data = "
              //<< v.mSize << '\n' << "Copying constructor." << '\n';
    try
    {
        if(mSize == 0)
        {
            throw  "Zero size of vector!" ;
        }
    }
    catch(const std::string msg)
    {
        std::cout << msg << '\n';
        return;
    }
    mData = new long double[mSize];
    //mData = (long double*)malloc(mSize*sizeof(long double));
    for(unsigned int i = 0; i < mSize; ++i)
        mData[i] = v.mData[i];
}

Vector::Vector(long double* m, unsigned int n) noexcept
{
    mSize = n;
    mData =  new long double[mSize];

    for(unsigned int i = 0; i < mSize; ++i)
        mData[i] = m[i];

}

Vector::Vector(std::initializer_list<long double> init)
    : mData(nullptr)
{
    mSize = init.size();

    mData = new long double[mSize];

    int i = 0;
    for(auto iter = init.begin(); iter < init.end(); ++iter, ++i)
    {
        mData[i] = *iter;
    }
}

Vector::~Vector()
{
    //std::cout << "Destructor Vector class.." << '\n';
    if(mData)
    {
        //std::cout << "Deleting member data(mData)!" << '\n';
        delete [] mData;
    }
    //std::cout << "Successfully deleting mData!" << '\n';
}

Vector& Vector::operator=(Vector const& v) noexcept
{
    //std::cout << "Assignment copy operator Vector class.." << '\n';

    if(this != &v)
    {
        if(mSize != v.mSize)
        {
            delete [] mData;
            mSize = v.mSize;
            mData = new long double[mSize];
        }
    }
    for(unsigned int i = 0; i < mSize; ++i)
        mData[i] = v.mData[i];

    return (*this);
}

Vector::Vector(Vector&& v) noexcept
    : mSize(0)
    , mData(nullptr)
{
    //std::cout << "Move constructor Vector class.." << '\n';
    mSize = v.mSize;

    for(unsigned int i = 0; i < mSize; ++i)
        mData[i] = v.mData[i];
    // This prevents the destructor from freeing memory multiple times
    for(unsigned int i = 0; i < mSize; ++i)
        v.mData[i] = 0.0;

    // Using above code or this below code
    v.mData = nullptr;

    v.mSize = 0;
}

Vector& Vector::operator=(Vector&& v) noexcept
{
    //std::cout << "Move constructor assignment Vector class.." << '\n';

    if(this != &v)
    {
        if(mSize != v.mSize)
        {
            delete [] mData;
            mSize = v.mSize;
            mData = new long double[mSize];
        }
    }
    for(unsigned int i = 0; i < mSize; ++i)
        mData[i] = v.mData[i];

    for(unsigned int i = 0; i < mSize; ++i)
        v.mData[i] = 0.0;
    v.mData = nullptr;

    v.mSize = 0;

    return (*this);
}

int Vector::length(const Vector& v) {return v.mSize ;}

const unsigned int Vector::Size() { return  mSize; }

std::ostream& operator<<(std::ostream& os, const Vector& obj)
{
    for(unsigned int i = 0; i < obj.mSize; ++i)
        os << obj.mData[i] << " ";
    return os;
}

std::istream& operator>>(std::istream& is, Vector& obj)
{
    for(unsigned int i = 0; i < obj.mSize; ++i)
        is >> obj.mData[i];
    return is;
}

bool Vector::operator[](unsigned int idx)
{
    try
    {
        if(idx < 0 || idx > mSize)
        {
            throw std::out_of_range("error: out of range");
        }
    }
    catch(std::out_of_range& er)
    {
        std::cout << er.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }

    return true;
}

bool Vector::operator[](unsigned int idx) const
{
    try
    {
        if(idx < 0 || idx > mSize)
        {
            throw std::out_of_range("error: out of range");
        }
    }
    catch(std::out_of_range& er)
    {
        std::cout << er.what() << '\n';
        return false;
        exit(EXIT_FAILURE);
    }

    return true;
}

long double& Vector::operator()(unsigned int rhs) const
{
    try
    {
        if(rhs > mSize)
        {
            throw std::out_of_range("Error: out of range");
        }
    }
    catch(std::out_of_range& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }

    return mData[rhs - 1];
}

long double& Vector::operator()(unsigned int rhs)
{
    try
    {
        if(rhs > mSize)
        {
            throw std::out_of_range("Error: out of range");
        }
    }
    catch(std::out_of_range& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }

    return mData[rhs - 1];
}

/* Binary operator */
// Arithmetic operator
const Vector Vector::operator*(const long double& rhs) const
{
    Vector res(mSize);

    for(unsigned int i = 0; i < res.mSize; ++i)
        res.mData[i] = mData[i] * rhs;

    return res;
}

Vector operator*(const long double& lhs, const Vector& rhs)
{
    Vector res(rhs.mSize);


    for(unsigned int i = 0; i < res.mSize; ++i)
        res.mData[i] = lhs * rhs.mData[i];

    return res;
}

const Vector Vector::operator*(const long double& rhs)
{
    Vector res(mSize);

    for(unsigned int i = 0; i < res.mSize; ++i)
        res.mData[i] = mData[i] * rhs;

    return res;
}

const Vector Vector::operator+(const Vector& rhs) const
{
    Vector res;
    try
    {
        if(mSize != rhs.mSize)
            throw std::logic_error("Error: the size of two vectors is not appropriate!");
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }
    if(mSize == rhs.mSize)
    {
        for(unsigned int i = 0; i < res.mSize; ++i)
            res.mData[i] = mData[i] + rhs.mData[i];

    }
    return res;
}

const Vector Vector::operator+(const Vector& rhs)
{
    Vector res(mSize);
    try
    {
        if(mSize != rhs.mSize)
            throw std::logic_error("Error: the size of two vectors is not appropriate!");
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }
    if(mSize == rhs.mSize)
    {
        for(unsigned int i = 0; i < res.mSize; ++i)
            res.mData[i] = mData[i] + rhs.mData[i];

    }
    return res;
}

const Vector Vector::operator-(const Vector& rhs) const
{
    Vector res;
    try
    {
        if(mSize != rhs.mSize)
            throw std::logic_error("Error: the size of two vectors is not appropriate!");
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }
    if(mSize == rhs.mSize)
    {
        for(unsigned int i = 0; i < res.mSize; ++i)
            res.mData[i] = mData[i] - rhs.mData[i];

    }
    return res;
}

const Vector Vector::operator-(const Vector& rhs)
{
    Vector res(mSize);
    try
    {
        if(mSize != rhs.mSize)
            throw std::logic_error("Error: the size of two vectors is not appropriate!");
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }
    if(mSize == rhs.mSize)
    {
        res.mSize = mSize;

        for(unsigned int i = 0; i < res.mSize; ++i)
            res.mData[i] = mData[i] - rhs.mData[i];

    }
    return res;
}

// Multiplication between 2 vectors underlying dot product following the constant Vector class
long double Vector::operator*(const Vector& rhs) const
{
    Vector res(mSize);
    long double sum = 0.0;

    for(unsigned int i = 0; i < res.mSize; ++i)
    {
        res.mData[i] = mData[i] * rhs.mData[i];
        sum += res.mData[i];
    }

    return sum;
}

// Multiplication between 2 vectors underlying dot product following the Vector class
long double Vector::operator*(const Vector& rhs)
{
    Vector res(mSize);
    long double sum = 0.0;

    for(unsigned int i = 0; i < res.mSize; ++i)
    {
        res.mData[i] = mData[i] * rhs.mData[i];
        sum += res.mData[i];
    }

    return sum;
}

// Assignment operator
Vector& Vector::operator+=(const Vector& rhs)
{
    try
    {
        if(mSize != rhs.mSize)
            throw std::logic_error("Error: the size of two vectors is not appropriate!");
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }
    if(mSize == rhs.mSize)
    {
        for(unsigned int i = 0; i < mSize; ++i)
            this->mData[i] += rhs.mData[i];

    }
    return (*this);
}

Vector& Vector::operator-=(const Vector& rhs)
{
    try
    {
        if(mSize != rhs.mSize)
            throw std::logic_error("Error: the size of two vectors is not appropriate!");
        else
        {
            for(unsigned int i = 0; i < mSize; ++i)
                this->mData[i] -= rhs.mData[i];
        }
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        exit(EXIT_FAILURE);
    }
    return (*this);
}

long double Vector::operator*=(const Vector& rhs)
{
    mSize = rhs.mSize;
    long double sum = 0.0;

    for(unsigned int i = 0; i < mSize; ++i)
    {
        this->mData[i] *= rhs.mData[i];
        sum += this->mData[i];
    }

    return sum;
}

Vector& Vector::operator*=(const long double& rhs)
{
    for(unsigned int i = 0; i < mSize; ++i)
        this->mData[i] += mData[i] * rhs;

    return (*this);
}

// Logical operator
bool Vector::operator==(const Vector& rhs)
{
    try
    {
        if(this->mSize != rhs.mSize)
        {
            throw std::logic_error("Error: The 2 vectors have not share the same size!");
        }
        else if(this->mSize == rhs.mSize)
        {
            for(unsigned int i = 0; i < mSize; ++i)
            {
                if(this->mData[i] == rhs.mData[i]) return true;
                else return false;
            }
        }
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << "\n";
        return false;
    }
    return true;
}

bool Vector::operator==(const Vector& rhs) const
{
    try
    {
        if(this->mSize != rhs.mSize)
        {
            throw std::logic_error("Error: The 2 vectors have not share the same size!");
        }
        else if(this->mSize == rhs.mSize)
        {
            for(unsigned int i = 0; i < mSize; ++i)
            {
                if(this->mData[i] == rhs.mData[i]) continue;
                else return false;
            }
        }
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << "\n";
        return false;
    }
    return true;
}

bool Vector::operator!=(const Vector& rhs)
{
    try
    {
        if(this->mSize != rhs.mSize)
            throw std::length_error("Error: The 2 vectors have not share the same size!");
        else if(this->mSize == rhs.mSize)
        {
            for(unsigned int i = 0; i < mSize; ++i)
            {
                if(this->mData[i] != rhs.mData[i]) continue;
                else return false;
            }
        }
    }
    catch(std::length_error& err)
    {
        std::cout << err.what() << '\n';
        return false;
    }
    return true;
}

bool Vector::operator!=(const Vector& rhs) const
{
    try
    {
        if(this->mSize != rhs.mSize)
            throw std::logic_error("Error: The 2 vectors have not share the same size!");
        else if(this->mSize == rhs.mSize)
        {
            for(unsigned int i = 0; i < mSize; ++i)
            {
                if(this->mData[i] != rhs.mData[i]) continue;
                else return false;
            }
        }
    }
    catch(std::logic_error& err)
    {
        std::cout << err.what() << '\n';
        return false;
    }
    return true;
}

bool Vector::operator<=(const Vector& rhs)
{
    try
    {
        if(this->mSize != rhs.mSize)
            throw std::logic_error("Error: The 2 vectors have not share the same size!");
        else if(this->mSize == rhs.mSize)
        {
            for(unsigned int i = 0; i < mSize; ++i)
            {
                if(this->mData[i] <= rhs.mData[i]) return true;
                else return false;
            }
        }
    }
    catch(std::logic_error& err)
    {
        std::cout << err.what() << '\n';
        return false;
    }
    return true;
}

bool Vector::operator<=(const Vector& rhs) const
{
    try
    {
        if(this->mSize != rhs.mSize)
            throw std::logic_error("Error: The 2 vectors have not share the same size!");
        else if(this->mSize == rhs.mSize)
        {
            for(unsigned int i = 0; i < mSize; ++i)
            {
                if(this->mData[i] <= rhs.mData[i]) continue;
                    else return false;
            }
        }
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        return false;
    }
    return true;
}

bool Vector::operator>=(const Vector& rhs)
{
    try
    {
        if(this->mSize != rhs.mSize)
            throw std::logic_error("Error: The 2 vectors have not share the same size!");
        else if(this->mSize == rhs.mSize)
        {
            for(unsigned int i = 0; i < mSize; ++i)
            {
                if(this->mData[i] >= rhs.mData[i]) return true;
                    else return false;
            }
        }
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        return false;
    }
    return true;
}

bool Vector::operator>=(const Vector& rhs) const
{
    try
    {
        if(this->mSize != rhs.mSize)
            throw std::logic_error("Error: The 2 vectors have not share the same size!");
        else if(this->mSize == rhs.mSize)
        {
            for(unsigned int i = 0; i < mSize; ++i)
            {
                if(this->mData[i] >= rhs.mData[i]) return true;
                    else return false;
            }
        }
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        return false;
    }
    return true;
}

bool Vector::operator>(const Vector& rhs)
{
    try
    {
        if(this->mSize != rhs.mSize)
            throw std::logic_error("Error: The 2 vectors have not share the same size!");
        else if(this->mSize == rhs.mSize)
        {
            for(unsigned int i = 0; i < mSize; ++i)
            {
                if(this->mData[i] > rhs.mData[i]) continue;
                    else return false;
            }
        }
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        return false;
    }
    return true;
}

bool Vector::operator>(const Vector& rhs) const
{
    try
    {
        if(this->mSize != rhs.mSize)
            throw std::logic_error("Error: The 2 vectors have not share the same size!");
        else if(this->mSize == rhs.mSize)
        {
            for(unsigned int i = 0; i < mSize; ++i)
            {
                if(this->mData[i] > rhs.mData[i]) continue;
                    else return false;
            }
        }
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        return false;
    }
    return true;
}

bool Vector::operator<(const Vector& rhs)
{
    try
    {
        if(this->mSize != rhs.mSize)
            throw std::logic_error("Error: The 2 vectors have not share the same size!");
        else if(this->mSize == rhs.mSize)
        {
            for(unsigned int i = 0; i < mSize; ++i)
            {
                if(this->mData[i] < rhs.mData[i]) continue;
                    else return false;
            }
        }
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        return false;
    }
    return true;
}

bool Vector::operator<(const Vector& rhs) const
{
    try
    {
        if(this->mSize != rhs.mSize)
            throw std::logic_error("Error: The 2 vectors have not share the same size!");
        else if(this->mSize == rhs.mSize)
        {
            for(unsigned int i = 0; i < mSize; ++i)
            {
                if(this->mData[i] < rhs.mData[i]) continue;
                    else return false;
            }
        }
    }
    catch(std::logic_error& er)
    {
        std::cout << er.what() << '\n';
        return false;
    }
    return true;
}

void Vector::display() const
{
    for(std::size_t i = 0; i < this->mSize; ++i)
        std::cout << std::setw(10) << this->mData[i] << " " << '\n';
    std::cout << std::endl;
}
