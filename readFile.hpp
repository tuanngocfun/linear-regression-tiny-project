#ifndef READFILE_HPP
#define READFILE_HPP
#include "Vector_.hpp"
#include "Matrix_.hpp"
#include "LinearSystem.hpp"
#include "PosSymLinSystem.hpp"
#include <fstream>  // std::ifstream
#include <sstream>  // std::istringstream
#include <cctype>   // std::toupper
#include <string>   // std::stoi
#include <ctime>    // std::time, std::localtime
#include <stdlib.h> // defines putenv in POSIX
#include <algorithm>// std::for_each
#include <map>
#include <vector>
#include <iomanip>  // std::put_time
#include <chrono>   // std::chrono
#include <thread>   // std::this_thread
#include <iostream>

/*
struct TIME
{
    int tm_sec;
    int tm_min;
    int tm_hour;
    int tm_mday;
    int tm_month;
    int tm_year;
    int tm_wday;    // since weekday
};
*/
using namespace std::literals::chrono_literals;

class Dataset
{
public:
    Dataset(){}
    virtual ~Dataset(){}
    // 7 setters and 7 getters for numeric data
    void setMYCT(int value);
    const int& getMYCT() const;
    void setMMIN(int value);
    const int& getMMIN() const;
    void setMMAX(int value);
    const int& getMMAX() const;
    void setCACH(int value);
    const int& getCACH() const;
    void setCHMIN(int value);
    const int& getCHMIN() const;
    void setCHMAX(int value);
    const int& getCHMAX() const;
    void setPRP(int value);
    const int& getPRP() const;
    // 2 pair of setters and getters for information of the product's name
    void setVendorName(std::string vendor);
    std::string getVendorName() const;
    void setModelName(std::string model);
    std::string getModelName() const;
    // A pair of setter and getter for testing set data
    void setERP(int value);
    int getERP() const;
    // Counting number of row and column of the linear system
    const std::size_t getRowOfData() const;
    const std::size_t getColOfData() const;
    // Method to execute file
    void excecuteFile();
    void ExecuteFile();

    friend std::ostream& operator<<(std::ostream& out, const Dataset& d);
    friend std::istream& operator>>(std::istream& is, Dataset& d);

private:
    std::string _VendorName;
    std::string _ModelName;
    //int MYCT;           // 3. machine cycle time in nanoseconds (integer)
    //int MMIN;           // 4.minimum main memory in kilobytes (integer)
    //int MMAX;           // 5.maximum main memory in kilobytes (integer)
    //int CACH;           // 6.cache memory in kilobytes (integer)
    //int CHMIN;          // 7.minimum channels in units (integer)
    //int CHMAX;          // 8.maximum channels in units (integer)
    //int PRP;            // 9.published relative performance (integer)
    int ERP;              // 10. estimated relative performance from the original article (integer)
    std::map<int, int> dataValue;
};

// Dataset for coefficients of given linear system from text.file and solution vector
void Dataset::setMYCT(int value)
{
    dataValue.insert_or_assign(0, value);
}

const int& Dataset::getMYCT() const
{
    return dataValue.at(0);
}

void Dataset::setMMIN(int value)
{
    dataValue.insert_or_assign(1, value);
}

const int& Dataset::getMMIN() const
{
    return dataValue.at(1);
}

void Dataset::setMMAX(int value)
{
    dataValue.insert_or_assign(2, value);
}

const int& Dataset::getMMAX() const
{
    return dataValue.at(2);
}

void Dataset::setCACH(int value)
{
    dataValue.insert_or_assign(3, value);
}

const int& Dataset::getCACH() const
{
    return dataValue.at(3);
}

void Dataset::setCHMIN(int value)
{
    dataValue.insert_or_assign(4, value);
}

const int& Dataset::getCHMIN() const
{
    return dataValue.at(4);
}

void Dataset::setCHMAX(int value)
{
    dataValue.insert_or_assign(5, value);
}

const int& Dataset::getCHMAX() const
{
    return dataValue.at(5);
}

void Dataset::setPRP(int value)
{
    dataValue.insert_or_assign(6, value);
}

const int& Dataset::getPRP() const
{
    return dataValue.at(6);
}

// Information of product name
void Dataset::setVendorName(std::string value)
{
    _VendorName = value;
}

std::string Dataset::getVendorName() const
{
    return _VendorName;
}

void Dataset::setModelName(std::string value)
{
    _ModelName = value;
}

std::string Dataset::getModelName() const
{
    return _ModelName;
}

// testing set for evaluation
void Dataset::setERP(int value)
{
    ERP = value;
}

int Dataset::getERP() const
{
    return ERP;
}

const std::size_t Dataset::getColOfData() const
{
    return dataValue.size();
}

std::istream& operator>>(std::istream& is, Dataset& d)
{
    std::string myct, mmin, mmax, cach, chmin, chmax, prp, erp;
    if(std::getline(is, d._VendorName, ',') &&
       std::getline(is, d._ModelName, ',')  &&
       std::getline(is, myct, ',')  &&
       std::getline(is, mmin, ',')  &&
       std::getline(is, mmax, ',')  &&
       std::getline(is, cach, ',')  &&
       std::getline(is, chmin, ',') &&
       std::getline(is, chmax, ',') &&
       std::getline(is, prp, ',')   &&
       std::getline(is, erp, ','))
    {
        d.setMYCT(std::stoi(myct));
        d.setMMIN(std::stoi(mmin));
        d.setMMAX(std::stoi(mmax));
        d.setCACH(std::stoi(cach));
        d.setCHMIN(std::stoi(chmin));
        d.setCHMAX(std::stoi(chmax));
        d.setPRP(std::stoi(prp));
        d.setERP(std::stoi(erp));
    }
    return is;
}

std::ostream& operator<<(std::ostream& out, const Dataset& data)
{
    for(auto&& x : data.dataValue)
    {
        out << x.first << " -> " << x.second << "\n";
    }
    return out;
}

// Not want user to change the Linear System, also because copy constructor is now in private so in order to use that.
// I make it as reference so it wont be an error anymore, and since I use it as reference so I can still do copy constructor but with reference
void Dataset::excecuteFile()
{
    char rflag;         // Readiness flag
    // static LinearSystem dummy;      // Since we want our variable to be as global variable, but it will end in this program (end in run-time)

    std::cout << "\n\n";
    std::time_t result = std::time(nullptr);
    tm* local_tm = localtime(&result);
    std::cout << "Year: " << 1900 + local_tm->tm_year << '\n';
    std::cout << "Month: " << 1 + local_tm->tm_mon << '\n';
    std::cout << "Day: " << local_tm->tm_mday << '\n';
    std::cout << "Local time:   " << std::put_time(std::localtime(&result), "%c %Z") << '\n';
    std::cout << "UTC:          " << std::put_time(std::gmtime(&result), "%c %Z") << '\n';
    std::cout << "Everything ready? If yes, press y. " << '\n';
    std::cout << "Otherwise, press any other key. " << '\n';
    std::cin >> rflag;

    if( std::toupper(rflag) == 'Y' )
    {
        std::ios_base::sync_with_stdio(false);
        std::ifstream fileName;
        fileName.open("machine.data", std::ios_base::in);

        if(!fileName)
        {
            std::cerr << "Unable to open file machine.txt" << '\n';
            auto start = std::chrono::high_resolution_clock::now();
            std::this_thread::sleep_for(2000ms);
            auto end   = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> elapsed = end-start;
            std::cout << "Please wait for: " << elapsed.count() << "ms" <<'\n';
            exit(EXIT_FAILURE);
        }else{
            Dataset d;
            std::string value;
            std::vector<Dataset> dContainer;
            while(!fileName.eof())
            {
                if(fileName.good())
                {
                    //int i; // counting line/row of the data in the machine.txt file
                    while( std::getline(fileName, value, '\n') )
                    {
                        std::istringstream line(value);

                        // ignore the first data from the file
                        std::getline(line, value, ',');
                        d.setVendorName(value);

                        // ignore the second data from the file
                        std::getline(line, value, ',');
                        d.setModelName(value);

                        std::getline(line, value, ',');
                        d.setMYCT(std::stoi(value));

                        std::getline(line, value, ',');
                        d.setMMIN(std::stoi(value));

                        std::getline(line, value, ',');
                        d.setMMAX(std::stoi(value));

                        std::getline(line, value, ',');
                        d.setCACH(std::stoi(value));

                        std::getline(line, value, ',');
                        d.setCHMIN(std::stoi(value));

                        std::getline(line, value, ',');
                        d.setCHMAX(std::stoi(value));

                        std::getline(line, value, ',');
                        d.setPRP(std::stoi(value));

                        // temporarily object value now is for ERP
                        std::getline(line, value, ',');
                        d.setERP(std::stoi(value));

                        dContainer.push_back(d);
                    }
                }
            }
            fileName.close();
            /* Check for output of the class Dataset member
            for(std::vector<Dataset>::iterator it = dContainer.begin(); it != dContainer.end(); ++it)
            {
                std::cout << it->getVendorName() << ", " << it->getModelName() << ", "
                << it->getMYCT() << ", " << it->getMMIN() << ", " << it->getMMAX() <<
                it->getCACH() << ", " << it->getCHMIN() << ", " << it->getCHMAX() <<
                it->getPRP() << ", " << it->getERP() << std::endl;
            }*/
            /*
            // Lambda expression
            std::cout << "Output for dContainer: " << '\n';
            std::for_each(dContainer.begin(), dContainer.end(), [](Dataset i){ std::cout << i ;});*/
            /* //code works similar to lambda expression
            for(std::size_t i = 0; i < dContainer.size(); ++i)
                std::cout << "Output for dContainer: " << '\n' << dContainer[i] << '\n'; */

            //std::cout << "Output for dContainer.size(): " << '\n' << dContainer.size() << '\n';
            //std::cout << "Output for Dataset.size(): " << '\n' << d.getColOfData() << '\n';

            Matrix mat( static_cast<int>(dContainer.size()), static_cast<int>(d.getColOfData()) - 1 );
            Vector vec( static_cast<int>(dContainer.size()) );
/*
            std::cout << "Output for mat row: " << '\n' << mat.mNumRows << '\n';
            std::cout << "Output for mat column: " << '\n' << mat.mNumCols << '\n';
            std::cout << "Output for vector size: " << '\n' << vec.mSize << '\n';
*/
            for(int i = 1; i <= mat.mNumRows; ++i){
                for(int j = 1; ;)
                {
                    if(j < mat.mNumCols){
                        mat(i,j++) = dContainer[i-1].getMYCT();
                        mat(i,j++) = dContainer[i-1].getMMIN();
                        mat(i,j++) = dContainer[i-1].getMMAX();
                        mat(i,j++) = dContainer[i-1].getCACH();
                        mat(i,j++) = dContainer[i-1].getCHMIN();
                        mat(i,j++) = dContainer[i-1].getCHMAX();
                    }
                    else{
                        //std::cout << "Output for vector: " << '\n';
                        vec(i) = dContainer[i-1].getPRP();
                        break;
                    }
                }
            }

            //std::cout << std::endl;
            //std::cout << "Output for matrix: " << '\n' << mat << '\n';
            //std::cout << "Flag here: " << '\n';

            //Matrix _AMatrix( static_cast<int>(d.getColOfData()) - 1, static_cast<int>(dContainer.size()) );
            LinearSystem lin(mat, vec);
            //std::cout << "Output for Lin system constructor: " << '\n';
            //lin.display();

            //std::cout << "Output for posSymLin system constructor: " << '\n' << lin << '\n';
            Vector sol = lin.SolveAll();
            sol.display();
            // Calculate RMSE and deviation
            double _sum = 0.0;
            double _sum2 = 0.0;
            for(std::size_t i = 1; i <= 6; ++i)
            {
                _sum += std::pow( (sol(i) - vec(i)), 2);
                _sum2 += std::fabs(sol(i) - vec(i))/vec(i) * 100;
            }
            double RMSE = std::sqrt(_sum / vec.mSize);
            double averageDeviation = _sum2 / vec.mSize;
            std::cout << "RMSE: " << RMSE  << ", average deviation: "  << averageDeviation << std::endl;

            //std::cout << "Output for matrix: " << '\n' << mat << '\n';

        }
    }
}

#endif // READFILE_HPP
