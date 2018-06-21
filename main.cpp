#include <iostream>
#include <Eigen/Dense>
#include "src/example.cpp"

int main()
{
    std::cout << "Hello World" << std::endl;

    // Using eigen
    Eigen::MatrixXd M(2,2); // declare a 2x2 matrix of doubles
    Eigen::VectorXd x(2); // declare a 2d vector of doubles
    M(0,0) = 1;
    M(0,1) = 2;
    M(1,0) = 3;
    M(1,1) = 4;
    x(0) = 5;
    x(1) = 6;

    std::cout << "M:" << std::endl;
    std::cout << M << std::endl;

    std::cout << "x^T:" << std::endl;
    std::cout << x.transpose() << std::endl;

    std::cout << "inverse of M:" << std::endl;
    std::cout << M.inverse() << std::endl;

    std::cout << "Mx:" << std::endl;
    std::cout << M*x << std::endl;

    // Using the class found in src/example.hpp
    example_class X{5};
    std::cout << "example_class usage:" << std::endl;
    std::cout << X.example_member << std::endl;
    std::cout << X.example_member_function(6) << std::endl;


}
