#pragma once
// example header file
// this is where the signatures of functions and classes are stored.

double example_function(double);

class example_class
{
    // Keep everything public for now -- or equivalently use struct in place of class
    public:
        double example_member;

        example_class(double); // constructor
        double example_member_function(double);
};
