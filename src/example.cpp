#include "example.hpp"


double example_function(double x)
{
    double y = x * 2;
    return x + y;
}

example_class::example_class(double x)
{
    example_member = x;
}

double example_class::example_member_function(double x)
{
    double z = example_function(x);
    return z + x;
}
