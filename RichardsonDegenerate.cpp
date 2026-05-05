#include <iostream>
#include <complex>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>

#include "rg-degen-lib.h"

systemData rGSystem;

void init_data();

int main()
{
    auto start = std::chrono::high_resolution_clock::now();

    init_data();
    cEquations test = cEquations(rGSystem.lambda,rGSystem.occupation, rGSystem.degeneracies, rGSystem.energies, 0.1);
    test.function_i(rGSystem.lambda, 7);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> chronoDiff = end - start;
    std::cout << "Time elapsed: " << chronoDiff.count() << "s.\n";
}

void init_data()
{
    rGSystem.energies = { -2.0, 0.0, 2.0 };
    rGSystem.degeneracies = { 4,2,3 };
    rGSystem.occupation = { 2,0,1 };
    rGSystem.init_lambda();
}