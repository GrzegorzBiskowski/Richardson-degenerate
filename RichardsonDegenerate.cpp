#include <iostream>
#include <complex>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>

#include "rg-degen-lib.h"

systemData rGSystem;
double g = 0.01;

void init_data();

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::complex<double>> rapidities = { -2.1, -1.9, 2.0};
    init_data();
    cEquations test = cEquations(rGSystem.lambda,rGSystem.occupation, rGSystem.degeneracies, rGSystem.energies, g);
    test.newton_raphson();
    cPolynomial test2 = cPolynomial(rGSystem.lambda, rGSystem.energies, rGSystem.degeneracies, rGSystem.occupation, g);
    test2.root_finder(rapidities);
    for (int i = 0; i < rapidities.size(); i++)
    {
        std::cout << rapidities[i] << '\t';
    }
    std::cout << '\n';
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> chronoDiff = end - start;
    std::cout << "Time elapsed: " << chronoDiff.count() << "s.\n";
}

void init_data()
{
    rGSystem.energies = { -2.0, 0.0, 2.0  };
    rGSystem.degeneracies = { 2, 2, 2};
    rGSystem.occupation = { 2, 0, 1 };
    rGSystem.init_lambda();
}