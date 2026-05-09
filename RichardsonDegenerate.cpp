#include <iostream>
#include <complex>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <chrono>

#include "rg-degen-lib.h"

systemData rGSystem;

std::ofstream fEnergiesRe, fEnergiesIm, fEnergyTot, fLambda;
double g = 0.01;

void init_data();

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::complex<double>> rapidities;
    double totEnergy, totLambda;
    init_data();
    for (int i = 0; i < rGSystem.energies.size(); i++)
    {
        for (int j = 0; j < rGSystem.occupation[i]; j++)
        {
            rapidities.push_back(rGSystem.energies[i]);
        }
    }
    fEnergiesRe.open("EnergiesRe.txt");
    fEnergiesIm.open("EnergiesIm.txt");
    fEnergyTot.open("EnergyTot.txt");
    fLambda.open("Lambda.txt");
    do
    {
        std::cout << "g = " << g << '\n';
        fEnergiesRe << g << '\t';
        fEnergiesIm << g << '\t';
        fLambda << g << '\t';
        cEquations test = cEquations(rGSystem.lambda, rGSystem.occupation, rGSystem.degeneracies, rGSystem.energies, g);
        test.newton_raphson();
        totLambda = 0;
        for (int i = 0; i < rGSystem.energies.size(); i++)
        {
            for (int j = 0; j < rGSystem.degeneracies[i]; j++)
            {
                fLambda << real(rGSystem.lambda[i][j]) << '\t';
            }
            //totLambda += real(rGSystem.lambda[i][0]);
        }
        //std::cout << "totLambda = " << totLambda << '\n';
        fLambda << '\n';
        cPolynomial test2 = cPolynomial(rGSystem.lambda, rGSystem.energies, rGSystem.degeneracies, rGSystem.occupation, g);
        test2.root_finder(rapidities);
        totEnergy = 0;
        for (int i = 0; i < rapidities.size(); i++)
        {
            fEnergiesRe << real(rapidities[i]) << '\t';
            totEnergy += real(rapidities[i]);
            fEnergiesIm << imag(rapidities[i]) << '\t';
        }
        fEnergyTot << g << '\t' << totEnergy << '\n';
        g += 0.005;
        fEnergiesIm << '\n';
        fEnergiesRe << '\n';
    } while (g <= 1.5);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> chronoDiff = end - start;
    std::cout << "Time elapsed: " << chronoDiff.count() << "s.\n";

    fEnergiesIm.close();
    fEnergiesRe.close();
    fEnergyTot.close();
    fLambda.close();
}

void init_data()
{
    rGSystem.energies = { 2,4,6,8,10};
    rGSystem.degeneracies = {1,2,3,4,5};
    rGSystem.occupation = {1,2,3,3,3};
    rGSystem.init_lambda();
}