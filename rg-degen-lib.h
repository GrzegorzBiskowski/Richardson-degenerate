#ifndef RGDEGENLIB_H_
#define RGDEGENLIB_H_

#include <complex>
#include <vector>

class systemData
{
private:
public:
	std::vector<double> energies;
	std::vector<int> degeneracies;
	std::vector<int> occupation;
	std::vector<std::vector<std::complex<double>>> lambda;
	void init_lambda();
};


class cEquations
{
private:
	double& eG;
	std::vector<std::vector<std::complex<double>>>& eLambda;
	std::vector<int>& eOccupation;
	std::vector<int>& eDegeneracies;
	std::vector<double>& eEnergies;
public:
	cEquations(std::vector<std::vector<std::complex<double>>>& lambdaData, std::vector<int>& occData, std::vector<int>& degData, std::vector<double>& enData, double g)
		:eLambda(lambdaData), eOccupation(occData), eDegeneracies(degData), eEnergies(enData), eG(g) {
	}
	std::complex<double> function_i(std::vector<std::vector<std::complex<double>>> xVector, int iIndex);
	double binom(int n, int k);
};
#endif // !RGDEGENLIB_H_
