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
	int totalLambda = 0;
	std::vector<std::vector<std::complex<double>>>& eLambda;
	std::vector<int>& eOccupation;
	std::vector<int>& eDegeneracies;
	std::vector<double>& eEnergies;
public:
	cEquations(std::vector<std::vector<std::complex<double>>>& lambdaData, std::vector<int>& occData, std::vector<int>& degData, std::vector<double>& enData, double g)
		:eLambda(lambdaData), eOccupation(occData), eDegeneracies(degData), eEnergies(enData), eG(g) {
		for (int i = 0; i < eEnergies.size(); i++)
		{
			for (int j = 0; j < eDegeneracies[i]; j++)
			{
				totalLambda++;
			}
		}
	}
	std::complex<double> function_i(std::vector<std::vector<std::complex<double>>> xVector, int iIndex);
	double binom(int n, int k);
	void newton_raphson();
};

class LUDecomp
{
private:
	std::vector<std::vector<std::complex<double>>>& mData;
	std::vector<std::vector<int>> parMatrix;
	bool isDecomposed = false;
public:
	LUDecomp(std::vector<std::vector<std::complex<double>>>& matrix)
		: mData(matrix), parMatrix(matrix.size(), std::vector<int>(matrix.size())) {
		decompose();
	};
	void decompose();
	void solve(std::vector<std::complex<double>> dataIn, std::vector<std::complex<double>>& dataOut);
	void row_swap(int idx1, int idx2, std::vector<std::vector<std::complex<double>>>& matrix);
	void row_swap(int idx1, int idx2, std::vector<std::vector<int>>& matrix);
};

class cPolynomial
{
private:
	double& pG;
	int pDegree = 0;
	std::vector<std::vector<std::complex<double>>>& pLambda;
	std::vector<double>& pEnergies;
	std::vector<int>& pOccupations;
	std::vector<int>& pDegeneracies;
	std::vector<std::complex<double>> pCoefficients;
	std::complex<double> firstDer;
	std::complex<double> secondDer;
public:
	cPolynomial(std::vector<std::vector<std::complex<double>>>& lambdaData, std::vector<double>& energyData, std::vector<int>& degData, std::vector<int>& occData, double& g)
		:pLambda(lambdaData), pEnergies(energyData), pDegeneracies(degData), pOccupations(occData), pG(g) {
		init_proper_lambda();
		find_coefficients();
	}
	double binom(int n, int k);
	std::vector<std::vector<std::complex<double>>> properLambda;
	std::complex<double> bell_partials(int nMax, int bellIndex, std::vector<std::complex<double>> currentLambda);
	void init_proper_lambda();
	void find_coefficients();
	std::complex<double> value(std::complex<double> zVal, int currentDegree);
	void deflate(std::complex<double> root);
	void root_finder(std::vector<std::complex<double>>& roots);
};
#endif // !RGDEGENLIB_H_
