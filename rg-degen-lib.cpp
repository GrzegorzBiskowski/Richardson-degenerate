#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

#include "rg-degen-lib.h"

void systemData::init_lambda()
{
	std::vector<std::complex<double>> temp;
	for (int i = 0; i < energies.size(); i++)
	{
		temp.clear();
		for (int j = 0; j < degeneracies[i]; j++)
		{
			if (j == 0)
			{
				if (occupation[i] > 0)
				{
					temp.push_back(1);
				}
				else
				{
					temp.push_back(0);
				}
			}
			else
			{
				temp.push_back(0);
			}
		}
		lambda.push_back(temp);
	}
}

std::complex<double> cEquations::function_i(std::vector<std::vector<std::complex<double>>> xVector, int iIndex)
{
	int energyIndex = 0, equationIndex = 0, temp = 0;
	for (int i = 0; i < eDegeneracies.size(); i++)
	{
		temp += eDegeneracies[i];
		if (iIndex < temp)
		{
			energyIndex = i;
			equationIndex = iIndex + eDegeneracies[i] - temp+1;
			break;
		}
	}
	//std::cout << "energyIndex = " << energyIndex << '\t' << "equationIndex = " << equationIndex << '\n';
	std::complex<double> sum1 = 0, sum2 = 0, sum3 = 0, temp2;
	if (equationIndex < eDegeneracies[energyIndex])
	{
		temp2 = eG * (1 - eDegeneracies[energyIndex] / equationIndex) * xVector[energyIndex][equationIndex];
	}
	else
	{
		temp2 = 0.0;
	}
	//std::cout << "TEMP2 = " << temp2 << '\n';
	for (int k = 0; k < equationIndex; k++)
	{
		sum1 += binom(equationIndex - 1, k) * xVector[energyIndex][k] * xVector[energyIndex][equationIndex - k - 1];
	}
	//std::cout << "SUM1 = " << sum1 << '\n';
	for (int k = 1; k < equationIndex; k++)
	{
		for (int i = 0; i < eEnergies.size(); i++)
		{
			if (i != energyIndex)
			{
				sum2 += (eG * eDegeneracies[i] * tgamma(equationIndex)) / (pow(eEnergies[i] - eEnergies[energyIndex], k) * tgamma(equationIndex - k - 1)) * xVector[energyIndex][equationIndex - k];
			}
		}
	}
	//std::cout << "SUM2 = " << sum2 << '\n';
	for (int i = 0; i < eEnergies.size(); i++)
	{
		if (i != energyIndex)
		{
			sum3 += eG * eDegeneracies[i] * tgamma(equationIndex) * (xVector[energyIndex][0] - xVector[i][0]) / (pow(eEnergies[i] - eEnergies[energyIndex], equationIndex));
		}
	}
	//std::cout << "SUM3 = " << sum3 << '\n';
	std::complex<double> result = temp2 - xVector[energyIndex][equationIndex - 1] + sum1 - sum2 + sum3;
	return result;
}

double cEquations::binom(int n, int k)
{
	return tgamma(n + 1) / (tgamma(k + 1) * tgamma(n - k + 1));
}