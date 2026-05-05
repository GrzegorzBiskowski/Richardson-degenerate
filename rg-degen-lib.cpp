#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

#include "rg-degen-lib.h"

void systemData::init_lambda()
{
	int lambdaCount = 0;
	int maxDegen = 0;
	for (int i = 0; i < degeneracies.size(); i++)
	{
		if (degeneracies[i] > maxDegen)
		{
			maxDegen = degeneracies[i];
		}
		lambdaCount += degeneracies[i];
	}
	for (int i = 0; i < lambdaCount; i++)
	{
		if (i < energies.size())
		{
			if (occupation[i] == 0)
			{
				lambda.push_back(0);
			}
			else
			{
				lambda.push_back(1);
			}
		}
		else
		{
			lambda.push_back(0);
		}
	}
}

std::complex<double> cEquations::function_i(std::vector<std::complex<double>> xVector, int iIndex)
{
	int energyIndex = 0, equationIndex = 0, temp = 0;
	for (int i = 0; i < eDegeneracies.size(); i++)
	{
		std::cout << "i = " << i << '\n';
		temp += eDegeneracies[i];
		std::cout << "temp = " << temp << '\n';
		if (iIndex < temp)
		{
			energyIndex = i;
			equationIndex = iIndex + eDegeneracies[i] - temp;
			break;
		}
	}
	return 0;
}