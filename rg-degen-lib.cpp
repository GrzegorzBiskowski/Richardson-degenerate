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

void cEquations::newton_raphson()
{
	const std::complex<double> PREC = 1e-4 * sqrt(static_cast<std::complex<double>>(-1));
	std::complex<double> deltaNorm;
	//std::cout << "PREC = " << PREC << "\tIMAG(PREC) = " << imag(PREC) << '\n';
	int currentIndex, iteration = 0;
	std::vector<std::vector<std::complex<double>>> jacobian(totalLambda, std::vector<std::complex<double>>(totalLambda));
	std::vector<std::complex<double>> minusFVec(totalLambda);
	do
	{
		//std::cout << "iteration = " << iteration << '\n';
		currentIndex = 0;
		for (int i = 0; i < eEnergies.size(); i++)
		{
			for (int j = 0; j < eDegeneracies[i]; j++)
			{
				std::vector<std::vector<std::complex<double>>> dLambda = eLambda;
				dLambda[i][j] += PREC;
				for (int k = 0; k < totalLambda; k++)
				{
					jacobian[k][currentIndex] = imag(function_i(dLambda, k)) / imag(PREC);
				}
				minusFVec[currentIndex] = -function_i(eLambda, currentIndex);
				currentIndex++;
			}
		}
		LUDecomp decomposition = LUDecomp(jacobian);
		decomposition.solve(minusFVec, minusFVec);
		currentIndex = 0;
		deltaNorm = 0;
		for (int i = 0; i < eEnergies.size(); i++)
		{
			for (int j = 0; j < eDegeneracies[i]; j++)
			{
				eLambda[i][j] += minusFVec[currentIndex];
				deltaNorm += norm(minusFVec[currentIndex]);
				currentIndex++;
			}
		}
		iteration++;
	}while (abs(sqrt(deltaNorm)) > 1e-4);
	/*for (int i = 0; i < totalLambda; i++)
	{
		std::cout << real(function_i(eLambda, i)) << '\t';
	}*/
}

void LUDecomp::decompose()
{
	std::vector<std::complex<double>> tempCol(mData.size());
	std::vector<std::complex<double>> yVec(mData.size());
	std::vector<std::complex<double>> scaleVec(mData.size());
	std::complex<double> big = 0, temp;
	int bigIdx = 0;
	for (int i = 0; i < mData.size(); i++)
	{
		for (int j = 0; j < mData.size(); j++)
		{
			parMatrix[i][j] = 0;
		}
		parMatrix[i][i] = 1;
	}
	for (int i = 0; i < mData.size(); i++)
	{
		big = 0;
		for (int j = 0; j < mData.size(); j++)
		{
			temp = abs(mData[i][j]);
			if (abs(temp) > abs(big))
			{
				big = temp;
			}
		}
		scaleVec[i] = pow(big, -1);
		if (abs(big) == 0)
		{
			std::cout << "SINGULAR MATRIX!\n";
		}
	}
	for (int k = 0; k < mData.size(); k++)
	{
		big = 0;
		for (int i = k; i < mData.size(); i++)
		{
			temp = scaleVec[i] * abs(mData[i][k]);
			if (abs(temp) > abs(big))
			{
				big = temp;
				bigIdx = i;
			}
		}
		if (k != bigIdx)
		{
			for (int j = 0; j < mData.size(); j++)
			{
				row_swap(bigIdx, k, mData);
				row_swap(bigIdx, k, parMatrix);
			}
			temp = scaleVec[k];
			scaleVec[k] = scaleVec[bigIdx];
			scaleVec[bigIdx] = temp;
		}
		if (mData[k][k] == 0.0)
		{
			mData[k][k] = 1e-40;
		}
		for (int i = k + 1; i < mData.size(); i++)
		{
			temp = mData[i][k] /= mData[k][k];
			for (int j = k + 1; j < mData.size(); j++)
			{
				mData[i][j] -= temp * mData[k][j];
			}
		}
	}
	isDecomposed = true;
}

void LUDecomp::solve(std::vector<std::complex<double>> dataIn, std::vector<std::complex<double>>& dataOut)
{
	std::vector<std::complex<double>> yVec(dataIn.size());
	std::vector<std::complex<double>> tempVec(dataIn.size());
	std::vector<std::complex<double>> bVec(dataIn.size());
	std::complex<double> temp;
	for (int i = 0; i < dataIn.size(); i++)
	{
		bVec[i] = 0;
		for (int j = 0; j < dataIn.size(); j++)
		{
			bVec[i] += static_cast<std::complex<double>>(parMatrix[i][j]) * dataIn[j];
		}
	}
	yVec[0] = bVec[0];
	for (int i = 1; i < dataIn.size(); i++)
	{
		temp = 0;
		for (int j = 0; j < i; j++)
		{
			temp += mData[i][j] * yVec[j];
		}
		yVec[i] = bVec[i] - temp;
	}
	dataOut[dataIn.size() - 1] = yVec[dataIn.size() - 1] / mData[dataIn.size() - 1][dataIn.size() - 1];
	for (int i = dataIn.size() - 2; i >= 0; i--)
	{
		temp = 0;
		for (int j = i + 1; j < dataIn.size(); j++)
		{
			temp += mData[i][j] * dataOut[j];
		}
		dataOut[i] = (yVec[i] - temp) / mData[i][i];
	}
}

void LUDecomp::row_swap(int idx1, int idx2, std::vector<std::vector<std::complex<double>>>& matrix)
{
	std::complex<double> temp;
	for (int i = 0; i < matrix.size(); i++)
	{
		temp = matrix[idx1][i];
		matrix[idx1][i] = matrix[idx2][i];
		matrix[idx2][i] = temp;
	}
}

void LUDecomp::row_swap(int idx1, int idx2, std::vector<std::vector<int>>& matrix)
{
	int temp;
	for (int i = 0; i < matrix.size(); i++)
	{
		temp = matrix[idx1][i];
		matrix[idx1][i] = matrix[idx2][i];
		matrix[idx2][i] = temp;
	}
}