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
	//std::cout << "TEMP2 = " << temp2 << '\t';
	for (int k = 0; k < equationIndex; k++)
	{
		sum1 += binom(equationIndex - 1, k) * xVector[energyIndex][k] * xVector[energyIndex][equationIndex - k - 1];
	}
	//std::cout << "SUM1 = " << sum1 << '\t';
	for (int k = 1; k < equationIndex; k++)
	{
		for (int i = 0; i < eEnergies.size(); i++)
		{
			if (i != energyIndex)
			{
				sum2 += (eG * eDegeneracies[i] * tgamma(equationIndex)) / (pow(eEnergies[i] - eEnergies[energyIndex], k) * tgamma(equationIndex - k+1)) * xVector[energyIndex][equationIndex - k];
			}
		}
	}
	//std::cout << "SUM2 = " << sum2 << '\n';
	sum3 = 0;
	for (int i = 0; i < eEnergies.size(); i++)
	{
		if (i != energyIndex)
		{
			//std::cout << "eG = " << eG << " eDegeneracies[i] = " << eDegeneracies[i] << " tgamma(equationIndex) = " << tgamma(equationIndex) << " xVector[enIdx][0] = " << xVector[energyIndex][0] << '\n';
			//std::cout << "xVector[i][0] = " << xVector[i][0] << " eEnergies[i] = " << eEnergies[i] << " eEnergies[energyIndex] = " << eEnergies[energyIndex] << '\n';
			sum3 += eG * eDegeneracies[i] * tgamma(equationIndex) * (xVector[energyIndex][0] - xVector[i][0]) / (pow(eEnergies[i] - eEnergies[energyIndex], equationIndex));
			//std::cout << "(i) = " << i << " current sum3 = " << sum3 << '\n';
		}
	}
	//std::cout << '\n';
	//std::cout << "SUM3 = " << sum3 << '\n';
	std::complex<double> result = temp2 - xVector[energyIndex][equationIndex - 1] + sum1 + sum2 + sum3;
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
					if (currentIndex == 0)
					{
						minusFVec[k] = -function_i(eLambda, k);
						//std::cout << "minusFVec[" << k << "] = " << minusFVec[k] << '\n';
					}
				}
				/*std::cout << "JACOBIAN ROW FINISHED" << '\n';

				std::cout << "func_i = " << -minusFVec[currentIndex] << '\n';
				std::cout << "MINUS F VEC ELEMENT FINISHED" << '\n';*/
				currentIndex++;
			}
		}
		/*std::cout << "JACOBIAN:" << '\n';
		for (int i = 0; i < totalLambda; i++)
		{
			for (int j = 0; j < totalLambda; j++)
			{
				std::cout << real(jacobian[i][j]) << '\t';
			}
			std::cout << '\n';
		}
		std::cout << "MINUS F VEC: " << '\n';
		for (int i = 0; i < totalLambda; i++)
		{
			std::cout << real(minusFVec[i]) << '\t';
		}
		std::cout << '\n';*/
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

void cPolynomial::init_proper_lambda()
{
	//std::cout << "PROPER LAMBDA\n";
	std::vector<std::complex<double>> tempVec;
	//std::cout << "PROPER LAMBDA\n";
	for (int i = 0; i < pEnergies.size(); i++)
	{
		tempVec.clear();
		if (pOccupations[i] > 0)
		{
			for (int j = 0; j < pOccupations[i]; j++)
			{
				tempVec.push_back(pLambda[i][j]/pG);
			}
			properLambda.push_back(tempVec);
			/*for (int k = 0; k < tempVec.size(); k++)
			{
				std::cout << real(tempVec[k]) << '\t';
			}
			std::cout << '\n';*/
		}
		else
		{
			tempVec.push_back(0);
			properLambda.push_back(tempVec);
		}
	}
	//std::cout << "PROPER LAMBDA SIZE = " << properLambda.size() << '\n';
}

std::complex<double> cPolynomial::bell_partials(int nMax, int bellIndex, std::vector<std::complex<double>> currentLambda)
{
	std::vector<std::vector<std::complex<double>>> bellPartials;
	std::vector<std::complex<double>> tempVec;
	std::complex<double> sum;
	for (int i = 0; i <= nMax; i++)
	{
		tempVec.push_back(0);
	}
	for (int i = 0; i <= nMax; i++)
	{
		bellPartials.push_back(tempVec);
	}
	bellPartials[0][0] = 1;
	for (int n = 1; n <= nMax; n++)
	{
		for (int k = 1; k <= nMax; k++)
		{
			sum = 0;
			for (int i = 1; i <= n-k+1; i++)
			{
				sum += binom(n - 1, i - 1) * currentLambda[i - 1] * bellPartials[n - i][k - 1];
			}
			bellPartials[n][k] = sum;
		}
	}
	/*for (int i = 0; i <= nMax; i++)
	{
		for (int j = 0; j <= nMax; j++)
		{
			std::cout << real(bellPartials[i][j]) << '\t';
		}
		std::cout << '\n';
	}*/
	std::complex<double> result = 0;
	for (int i = 0; i <= nMax; i++)
	{
		result += bellPartials[bellIndex][i];
	}
	//std::cout << "RESULT = " << result << '\n';
	return result;
}

void cPolynomial::find_coefficients()
{
	int pairCount = 0, kIndex = 1, energyIndex = 0;
	for (int i = 0; i < pEnergies.size(); i++)
	{
		pairCount += pOccupations[i];
	}
	std::vector<std::vector<std::complex<double>>> aMatrix(pairCount, std::vector<std::complex<double>>(pairCount));
	std::vector<std::complex<double>> bVector(pairCount);
	double insideProd;
	for (int i = 0; i < pairCount; i++)
	{
		if (kIndex > pOccupations[energyIndex])
		{
			kIndex = 1;
			energyIndex++;
			i--;
			continue;
		}
		//std::cout << "energyIndex = " << energyIndex << '\t' << "kIndex = " << kIndex << '\n';
		for (int m = 0; m < pairCount; m++)
		{
			insideProd = 1.0;
			for (int z = m; z > m - kIndex; z--)
			{
				insideProd *= z;
			}
			//std::cout << "bell_partials = " << bell_partials(pOccupations[energyIndex], kIndex, properLambda[energyIndex]) << '\n';
			//std::cout << "pow(pEnergies[energyIndex],m) = " << pow(pEnergies[energyIndex], m) << '\n';
			aMatrix[i][m] = bell_partials(pOccupations[energyIndex], kIndex, properLambda[energyIndex]) * pow(pEnergies[energyIndex], m) - insideProd*pow(pEnergies[energyIndex],m-kIndex);
		}
		insideProd = 1;
		for (int z = pairCount; z > pairCount - kIndex; z--)
		{
			insideProd *= z;
		}
		//std::cout << "insideProd = " << insideProd << '\n';
		//std::cout << "energyIndex = " << energyIndex << '\n';
		bVector[i] = insideProd * pow(pEnergies[energyIndex], pairCount - kIndex) - bell_partials(pOccupations[energyIndex], kIndex, properLambda[energyIndex]) * pow(pEnergies[energyIndex], pairCount);
		kIndex++;
	}
	/*std::cout << "AMATRIX:\n";
	for (int i = 0; i < pairCount; i++)
	{
		for (int j = 0; j < pairCount; j++)
		{
			std::cout << real(aMatrix[i][j]) << '\t';
		}
		std::cout << '\n';
	}
	std::cout << "BVECTOR:\n";
	for (int i = 0; i < pairCount; i++)
	{
		std::cout << real(bVector[i]) << '\t';
	}
	std::cout << '\n';*/
	LUDecomp coeffSolver = LUDecomp(aMatrix);
	coeffSolver.solve(bVector, bVector);
	pCoefficients.push_back(1);
	for (int i = 1; i <= pairCount; i++)
	{
		pCoefficients.push_back(bVector[pairCount-i]);
	}
	/*for (int i = 0; i < pairCount; i++)
	{
		pCoefficients.push_back(bVector[i]);
	}*/
	/*std::cout << "COEFFICIENTS:\n";
	for (int i = 0; i <= pairCount; i++)
	{
		std::cout << real(pCoefficients[i]) << '\t';
	}
	std::cout << '\n';*/
}

double cPolynomial::binom(int n, int k)
{
	return tgamma(n + 1) / (tgamma(k + 1) * tgamma(n - k + 1));
}

std::complex<double> cPolynomial::value(std::complex<double> zVal, int currentDegree)
{
	std::complex<double> cVal = pCoefficients[0];
	firstDer = static_cast<std::complex<double>>(currentDegree) * pCoefficients[0];
	secondDer = static_cast<std::complex<double>>(currentDegree * (currentDegree - 1)) * pCoefficients[0];
	for (int i = 1; i <= currentDegree; i++)
	{
		cVal = cVal * zVal + pCoefficients[i];
		if (i <= currentDegree - 1)
		{
			firstDer = firstDer * zVal + static_cast<std::complex<double>>(currentDegree - i) * pCoefficients[i];
		}
		if (i <= currentDegree - 2)
		{
			secondDer = secondDer * zVal + static_cast<std::complex<double>>((currentDegree - i) * (currentDegree - i - 1)) * pCoefficients[i];
		}
	}
	return cVal;
}

void cPolynomial::deflate(std::complex<double> root)
{
	for (int i = 1; i <= pDegree; i++)
	{
		pCoefficients[i] += pCoefficients[i - 1] * root;
	}
}

void cPolynomial::root_finder(std::vector<std::complex<double>>& roots)
{
	pDegree = roots.size();
	int total = roots.size();
	std::complex<double> initZ, val, aChange, root, gVal, hVal;
	for (int i = 0; i < total; i++)
	{
		initZ = roots[i];
		do
		{
			val = value(initZ, pDegree);
			if (abs(val) < 1e-10)
			{
				roots[i] = initZ;
				break;
			}
			gVal = firstDer / val;
			hVal = pow(gVal, 2) - secondDer / val;
			root = sqrt(static_cast<std::complex<double>>(pDegree - 1) * static_cast<std::complex<double>>(pDegree) * hVal - pow(gVal, 2));
			if (abs(gVal + root) > abs(gVal - root))
			{
				aChange = static_cast<std::complex<double>>(pDegree) / (gVal + root);
			}
			else
			{
				aChange = static_cast<std::complex<double>>(pDegree) / (gVal - root);
			}
			initZ -= aChange;
		} while (abs(aChange) > 1e-7);
		roots[i] = initZ;
		deflate(initZ);
		pDegree -= 1;
	}
}