#include <set>
#include <algorithm>
#include <cstring>
#include <climits>
#include <cmath>
#include <numeric>
#include "connection.h"

int Connection::getNumCons()
{
	return totCon;
}

float Connection::getAvgNumCons()
{
	return avgNumCon;
}

float Connection::getRelRecipCons()
{
	return fracRecip;
}

void Connection::resetCons()
{
	totCon    = 0;
	avgNumCon = 0.0;
	fracRecip = 0.0;	
	
	initConArrs();
}

Connection::Connection(CellType &srcCell, CellType &destCell, ConParams &conParams)
{
	// debatable whether we should have local copies of these cells
	this->srcCell  = srcCell;
	this->destCell = destCell;

	initConParams(conParams);
	allocateConArrs();
	initConArrs();
}

Connection::Connection(const Connection &otherCon)
{
	deepCopy(otherCon);
}

Connection::~Connection()
{
	// 2D Array deletion taken care of in destructor
	delete[] destNumConArr;
	delete[] srcNumConArr;
}

// based off of connectGOGODecayP in innetconnectivitystate
void Connection::establishConnectionDecay(int maxNumAttempts,
	unsigned int seed32, float (*probFunc2d)(float, float, int, int))
{
	int* spanArrX   = new int[spanSrcOnDestX + 1];
	int* spanArrY   = new int[spanSrcOnDestY + 1];
	int* spanIndArr = new int[numPCon];
	int* xCoordArr  = new int[numPCon];
	int* yCoordArr  = new int[numPCon];
	float* pConArr  = new float[numPCon];
	Dynamic2DArray<int> srcConBoolArr;
	srcConBoolArr.reshape(srcCell.getNumTotCells(), srcCell.getNumTotCells());

	// TODO: test whether we need to do this
	std::fill(spanArrX, spanArrX + (spanSrcOnDestX + 1), 0);
	std::fill(spanArrY, spanArrY + (spanSrcOnDestY + 1), 0);
	std::fill(spanIndArr, spanIndArr + numPCon, 0);	
	std::fill(xCoordArr, xCoordArr + numPCon, 0);
	std::fill(yCoordArr, yCoordArr + numPCon, 0);
	std::fill(pConArr, pConArr + numPCon, 0.0);
	srcConBoolArr.fill(0);

	for (int i = 0; i < spanSrcOnDestX + 1; i++)
	{
		spanArrX[i] = i - spanSrcOnDestX / 2;
	}

	for (int i = 0; i < spanSrcOnDestY + 1; i++)
	{
		spanArrY[i] = i - spanSrcOnDestY / 2;
	}

	for (int i = 0; i < numPCon; i++)
	{
		// TODO: check when the two span sizes are not the same
		xCoordArr[i] = spanArrX[i % (spanSrcOnDestX + 1)];
		yCoordArr[i] = spanArrY[i / (spanSrcOnDestX + 1)];
	}

	// Below is based on PDF that have two param like Gaussians.
	// Definitely not generalizable to PDF that have more or less params.	
	for (int i = 0; i < numPCon; i++)
	{
		pConArr[i] = probFunc2d(ampl, stdDev,
			xCoordArr[i], yCoordArr[i]);
	}
	
	// remove self-connections
	for (int i = 0; i < numPCon; i++)
	{
		if ((xCoordArr[i] == 0) && (yCoordArr[i] == 0)) pConArr[i] = 0.0;
	}

	// initialize the span indices array	
	for (int i = 0; i < numPCon; i++) spanIndArr[i] = i;

	// the positions we will be using to populate our connectivity arrays	
	int srcPosX, srcPosY;
	int destPosX, destPosY;
	int destIndex;	

	// initialize the random number generator counter and key
	rngx32::ctr_type c = {{0, 0}};
	rngx32::key_type k = {{seed32}};	
	r123::MicroURNG<rngx32> gen(c.incr(), k);

	for (int attempt = 1; attempt <= numConToMake; attempt++)
	{
		for (int i = 0; i < srcCell.getNumTotCells(); i++)
		{
			srcPosX = i % srcCell.getNumCellsX();
			srcPosY = i / srcCell.getNumCellsX();
	
			std::shuffle(spanIndArr, spanIndArr + numPCon, gen);

			for (int j = 0; j < numPCon; j++)
			{
				destPosX = srcPosX + xCoordArr[spanIndArr[j]];
				destPosY = srcPosY + yCoordArr[spanIndArr[j]];

				destPosX = (destPosX % srcCell.getNumCellsX() + srcCell.getNumCellsX())
					% srcCell.getNumCellsX();
				destPosY = (destPosY % srcCell.getNumCellsY() + srcCell.getNumCellsY())
					% srcCell.getNumCellsY();

				destIndex = destPosY * srcCell.getNumCellsX() + destPosX;

				// regular connections
				if (/*randReal(c, k) >= 1 - pConArr[spanIndArr[j]]
					&& */(srcConBoolArr(i, destIndex) == 0)
					&& srcNumConArr[i] < attempt
					&& destNumConArr[destIndex] < attempt)
				{
					srcConArr(i, srcNumConArr[i]) = destIndex;
					srcNumConArr[i]++;	

					destConArr(destIndex, destNumConArr[destIndex]) = i;
					destNumConArr[destIndex]++;	

					srcConBoolArr(i, destIndex) = 1;
					
					if (recipCons && randReal(c, k) >= 1 - pRecip 
						&& (srcConBoolArr(destIndex, i) == 0)
						&& srcNumConArr[destIndex] < attempt
						&& destNumConArr[i] < attempt)
					{
						srcConArr(destIndex, srcNumConArr[destIndex]) = i;
						srcNumConArr[destIndex]++;	

						destConArr(i, destNumConArr[i]) = destIndex;
						destNumConArr[i]++;	
	
						srcConBoolArr(destIndex, i) = 1;		
					}	
				}
				// using a decay introduces a bias, which is what the other
				// code in innetconnectivity did
			}
		}
	}

	// check the bounds here
	for (int i = 0; i < destCell.getNumTotCells(); i++)
	{
			totCon += destNumConArr[i];
	}

	avgNumCon = totCon / float(destCell.getNumTotCells());
	if (recipCons) calcFracRecip();	

	delete[] spanArrX; 
	delete[] spanArrY;
	delete[] spanIndArr;
	delete[] xCoordArr;
	delete[] yCoordArr;
	delete[] pConArr;
}

void Connection::establishConnectionCommon(bool needUnique,
	int normNumAttempts, int maxNumAttempts, unsigned int seed32)
{
	bool *srcConnected = new bool[srcCell.getNumTotCells()];

	rngx32::ctr_type c = {{0,0}};
	rngx32::key_type k = {{seed32}};	

	int srcInd;
	int srcPosX;
	int srcPosY;
	int tempDestPosX;
	int tempDestPosY;
	int derivedDestInd;
	int tempDestNumConLim;

	int numSrcConnected;


	for (int i = 0; i < numConToMake; i++)
	{
		numSrcConnected = 0;
		memset(srcConnected, false, srcCell.getNumTotCells() * sizeof(bool));
		// our randInt function returns value in half-open set [a, b)
		while (numSrcConnected < srcCell.getNumTotCells())
		{
			srcInd = randInt(c, k, 0, srcCell.getNumTotCells()); 
			if (!srcConnected[srcInd])
			{
				srcPosX = srcInd % srcCell.getNumCellsX();
				srcPosY = srcInd / srcCell.getNumCellsX();	
				// fully redundant at this point (04/22/22)
				tempDestNumConLim = numConToMake;			
				int attempts;
				for (attempts = 0; attempts < maxNumAttempts; attempts++)
				{
					if (attempts == normNumAttempts) tempDestNumConLim = numConToMake; 

					tempDestPosX = srcPosX;
					tempDestPosY = srcPosY;

					tempDestPosX += (int)std::round((randReal(c, k) - 0.5) * (spanSrcOnDestX + 1));		
					tempDestPosY += (int)std::round((randReal(c, k) - 0.5) * (spanSrcOnDestY + 1));

					tempDestPosX = (tempDestPosX % srcCell.getNumCellsX() + srcCell.getNumCellsX())
						% srcCell.getNumCellsX();	
					tempDestPosY = (tempDestPosY % srcCell.getNumCellsY() + srcCell.getNumCellsY())
						% srcCell.getNumCellsY();	
					
					derivedDestInd = tempDestPosY * srcCell.getNumCellsX() + tempDestPosX;	

					bool unique = true;

					if (needUnique)
					{
						for (int j = 0; j < i; j++)
						{
							if (derivedDestInd == srcConArr(srcInd, j))
							{
								unique = false;
								break;
							}
						}
					}
					
					if (unique)
					{
						if (destNumConArr[derivedDestInd] < tempDestNumConLim
								/*&& srcNumConArr[srcInd] < tempDestNumConLim */)
						{
							destConArr(derivedDestInd, destNumConArr[derivedDestInd]) = srcInd;
							destNumConArr[derivedDestInd]++;

							srcConArr(srcInd, srcNumConArr[srcInd]) = derivedDestInd;
							srcNumConArr[srcInd]++;
							break;
						}
					}
				}

				srcConnected[srcInd] = true;
				numSrcConnected++;	
			}
		}	
	}

	for (int i = 0; i < destCell.getNumTotCells(); i++)
	{
			totCon += destNumConArr[i];
	}
	avgNumCon = totCon / (float)destCell.getNumTotCells();
	
	delete[] srcConnected;
}

void Connection::deepCopy(const Connection &otherCon)
{
	// copy the cells
	this->srcCell  = otherCon.srcCell;
	this->destCell = otherCon.destCell;

	// copy the dest and src arrs (note: overloaded assignment for dynamic2darrs)
	this->destConArr = otherCon.destConArr;
	this->srcConArr  = otherCon.srcConArr;

	// three things to do for non-vector arrays:
	// 1) initialize connectivity params
	// 2) allocate the connectivity arrays
	// 3) initialize the connectivity arrays
	initConParams(otherCon);
	allocateConArrs(otherCon);
	initConArrs(otherCon);
}

void Connection::initConParams(ConParams &conParams)
{
	spanSrcOnDestX  = std::stoi(conParams["spanSrcOnDestX"]);
	spanSrcOnDestY  = std::stoi(conParams["spanSrcOnDestY"]);
	numPCon		    = std::stoi(conParams["numPCon"]);
	numConToMake 	= std::stoi(conParams["numConToMake"]);
	
	ampl		    = std::stof(conParams["ampl"]); 
	stdDev		    = std::stof(conParams["stdDev"]);

	pRecip	        = std::stof(conParams["pRecip"]);
	// they could've made a stob but noooo we have to do it via other means
	// also not a great solution as strcasecmp is non-cross platform
	recipCons	    = (strcasecmp("true", conParams["recipCons"].c_str()) == 0);
	pRecipLowerBase = std::stof(conParams["pRecipLowerBase"]);
	reduceBaseRecip = (strcasecmp("true", conParams["reduceBaseRecip"].c_str()) == 0);

	//normNumGOGOIn  = conParams.conParamMap["normNumGOGOIn"];
	//normNumGOGOOut = conParams.conParamMap["normNumGOGOOut"];
	//maxNumGOGOIn   = conParams.conParamMap["maxNumGOGOIn"];
	//maxNumGOGOOut  = conParams.conParamMap["maxNumGOGOOut"];

	totCon  		= 0;
	avgNumCon 		= 0.0;
	fracRecip     	= 0.0;	
}

void Connection::initConParams(const Connection &otherCon)
{
	this->spanSrcOnDestX  = otherCon.spanSrcOnDestX;
	this->spanSrcOnDestY  = otherCon.spanSrcOnDestY;
	this->numPCon		  = otherCon.numPCon;
	this->numConToMake 	  = otherCon.numConToMake;

	this->ampl		      = otherCon.ampl; 
	this->stdDev		  = otherCon.stdDev;

	this->pRecip	      = otherCon.pRecip;

	// they could've made a stob but noooo we have to do it via other means
	// also not a great solution as strcasecmp is non-cross platform
	this->recipCons	      = otherCon.recipCons;
	this->pRecipLowerBase = otherCon.pRecipLowerBase;
	this->reduceBaseRecip = otherCon.reduceBaseRecip;

	//normNumGOGOIn  = conParams.conParamMap["normNumGOGOIn"];
	//normNumGOGOOut = conParams.conParamMap["normNumGOGOOut"];
	//maxNumGOGOIn   = conParams.conParamMap["maxNumGOGOIn"];
	//maxNumGOGOOut  = conParams.conParamMap["maxNumGOGOOut"];

	this->totCon  		  = otherCon.totCon;
	this->avgNumCon 	  = otherCon.avgNumCon;
	this->fracRecip       = otherCon.fracRecip;	
}

void Connection::allocateConArrs()
{
	// allocate all arrays (reshape the 2D Ones)
	destNumConArr = new int[destCell.getNumTotCells()];
	destConArr.reshape(destCell.getNumTotCells(), numConToMake);
	
	srcNumConArr  = new int[srcCell.getNumTotCells()];
	srcConArr.reshape(srcCell.getNumTotCells(), numConToMake);
}

void Connection::allocateConArrs(const Connection &otherCon)
{
	destNumConArr = new int[otherCon.destCell.getNumTotCells()];
	srcNumConArr  = new int[otherCon.srcCell.getNumTotCells()];
}

void Connection::initConArrs()
{
	// initialize all the global connection arrays 
	std::fill(destNumConArr, destNumConArr + destCell.getNumTotCells(), 0);
	destConArr.fill(INT_MAX);

	std::fill(srcNumConArr, srcNumConArr + srcCell.getNumTotCells(), 0);
	srcConArr.fill(INT_MAX);
}

void Connection::initConArrs(const Connection &otherCon)
{
	for (int i = 0; i < otherCon.destCell.getNumTotCells(); i++)
	{
		this->destNumConArr[i] = otherCon.destNumConArr[i];
	}

	for (int i = 0; i < otherCon.srcCell.getNumTotCells(); i++)
	{
		this->srcNumConArr[i] = otherCon.srcNumConArr[i];
	}
}

bool compareArrays(const Connection &thisCon, const Connection &otherCon)
{
	bool equal = true;	
	// check the dest num arrays for equality
	if (thisCon.destCell.getNumTotCells() == otherCon.destCell.getNumTotCells())
	{
		for(int i = 0; i < thisCon.destCell.getNumTotCells(); i++)
		{
			if (thisCon.destNumConArr[i] != otherCon.destNumConArr[i])
			{
				equal = false;
				break;	
			}
		}	
	}
	else {equal = false;}

	// check the src num arrays for equality
	if (equal && thisCon.srcCell.getNumTotCells() == otherCon.srcCell.getNumTotCells())
	{
		for(int i = 0; i < thisCon.srcCell.getNumTotCells(); i++)
		{
			if (thisCon.srcNumConArr[i] != otherCon.srcNumConArr[i])
			{
				equal = false;
				break;	
			}
		}	
	}
	else {equal = false;}

	// check the 2D connectivity arrays for equality
	equal = equal && (thisCon.destConArr == otherCon.destConArr);
	equal = equal && (thisCon.srcConArr == otherCon.srcConArr);

	return equal;	
}

bool compareParams(const Connection &thisCon, const Connection &otherCon)
{
	bool equal = true;

	equal = equal && (thisCon.destCell == otherCon.destCell); 
	equal = equal && (thisCon.srcCell == otherCon.srcCell); 

	equal = equal && (thisCon.spanSrcOnDestX == otherCon.spanSrcOnDestX);
	equal = equal && (thisCon.spanSrcOnDestY == otherCon.spanSrcOnDestY);
	equal = equal && (thisCon.numPCon == otherCon.numPCon);
	equal = equal && (thisCon.numConToMake == otherCon.numConToMake);

	// using factor of 1000 times machine epsilon to determine whether we are equal 
	// something to play with if you ever use these equality operations
	equal = equal && nearly_equal(thisCon.ampl, otherCon.ampl, 1000);
	equal = equal && nearly_equal(thisCon.stdDev, otherCon.stdDev, 1000);
	equal = equal && nearly_equal(thisCon.pRecip, otherCon.pRecip, 1000);
	equal = equal && nearly_equal(thisCon.pRecipLowerBase, otherCon.pRecipLowerBase, 1000);
	equal = equal && nearly_equal(thisCon.reduceBaseRecip, otherCon.reduceBaseRecip, 1000);

	// also iffy on these factors of epsilon...
	equal = equal && nearly_equal(thisCon.avgNumCon, otherCon.avgNumCon, 100000);
	equal = equal && nearly_equal(thisCon.fracRecip, otherCon.fracRecip, 10000);

	equal = equal && (thisCon.recipCons && otherCon.recipCons);
	equal = equal && (thisCon.totCon && otherCon.totCon);

	return equal;
}


void Connection::calcFracRecip()
{
	//WARNING: ONLY VALID IF numSrcCells() == numDestCells()!!! 04/22/22
	int numRecip = 0;
	for (int i = 0; i < srcCell.getNumTotCells(); i++)
	{
		for (int j = 0; j < destNumConArr[i]; j++)
		{
			for (int k = 0; k < srcNumConArr[i]; k++)
			{
				if ((destConArr(i, j) == srcConArr(i, k))
						&& (destConArr(i, j) != INT_MAX)
						&& (srcConArr(i, k) != INT_MAX))
				{
					numRecip++;					
				}
			}
		}
	}
	fracRecip = numRecip / (float)totCon;	
}

std::string Connection::toString()
{
	std::string outString = "";
	outString += "[ totCon: " + std::to_string(totCon) + " ]\n";
	outString += "[ avgNumCon: " + std::to_string(avgNumCon) + " ]\n";
	if (recipCons)
	{
		outString += "[ fracRecip: " + std::to_string(fracRecip) + " ]\n";
	}
	return outString;
}

void Connection::toFile(const std::string outFileName)
{
	std::ofstream outFile(outFileName, std::ios::out | std::ios::binary);
	if (!outFile.is_open())
	{
		std::cerr << "[ERROR]: cannot find file " << outFileName 
			<< ". Aborting all write operations." << std::endl;
		exit(1);
	}
	outFile.write(reinterpret_cast<const char*>(destNumConArr),
		destCell.getNumTotCells() * sizeof(destNumConArr[0]));
	outFile.write(reinterpret_cast<const char*>(destConArr.data()),
		destConArr.size() * sizeof(destConArr(0, 0)));
	outFile.write(reinterpret_cast<const char*>(srcNumConArr),
		srcCell.getNumTotCells() * sizeof(srcNumConArr[0]));
	outFile.write(reinterpret_cast<const char*>(srcConArr.data()),
		srcConArr.size() * sizeof(srcConArr(0, 0)));

	outFile.close();
}

bool operator==(const Connection &thisCon, const Connection &otherCon)
{
	return compareArrays(thisCon, otherCon)
		&& compareParams(thisCon, otherCon);	
}

bool operator!=(const Connection &thisCon, const Connection &otherCon)
{
	return !(thisCon == otherCon);
}

Connection &Connection::operator=(const Connection &otherCon)
{
	if (this != &otherCon)
	{
		delete[] destNumConArr;
		delete[] srcNumConArr;	
		deepCopy(otherCon);	
	}
	return *this;
}

unsigned int randInt(rngx32::ctr_type &c, rngx32::key_type &k,
	int a, int b)
{
	c[0]++; // increment counter
	rngx32 g;
	rngx32::ctr_type r = g(c, k);
	return r.v[0] % (b - a) + a;	
}

float randReal(rngx32::ctr_type &c, rngx32::key_type &k,
	float a, float b)
{
	c[0]++; // increment counter	
	rngx32::ctr_type r = threefry2x32(c, k);
	float returnVal = r123::u01<float>(r.v[0]);
	if (a == 0.0 && b == 1.0) return returnVal;	
	return (b - a) * returnVal + a;	
}

bool nearly_equal(float a, float b, int factor /* factor of nextafter */)
{
	float min_a = a - (a - std::nextafter(a, std::numeric_limits<float>::lowest())) * factor;
	float max_a = a + (std::nextafter(a, std::numeric_limits<float>::max()) - a) * factor;

	return (min_a <= b) && (max_a >= b);
}

