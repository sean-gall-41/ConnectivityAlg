#ifndef CONNECTION_H_
#define CONNECTION_H_

#include <Random123/threefry.h>
#include <Random123/uniform.hpp>
#include <Random123/MicroURNG.hpp>
#include "dynamic2darray.h"
#include "conparams.h"
#include "celltype.h"

typedef r123::Threefry2x32 rngx32;

// TODO: utilize the algorithms suggested by rng123
unsigned int randInt(rngx32::ctr_type &c, rngx32::key_type &k,
	int a, int b);

float randReal(rngx32::ctr_type &c, rngx32::key_type &k,
	float a = 0.0, float b = 1.0);

bool nearly_equal(float a, float b, int factor /* factor of epsilon */);

class Connection
{
public:
	// obtain the correct connectivity parameters given the IDs of the cells
	Connection(CellType &srcCell, CellType &destCell, ConParams &conParams);
	Connection(const Connection &otherCon);
	~Connection();

	int getNumCons();
	float getAvgNumCons();
	float getRelRecipCons();
	void resetCons();

	void establishConnectionDecay(int maxNumAttempts,
		unsigned int seed32, float (*probFunc)(float, float, int, int));

	void establishConnectionCommon(bool needUnique,
		int normNumAttempts, int maxNumAttempts, unsigned int seed32);

	template<typename type>
	friend void shuffle(type *inArr, int arrSize, rngx32::ctr_type &c, rngx32::key_type &k);

	std::string toString();
	void toFile(const std::string outFile);	

	friend bool operator==(const Connection &thisCon, const Connection &otherCon);
	friend bool operator!=(const Connection &thisCon, const Connection &otherCon);
	Connection &operator=(const Connection &otherCon);

private:

// ======================== BEGIN GENERAL SRC DEST PARAMS ===============================

	// feel like these should be pointers?
	CellType srcCell;
	CellType destCell;

	// geometry
	int spanSrcOnDestX;
	int spanSrcOnDestY;
	int numPCon;
	// there should be two variables shouldn't there, like below?
	int numConToMake; // new var (04/21/2022)	
	//int numSrcCon;  // unused so far (04/21/2022)
	//int numDestCon; // unused so far (04/21/2022)

	// spatial probability
	float ampl;
	float stdDev;

	// reciprocal connections
	float pRecip;
	bool recipCons;
	float pRecipLowerBase;
	float reduceBaseRecip;

	// summary descriptive stats
	int totCon;
	float avgNumCon;
	float fracRecip;

	// input arrays
	int *destNumConArr;
	Dynamic2DArray<int> destConArr;

	//output arrays
	int *srcNumConArr;
	Dynamic2DArray<int> srcConArr;

// ========================== END GENERAL SRC DEST PARAMS ===============================

// ======================== BEGIN GENERAL PRIVATE METHODS =============================== 

	void deepCopy(const Connection &otherCon);
	void initConParams(ConParams &conParams);
	void initConParams(const Connection &otherCon);
	void allocateConArrs();
	void allocateConArrs(const Connection &otherCon);
	void initConArrs();
	void initConArrs(const Connection &otherCon);
	friend bool compareArrays(const Connection &thisCon, const Connection &otherCon);
	friend bool compareParams(const Connection &thisCon, const Connection &otherCon);
	void removeSelfConnection();
	void calcFracRecip();
};

#endif /* CONNECTION_H_ */

