#include "cellparams.h"
#include "celltype.h"
#include "conparams.h"
#include "connection.h"

const std::string ROOT 				   = "../../";
const std::string DATA_DIR 			   = ROOT + "data/";
const std::string INPUT_DIR 		   = DATA_DIR + "inputs/";
const std::string CELL_PARAM_FILE	   = INPUT_DIR + "cellParams.json";
const std::string CON_PARAM_FILE 	   = INPUT_DIR + "conParams.json";

float gaussian2d(float ampl, float stdDev, int xCoord, int yCoord);

int main()
{
	CellParams golgiParams(CELL_PARAM_FILE, "GO");
	CellType golgi(golgiParams);

	ConParams gogoConParams(CON_PARAM_FILE, "GOGO");
	
	Connection gogoConnection(golgi, golgi, gogoConParams);
	gogoConnection.establishConnectionDecay(12, 0, gaussian2d);

	//TODO: Here is where we will make our assertions.
	std::cout << "Method 'establishConnectionDecay' results:" << std::endl; 
	std::cout << "Number of Golgi-Golgi connections: " 
		<< gogoConnection.getNumCons() << std::endl;
	std::cout << "Average Number of Golgi-Golgi connections: "
		<< gogoConnection.getAvgNumCons() << std::endl;
	std::cout << "Percent reciprocal Golgi-Golgi connections: "
		<< (gogoConnection.getRelRecipCons() * 100) << std::endl;

	return 0;	
}

float gaussian2d(float ampl, float stdDev, int xCoord, int yCoord)
{ 
	float xCoordFloat = (xCoord * xCoord) / (2 * stdDev * stdDev);
	float yCoordFloat = (yCoord * yCoord) / (2 * stdDev * stdDev);
	return ampl * exp(-(xCoordFloat + yCoordFloat));
}

