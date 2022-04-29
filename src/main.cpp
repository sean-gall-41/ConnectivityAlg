#include <math.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <Random123/threefry.h>
#include <Random123/uniform.hpp>
#include "connection.h"
#include "params.h"
#include "conparams.h"
#include "contypeparams.h"
#include "celltype.h"

// Keep in mind root is relative to the executable (yuck)
const std::string ROOT 				  = "../../";
const std::string DATA_DIR 			  = ROOT + "data/";
const std::string INPUT_DIR 		  = DATA_DIR + "inputs/";
const std::string CELL_PARAM_FILE	  = INPUT_DIR + "cell_params.json";
const std::string CON_PARAM_FILE 	  = INPUT_DIR + "con_params.json";
const std::string CON_TYPE_PARAM_FILE = INPUT_DIR + "con_type_params.json";

float gaussian2d(float ampl, float stdDev, int xCoord, int yCoord);

int main()
{
	/* TODO:
	 * create the connectivity state
	 * print the state itself
	*/

	// define Cells	
	CellParams golgiParams(CELL_PARAM_FILE, "GO");
	CellType golgi(golgiParams);
	std::cout << golgi << std::endl;

	// define connection parameters
	ConParams gogoConParams(CON_PARAM_FILE, "GOGO");
	std::cout << gogoConParams << std::endl;

	// define connection type parameters 
	ConTypeParams gogoConTypeParams(CON_TYPE_PARAM_FILE);
	std::cout << gogoConTypeParams << std::endl;

	// create connection
	Connection gogoConnection(golgi, golgi, gogoConParams);
	gogoConnection.establishConnectionDecay(std::stoi(gogoConTypeParams["maxConAttempts"]),
		std::stoi(gogoConTypeParams["randSeed"]), gaussian2d);

	std::cout << "Method 'establishConnectionDecay' results:" << std::endl; 
	std::cout << "Number of Golgi-Golgi connections: " 
		<< gogoConnection.getNumCons() << std::endl;
	std::cout << "Average Number of Golgi-Golgi connections: "
		<< gogoConnection.getAvgNumCons() << std::endl;
	std::cout << "Percent reciprocal Golgi-Golgi connections: "
		<< (gogoConnection.getRelRecipCons() * 100) << std::endl;

	//ConParams pfgoConParams(CON_PARAM_FILE, "PFGO");
	//ConParams aagoConParams(CON_PARAM_FILE, "AAGO");

	return 0;
}

float gaussian2d(float ampl, float stdDev, int xCoord, int yCoord)
{ 
	float xCoordFloat = (xCoord * xCoord) / (2 * stdDev * stdDev);
	float yCoordFloat = (yCoord * yCoord) / (2 * stdDev * stdDev);
	return ampl * exp(-(xCoordFloat + yCoordFloat));
}

