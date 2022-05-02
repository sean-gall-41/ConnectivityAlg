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
float uniform(float dummy_1, float dummy_2, int dummy_3, int dummy_4);

int main()
{
	/* TODO:
	 * create the connectivity state
	 * print the state itself
	*/

	// define Cells	
	CellParams golgiParams(CELL_PARAM_FILE, "GO");
	CellType golgi(golgiParams);

	CellParams granuleParams(CELL_PARAM_FILE, "GR");
	CellType granule(granuleParams);

	// define connection parameters
	ConParams gogoConParams(CON_PARAM_FILE, "GOGO");
	ConParams pfgoConParams(CON_PARAM_FILE, "PFGO");
	ConParams aagoConParams(CON_PARAM_FILE, "AAGO");

	// define connection type parameters 
	ConTypeParams conTypeParams(CON_TYPE_PARAM_FILE);

	// establish golgi-golgi connections
	Connection gogoConnection(golgi, golgi, gogoConParams);
	gogoConnection.establishConnectionJoe(
		std::stoi(gogoConParams["maxNumAttempts"]),
		std::stoi(conTypeParams["randSeed"]),
		uniform);

	std::cout << "=== Golgi - Golgi Connections ===" << std::endl;

	std::cout << "Method 'establishConnectionJoe' results:" << std::endl; 
	std::cout << "Number of Golgi-Golgi connections: " 
		<< gogoConnection.getNumCons() << std::endl;
	std::cout << "Average Number of Golgi-Golgi connections: "
		<< gogoConnection.getAvgNumCons() << std::endl;
	std::cout << "Percent reciprocal Golgi-Golgi connections: "
		<< (gogoConnection.getRelRecipCons() * 100) << std::endl;

	// establish parallel-fiber-golgi connections
	Connection pfgoConnection(granule, golgi, pfgoConParams);	
	pfgoConnection.establishConnectionJoe(
		std::stoi(pfgoConParams["maxNumAttempts"]),
		std::stoi(conTypeParams["randSeed"]),
		uniform);

	std::cout << std::endl;
	
	std::cout << "=== Parallel Fiber - Golgi Connections ===" << std::endl;

	std::cout << "Method 'establishConnectionJoe' results:" << std::endl; 
	std::cout << "Number of Parallel Fiber - Golgi connections: " 
		<< pfgoConnection.getNumCons() << std::endl;
	std::cout << "Average Number of parallel fiber - Golgi connections: "
		<< pfgoConnection.getAvgNumCons() << std::endl;

	// TODO: add extra input variable for multiple values relating
	// to num connections to make
	Connection aagoConnection(granule, golgi, aagoConParams);
	aagoConnection.establishConnectionJoe(
		std::stoi(aagoConParams["maxNumAttempts"]),
		std::stoi(conTypeParams["randSeed"]),
		uniform);

	std::cout << std::endl;
	
	std::cout << "=== ascending axon - Golgi Connections ===" << std::endl;

	std::cout << "Method 'establishConnectionJoe' results:" << std::endl; 
	std::cout << "Number of ascending axon - Golgi connections: " 
		<< pfgoConnection.getNumCons() << std::endl;
	std::cout << "Average Number of ascending axon - Golgi connections: "
		<< pfgoConnection.getAvgNumCons() << std::endl;

	return 0;
}


float gaussian2d(float ampl, float stdDev, int xCoord, int yCoord)
{ 
	float xCoordFloat = (xCoord * xCoord) / (2 * stdDev * stdDev);
	float yCoordFloat = (yCoord * yCoord) / (2 * stdDev * stdDev);
	return ampl * exp(-(xCoordFloat + yCoordFloat));
}

float uniform(float dummy_1, float dummy_2, int dummy_3, int dummy_4)
{
	return 1.0;
}

