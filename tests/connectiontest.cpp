#include <assert.h>
#include "cellparams.h"
#include "celltype.h"
#include "conparams.h"
#include "contypeparams.h"
#include "connection.h"


const std::string ROOT 				   = "../../";
const std::string DATA_DIR 			   = ROOT + "data/";
const std::string INPUT_DIR 		   = DATA_DIR + "inputs/";
const std::string OUTPUT_DIR		   = DATA_DIR + "outputs/";
const std::string CELL_PARAM_FILE	   = INPUT_DIR + "cell_params.json";
const std::string CON_PARAM_FILE 	   = INPUT_DIR + "con_params_test_1.json";
const std::string CON_TYPE_PARAM_FILE  = INPUT_DIR + "con_type_params.json";
const std::string CON_DECAY_OUT_FILE   = OUTPUT_DIR + "con_decay.bin";
const std::string CON_COMMON_OUT_FILE  = OUTPUT_DIR + "con_common.bin";


float gaussian2d(float ampl, float stdDev, int xCoord, int yCoord);
float uniform(float dummy_1, float dummy_2, int dummy_3, int dummy_4);


int main()
{
	CellParams golgiParams(CELL_PARAM_FILE, "GO");

	CellType golgi(golgiParams);
	std::cout << "=== Golgi Cell ===" << std::endl;
	std::cout << golgi << std::endl;	

	ConParams gogoConParams(CON_PARAM_FILE, "GOGO");
	std::cout << "=== gogo connection params ===" << std::endl;
	std::cout << gogoConParams << std::endl;

	ConTypeParams gogoConTypeParams(CON_TYPE_PARAM_FILE);
	std::cout << "=== gogo connection *type* params ===" << std::endl;
	std::cout << gogoConTypeParams << std::endl;

	Connection gogoConnection(golgi, golgi, gogoConParams);
	gogoConnection.establishConnectionDecay(
		std::stoi(gogoConTypeParams["maxConAttempts"]),
		std::stoi(gogoConTypeParams["randSeed"]),
		uniform);

	std::cout << "Method 'establishConnectionDecay' results:" << std::endl; 
	std::cout << "Number of Golgi-Golgi connections: " 
		<< gogoConnection.getNumCons() << std::endl;
	std::cout << "Average Number of Golgi-Golgi connections: "
		<< gogoConnection.getAvgNumCons() << std::endl;
	std::cout << "Percent reciprocal Golgi-Golgi connections: "
		<< (gogoConnection.getRelRecipCons() * 100) << std::endl;

	gogoConnection.toFile(CON_DECAY_OUT_FILE);

	gogoConnection.resetCons();	
	gogoConnection.establishConnectionCommon(
		(strcasecmp("true", gogoConTypeParams["needUnique"].c_str()) == 0),
		std::stoi(gogoConTypeParams["normConAttempts"]),
		std::stoi(gogoConTypeParams["maxConAttempts"]),
		std::stoi(gogoConTypeParams["randSeed"]));
	
	std::cout << std::endl;

	std::cout << "Method 'establishConnectionCommon' results:" << std::endl; 
	std::cout << "Number of Golgi-Golgi connections: " 
		<< gogoConnection.getNumCons() << std::endl;
	std::cout << "Average Number of Golgi-Golgi connections: "
		<< gogoConnection.getAvgNumCons() << std::endl;
	std::cout << "Percent reciprocal Golgi-Golgi connections: "
		<< (gogoConnection.getRelRecipCons() * 100) << std::endl;

	gogoConnection.toFile(CON_COMMON_OUT_FILE);

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

