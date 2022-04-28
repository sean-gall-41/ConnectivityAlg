#include <assert.h>
#include "cellparams.h"
#include "celltype.h"
#include "conparams.h"
#include "connection.h"

const std::string ROOT 				    = "../../";
const std::string DATA_DIR 			    = ROOT + "data/";
const std::string INPUT_DIR 		    = DATA_DIR + "inputs/";
const std::string CELL_PARAM_FILE	    = INPUT_DIR + "cell_params.json";
const std::string CON_PARAM_FILE 	    = INPUT_DIR + "con_params_test_3.json";

float gaussian2d(float ampl, float stdDev, int xCoord, int yCoord);
float uniform(float dummy_1, float dummy_2, int dummy_3, int dummy_4);

int main()
{
	CellParams golgiParams(CELL_PARAM_FILE, "GO");
	CellType golgi(golgiParams);

	ConParams gogoConParams(CON_PARAM_FILE, "GOGO");
	
	Connection gogoConnection(golgi, golgi, gogoConParams);
	gogoConnection.establishConnectionDecay(20, 0, gaussian2d);

	if (CON_PARAM_FILE == "con_params_test_1.json")
	{
		assert(gogoConnection.getNumCons() == golgi.getNumTotCells());	
		assert(nearly_equal(gogoConnection.getAvgNumCons(), 1.0, 10000));
		assert(nearly_equal(gogoConnection.getRelRecipCons(), 1.0, 10000));
	}	
	else if (CON_PARAM_FILE == "con_params_test_2.json")
	{
		assert(gogoConnection.getNumCons() == golgi.getNumTotCells());	
		assert(nearly_equal(gogoConnection.getAvgNumCons(), 1.0, 10000));
		assert(nearly_equal(gogoConnection.getRelRecipCons(), 0.0, 10000));
	}	
	else /* for now */
	{
		std::cout << "Method 'establishConnectionDecay' results:" << std::endl; 
		std::cout << "Number of Golgi-Golgi connections: " 
			<< gogoConnection.getNumCons() << std::endl;
		std::cout << "Average Number of Golgi-Golgi connections: "
			<< gogoConnection.getAvgNumCons() << std::endl;
		std::cout << "Percent reciprocal Golgi-Golgi connections: "
			<< (gogoConnection.getRelRecipCons() * 100) << std::endl;
	}

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

