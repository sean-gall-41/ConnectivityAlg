#ifndef CELL_PARAMS_H_
#define CELL_PARAMS_H_

#include "params.h"

class CellParams : Params
{
public:
	CellParams(){}
	CellParams(const std::string inFile, std::string cellID) : Params(inFile)
	{
		cellParamMap = objectParams[cellID].get<std::map<std::string, std::string>>();	
	}

	std::string &operator[](const std::string key) {return cellParamMap[key];}
	// const unsigned int &operator[](const std::string key) const {return cellParamMap[key];}

private:
	std::map<std::string, std::string> cellParamMap;
};

#endif /* CELL_PARAMS_H_ */

