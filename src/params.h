#ifndef PARAMS_H_
#define PARAMS_H_

#include <string>
#include <map>
#include <fstream>
#include <sstream> 
#include <iostream>
#include <cstdlib>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class Params
{
public:
	Params();
	Params(const std::string inFile);
	
	std::string toString();	

	friend std::ostream &operator<<(std::ostream &os, Params &param);

protected:
	json objectParams;
};

#endif /* PARAMS_H_ */
