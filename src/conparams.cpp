#include "conparams.h"

ConParams::ConParams() {};
ConParams::ConParams(const std::string inFile, std::string conID) : Params(inFile)
{
	// usage of the array subscript operator for json objects 
	// and usage of json's get method to transform json -> map object
	conParamMap = objectParams[conID].get<std::map<std::string, std::string>>();

}

std::string ConParams::toString()
{
	std::string outString = "";
	for (auto i = conParamMap.begin(); i != conParamMap.end(); i++)
	{
		outString += "[ " + i->first + " : " + i->second + " ]\n";	
	}
	return outString;
}

std::ostream &operator<<(std::ostream &os, ConParams &conParams)
{
	return os << conParams.toString();
}


