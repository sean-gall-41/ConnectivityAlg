#include "contypeparams.h"

ConTypeParams::ConTypeParams() {}

ConTypeParams::ConTypeParams(const std::string inFile) : Params(inFile)
{
	conTypeParamMap = objectParams.get<std::map<std::string, std::string>>();
}

std::string ConTypeParams::toString()
{
	std::string outString = "";
	for (auto i = conTypeParamMap.begin(); i != conTypeParamMap.end(); i++)
	{
		outString += "[ " + i->first + " : " + i->second + " ]\n";	
	}
	return outString;
}

std::string &ConTypeParams::operator[](const std::string key)
{
	return conTypeParamMap[key];
}

std::ostream &operator<<(std::ostream &os, ConTypeParams &conTypeParams)
{
	return os << conTypeParams.toString();
}
