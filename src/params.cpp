#include "params.h"

Params::Params(){}

Params::Params(const std::string inFile)
{
	std::ifstream iParamBuff(inFile.c_str());

	if (!iParamBuff.is_open())
	{
		std::cerr << "File '" << inFile << "' cannot be opened." << std::endl;
		// TODO: extend error class, throw it
		exit(1);
	}
	objectParams = json::parse(iParamBuff);

	iParamBuff.close();
}

std::string Params::toString()
{
	// if want to prettify, go ahead
	return objectParams.dump();
}

std::ostream &operator<<(std::ostream &os, Params &param)
{
	return os << param.toString();
}
