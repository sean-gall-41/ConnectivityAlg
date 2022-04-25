#ifndef CON_PARAMS_H_
#define CON_PARAMS_H_

#include "params.h"

class ConParams : public Params
{
public:
	ConParams();
	ConParams(const std::string inFile, std::string conID);

	std::string toString();
	friend std::ostream &operator<<(std::ostream &os, ConParams &conParams);

	std::string &operator[](const std::string key) {return conParamMap[key];}
	// const unsigned int &operator[](const std::string key) const {return conParamMap[key];}

private:
		std::map<std::string, std::string> conParamMap;
};

#endif /* CON_PARAMS_H_ */
