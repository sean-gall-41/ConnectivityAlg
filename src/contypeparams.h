#ifndef CONTYPEPARAMS_H_
#define CONTYPEPARAMS_H_

#include "params.h"

class ConTypeParams : public Params
{
public:
	ConTypeParams();
	ConTypeParams(const std::string inFile);

	std::string toString();

	std::string &operator[](const std::string key);
	friend std::ostream &operator<<(std::ostream &os, ConTypeParams &conTypeParams);

protected:
	std::map<std::string, std::string> conTypeParamMap;
};

#endif /* CONTYPEPARAMS_H_ */

