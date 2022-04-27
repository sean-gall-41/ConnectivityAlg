#ifndef CBMSTATE_H_
#define CBMSTATE_H_

#include <fstream>
#include <iostream>
#include <time.h>
#include <limits.h>

#include "params.h"

class State
{
public:
	State(std::fstream &actPFile, std::fstream &conPFile, std::fstream &simParams);
	~State();

private:
	Params *actParams;
	Params *conParams;
	Params *simParams;
};

#endif /* CBMSTATE_H_ */

