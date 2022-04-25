#ifndef ASSERT_W_MSG_H_
#define ASSERT_W_MSG_H_

#include <iostream>

#define assertWMsg(msg, ...) do { \
	if(!(__VA_ARGS__)) { \
		std::cerr << msg << std::endl; \
	} \
} while(0)

#endif /* ASSERT_W_MSG_H_*/
