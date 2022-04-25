BUILD_DIR ?= ./build
SRC_DIR ?= ./src
TARGET ?= $(BUILD_DIR)/connectivityAlg 

RAND_INCLUDE := ../random123/include
JSON_INCLUDE := ../json/include
SRCS := $(shell find $(SRC_DIR) -name *.cpp | xargs -I {} basename {})
OBJS := $(SRCS:%.cpp=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIR) -type d)
INC_DIRS += $(RAND_INCLUDE)
INC_DIRS += $(JSON_INCLUDE)
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CXX := g++ 
CPPFLAGS ?= -g -Wall -Wextra -std=c++11 -pedantic
CXXFLAGS ?= $(INC_FLAGS) -MMD -MP 

LD := g++
LD_FLAGS ?= -g -Wall -Wextra

MKDIR := mkdir -p
RM    := rm -rf

$(TARGET): $(OBJS)
	$(LD) $(LD_FLAGS) $(OBJS) -o $@ 

# Pattern rules
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp 
	$(MKDIR) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(BUILD_DIR)

-include $(DEPS)

