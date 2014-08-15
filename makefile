################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# All of the sources participating in the build are defined here
include subdir.mk
include ./sofia-ml/sofia.mk

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: CNFalign

# Tool invocations
CNFalign: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
#	icpc -axS -o "CNFpred" -fp-model fast=2 -march=core2  -fast $(OBJS) $(USER_OBJS) $(LIBS)
	g++ -static -o "CNFalign" -fno-strict-aliasing -O3 -Wall $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '
	#-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS)

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS)
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
