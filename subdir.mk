################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
./ScoreMatrix.cpp \
./CNFalign_Basic.cpp \
./CNFalign_Feat.cpp \
./CNFalign_Calc.cpp \
./getopt.cpp \
./mt19937ar.cpp \
./seq.cpp \
./template.cpp \
./profile.cpp \
./hhpred_util.cpp \
./hhpred_hmm.cpp \
./hhpred_hit.cpp \
./CNFalign_Util.cpp \
./CNF_constants.cpp \
./hcsearch.cpp 
#./CNFalign.cpp 

OBJS += \
./ScoreMatrix.o \
./CNFalign_Basic.o \
./CNFalign_Feat.o \
./CNFalign_Calc.o \
./getopt.o \
./mt19937ar.o \
./seq.o \
./template.o \
./profile.o \
./hhpred_util.o \
./hhpred_hmm.o \
./hhpred_hit.o \
./CNFalign_Util.o \
./CNF_constants.o \
./hcsearch.o
#./CNFalign.o

CPP_DEPS += \
./ScoreMatrix.d \
./CNFalign_Basic.d \
./CNFalign_Feat.d \
./CNFalign_Calc.d \
./mt19937ar.d \
./getopt.d \
./seq.d \
./template.d \
./profile.d \
./hhpred_util.d \
./hhpred_hmm.d \
./hhpred_hit.d \
./CNFalign_Util.d \
./CNF_constants.d \
./hcsearch.o \
./hcsearch.d
#./CNFalign.d

# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.cpp
	@echo 'Building file: $<'
	g++ -O3 -static -ffast-math -c -fmessage-length=0 -fno-strict-aliasing -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


