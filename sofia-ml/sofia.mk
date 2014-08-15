################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
./sofia-ml/sf-data-set.cc \
./sofia-ml/sf-hash-inline.cc \
./sofia-ml/sf-hash-weight-vector.cc \
./sofia-ml/sf-sparse-vector.cc \
./sofia-ml/sf-weight-vector.cc \
./sofia-ml/sofia-ml-methods.cc

OBJS += \
./sofia-ml/sf-data-set.o \
./sofia-ml/sf-hash-inline.o \
./sofia-ml/sf-hash-weight-vector.o \
./sofia-ml/sf-sparse-vector.o \
./sofia-ml/sf-weight-vector.o \
./sofia-ml/sofia-ml-methods.o

CPP_DEPS += \
./sofia-ml/sf-data-set.d \
./sofia-ml/sf-hash-inline.d \
./sofia-ml/sf-hash-weight-vector.d \
./sofia-ml/sf-sparse-vector.d \
./sofia-ml/sf-weight-vector.d \
./sofia-ml/sofia-ml-methods.d

# Each subdirectory must supply rules for building sources it contributes
%.o: ./%.cc
	@echo 'Building file: $<'
	g++ -O3 -static -ffast-math -c -fmessage-length=0 -fno-strict-aliasing -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '
