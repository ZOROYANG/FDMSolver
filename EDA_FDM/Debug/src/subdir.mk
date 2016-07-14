################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Configure.cpp \
../src/FDMSolver.cpp \
../src/main.cpp \
../src/matrix.cpp 

OBJS += \
./src/Configure.o \
./src/FDMSolver.o \
./src/main.o \
./src/matrix.o 

CPP_DEPS += \
./src/Configure.d \
./src/FDMSolver.d \
./src/main.d \
./src/matrix.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


