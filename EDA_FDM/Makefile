lib = -lumfpack -lopenblas -pthread\
      -lamd  -lcholmod -lcxsparse  \
      -lsuitesparseconfig -lcolamd \
      -lcamd -lccolamd -lmetis -lrt
sub = main.o matrix.o FDMSolver.o Configure.o Debugger.o
subjects = $(sub:%.o=obj/%.o)

vpath = include:src:obj
vpath %.cpp src
vpath %.h include
vpath %.o obj
OBJDIR = ./obj

path_lib=/home/ming/lib/zhangbl/lib
path_include=/home/ming/lib/zhangbl/include

CXX = g++
FLAGS = -std=c++11 -O3 -Wall
strip = strip -g --strip-unneeded bin/FRW_FDM

bin/FRW_FDM:$(sub)
	@echo Linking FRW_FDM
	@$(CXX) $(FLAGS) -o bin/FRW_FDM $(subjects) -static -L$(path_lib) $(lib)
	@$(strip)
	
%.o:%.cpp
	@echo Compiling $(@F)...
	@$(CXX) $(FLAGS) -c -o $(OBJDIR)/$@ -I$(path_include) $< -I include
.PHONY : clean
clean:
	-rm -f $(subjects)
	