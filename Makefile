#makefile 

CC   =   g++
UCFLAGS = -O3 -fopenmp -Wall -gstabs+
RUCFLAGS := $(shell root-config --cflags) -I./include/ -I${CMSSW_BASE}/src/ -I${CMSSW_RELEASE_BASE}/src/ 
LIBS :=  -lgomp $(shell root-config --libs) -lTreePlayer  -lTMVA -lRooFit -lRooFitCore -L${CMSSW_BASE}/lib/${SCRAM_ARCH} -L${CMSSW_RELEASE_BASE}/lib/${SCRAM_ARCH} -lHiggsAnalysisGBRLikelihood -lCondFormatsEgammaObjects  
GLIBS := $(shell root-config --glibs)

vpath %.cpp ./src

SRCPP = main.cpp\
 	Utilities.cpp\
	SemiparametricGBRMaker.cpp

OBJCPP = $(patsubst %.cpp,obj/%.o,$(SRCPP))

all : regression.exe 
	#obj/libDictionary_C.so

obj/%.o : %.cpp
	@echo "> compiling $*"
	@mkdir -p obj/
	@$(CC) -c $< $(UCFLAGS) $(RUCFLAGS) -o $@

regression.exe : $(OBJCPP)
	@echo "> linking"
	@$(CC) $^ $(ACLIBS) $(LIBS) $(GLIBS)  -o $@

clean:
	@echo "> Cleaning object files"
	@rm  -f obj/*.o
        
cleanall: clean
	@echo "> Cleaning dictionary"
	@rm -f obj/libDictionary_C.so
	@echo "> Cleaning executable"
	@rm -f regression.exe

#obj/libDictionary_C.so: ./include/libDictionary.C
	#@echo "> Generating dictionary"
	#@cd include && root -b -q libDictionary.C++
	#@mv ./include/libDictionary_C.so ./obj/