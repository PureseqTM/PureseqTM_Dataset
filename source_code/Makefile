

# make the subdir first
SHELL :=/bin/bash
CXX = g++
FLAG = -O3

all:
	#----- for 2-state label -----#
	$(CXX) GFF_to_TransMemb.cpp $(FLAG) -o GFF_to_TransMemb
	#----- for evaluation --#
	$(CXX) TM2_to_Evaluation.cpp $(FLAG) -o TM2_to_Evaluation
	#-> move to ../util
	mv GFF_to_TransMemb TM2_to_Evaluation ../util

# clean all the temp files
clean:

