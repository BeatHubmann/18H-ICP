CXX = g++-8
CXXFLAGS = -std=c++1z -O3 -march=native -funroll-loops 

.PHONY: all
all: main 

main: main.cpp latticeview.h

.PHONY: clean
clean:
	@rm -rf task *.dat *.eps *.aux *.log *.synctex *.ppm 