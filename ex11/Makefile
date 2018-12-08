CXX = g++-8
CXXFLAGS = -std=c++1z -O3 -march=native -funroll-loops -ffast-math
CXXFLAGS += -I/usr/local/include/eigen3

.PHONY: all
all: main output figures 

main: main.cpp 

.PHONY: output
output: main
	./main 50 0.25 0.75 

.PHONY: figures
figures: output plot.py
	python3 plot.py

.PHONY: clean
clean:
	@rm -rf main *.dat *.aux *.log *.synctex *.txt 