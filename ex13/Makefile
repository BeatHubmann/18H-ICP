CXX = g++-8
CXXFLAGS = -std=c++17 -O3 -march=native -funroll-loops -ffast-math

.PHONY: all
all: main output movie 

main: main.cpp 

.PHONY: output
output: main
	./main 0.1 0.01 10 1000 

.PHONY: movie
movie: output plot.py
	python3 plot.py

.PHONY: clean
clean:
	@rm -rf main *.dat *.aux *.log *.synctex *.txt *.raw  