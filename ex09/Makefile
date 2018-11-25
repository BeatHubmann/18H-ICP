CXX = g++-8
CXXFLAGS = -std=c++1z -O3 -march=native -funroll-loops -ffast-math
CXXFLAGS += -I/usr/local/include/eigen3


.PHONY: all
all: main tasks.dat figures 

main: main.cpp 

tasks.dat: main
	./main 2 2 1.7 2.3 > tasks.dat

.PHONY: figures
figures: tasks.dat figures.gplt
	gnuplot figures.gplt

.PHONY: clean
clean:
	@rm -rf main *.dat *.aux *.log *.synctex 