CXX = g++-8
CXXFLAGS = -std=c++1z -O3 -march=native -funroll-loops -fopenmp

.PHONY: all
all: example task1 image2 task3

example: example.cpp

task1: task1.cpp

figure1: task1
	./task1
	convert task1.png task1.eps

task2: task2.cpp

figure2: task2
	./task2 0.6 100 10 1
	convert task2_*.png task2.gif
	convert task2_002.png task2_002.eps
	convert task2_100.png task2_100.eps

task3: task3.cpp

figure3: task3
	./task3  20   10 42 1000 >   task3_10.dat
	./task3  20   25 42 1000 >   task3_25.dat
	./task3  20   50 42 1000 >   task3_50.dat
	./task3  20  100 42 1000 >  task3_100.dat
	./task3  20  250 42 1000 >  task3_250.dat
	./task3  20  500 42 1000 >  task3_500.dat
	./task3  20 1000 42 1000 > task3_1000.dat
	gnuplot task3.gplt

.PHONY: clean
clean:
	@rm -rf example task1 task2 task3 *.ppm *.png *.gif *.dat *.eps