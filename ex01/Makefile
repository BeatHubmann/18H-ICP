CXX=g++-8
CXXFLAGS = -std=c++11

.PHONY: all
all: task1 task2 task3

task2.out: task2.cc
	$(CXX) $< -o $@

.PHONY: task2
task2: task2.out figure2.gplt
	./$< 2836 127773 1 200
	gnuplot figure2.gplt

.PHONY: clean
clean:
	@rm -rf *.out *.dat *.eps
