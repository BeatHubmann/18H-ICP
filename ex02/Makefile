CXX = g++-8
CXXFLAGS = -std=c++11 -O3 -march=native -funroll-loops

.PHONY: all
all: example

example: example.cpp
	$(CXX) $(CXXFLAGS) -o example.out $<


.PHONY: clean
clean:
	@rm *.out *.ppm *.png