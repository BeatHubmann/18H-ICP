CXX = g++-8
CXXFLAGS = -std=c++11

.PHONY: all
all: task1 task2 task3

task1.out: task1.cc
	$(CXX) $< -o $@

.PHONY: task1
task1: task1.out figure1.gplt
	./$< 3 31 1 200
	./$< 2836 127773 1 200
	gnuplot -e "input_2D='rand_num_c3_p31_2D.dat'; output_2D='figure1a.eps';\
	 input_3D='rand_num_c3_p31_3D.dat'; output_3D='figure1b.eps'" figure1.gplt
	gnuplot -e "input_2D='rand_num_c2836_p127773_2D.dat'; output_2D='figure1c.eps';\
	input_3D='rand_num_c2836_p127773_3D.dat'; output_3D='figure1d.eps'" figure1.gplt

task2.out: task2.cc
	$(CXX) $< -o $@

.PHONY: task2
task2: task2.out figure2.gplt
	./$< 2836 127773 1 200
	gnuplot figure2.gplt

task3.out: task3.cc
	$(CXX) $< -o $@

.PHONY: task3
task3: task3.out 

.PHONY: clean
clean:
	@rm -rf *.out *.dat *.eps
