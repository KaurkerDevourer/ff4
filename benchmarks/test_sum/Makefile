CC=g++
FLAGS=-std=c++20 -O3 -march=native

benchmake: main.cpp ../benchmarking.h 
	$(CC) $(FLAGS)  main.cpp ../benchmarking.h -o ../bin/benchmark_sum