CC=g++
FLAGS=-std=c++20 -O3 -march=native

buchberger: main.cpp ../../lib/algo/buchberger.h ../../lib/algo/buchberger.cpp ../../lib/util/math.h ../../lib/util/math.cpp ../../lib/util/monomial.h ../../lib/util/monomial.cpp ../../lib/util/polynomial.h ../../lib/util/polynomial.cpp ../../lib/util/rational.h ../../lib/util/rational.cpp 
	$(CC) $(FLAGS) main.cpp ../../lib/algo/buchberger.h ../../lib/algo/buchberger.cpp ../../lib/util/math.h ../../lib/util/math.cpp ../../lib/util/monomial.h ../../lib/util/monomial.cpp ../../lib/util/polynomial.h ../../lib/util/polynomial.cpp ../../lib/util/rational.h ../../lib/util/rational.cpp -o ../bin/benchmark_buchberger