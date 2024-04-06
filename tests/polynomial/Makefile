CC=g++
FLAGS=-std=c++20 -g -rdynamic -march=native

ational: main.cpp ../../lib/util/term.h ../../lib/util/term.cpp ../../lib/util/rational.h ../../lib/util/rational.cpp  ../../lib/util/monomial.h  ../../lib/util/monomial.cpp ../../lib/util/polynomial.h  ../../lib/util/polynomial.cpp ../testing.h
	$(CC) $(FLAGS) main.cpp ../../lib/util/term.h ../../lib/util/term.cpp ../../lib/util/rational.h ../../lib/util/rational.cpp ../../lib/util/monomial.h  ../../lib/util/monomial.cpp ../../lib/util/polynomial.h  ../../lib/util/polynomial.cpp ../testing.h -o ../bin/test_polynomial
