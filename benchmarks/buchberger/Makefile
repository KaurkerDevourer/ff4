CC=g++
FLAGS=-std=c++20 -O3 -march=native

buchberger: main.cpp ../../lib/algo/buchberger.h ../../lib/algo/buchberger.cpp ../../lib/util/monomial.h ../../lib/util/monomial.cpp ../../lib/util/polynomial.h ../../lib/util/polynomial.cpp ../../lib/util/rational.h ../../lib/util/rational.cpp ../../lib/util/term.h ../../lib/util/term.cpp ../../GroebnerBasis/GroebnerLib/srcs/Lib.cpp ../../GroebnerBasis/GroebnerLib/srcs/Rational.cpp ../../GroebnerBasis/GroebnerLib/srcs/Term.cpp
	$(CC) $(FLAGS) main.cpp ../../lib/algo/buchberger.h ../../lib/algo/buchberger.cpp ../../lib/util/monomial.h ../../lib/util/monomial.cpp ../../lib/util/term.h ../../lib/util/term.cpp ../../lib/util/polynomial.h ../../lib/util/polynomial.cpp ../../lib/util/rational.h ../../lib/util/rational.cpp  -I ../../GroebnerBasis/GroebnerLib/includes/ ../../GroebnerBasis/GroebnerLib/srcs/Lib.cpp ../../GroebnerBasis/GroebnerLib/srcs/Rational.cpp ../../GroebnerBasis/GroebnerLib/srcs/Term.cpp -o ../bin/benchmark_buchberger
