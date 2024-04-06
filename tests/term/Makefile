CC=g++
FLAGS=-std=c++20 -g -rdynamic -march=native

ational: main.cpp ../../lib/util/term.h ../../lib/util/term.cpp ../../lib/util/rational.h ../../lib/util/rational.cpp ../testing.h
	$(CC) $(FLAGS) main.cpp ../../lib/util/term.h ../../lib/util/term.cpp ../../lib/util/rational.h ../../lib/util/rational.cpp ../testing.h -o ../bin/test_term
