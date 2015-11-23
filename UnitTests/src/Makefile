CC=g++
CFLAGS=-c 
LDFLAGS= -g3 -I ../include -I ../../eigen/ -L../
SOURCES=../../Polynomial.cpp ../../SylvesterMatrix.cpp ../../Solution.cpp ../../BivariatePolynomial.cpp ../../ProblemSolver.cpp ../../SylvesterPolynomial.cpp ../../Matrix.cpp main.cpp
OBJECTS=$(SOURCES:.h.c=.o)
EXECUTABLE=exeTests

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@  -llapacke -lgtest -lpthread

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm exeTests

