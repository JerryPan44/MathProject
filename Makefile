CC=g++
CFLAGS=-c
LDFLAGS= -g3 -I eigen/
SOURCES=main.cpp Parse.cpp Solution.cpp BivariatePolynomial.cpp SylvesterMatrix.cpp Polynomial.cpp SylvesterPolynomial.cpp Matrix.cpp ProblemSolver.cpp Interpolation.cpp
OBJECTS=$(SOURCES:.h.c=.o)
EXECUTABLE=exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ -llapacke

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm exe




