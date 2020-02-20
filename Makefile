CC=g++  -std=c++11
CFLAGS=-c -Wpedantic -Ofast
LDFLAGS=
SOURCES=cluster.cpp my_vector.cpp curve.cpp curve_point.cpp cluster_object.cpp ht.cpp h_funs.cpp g_funs.cpp NNpair.cpp grid.cpp curve_ht.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=cluster
TESTSOURCE=ourunit.cpp
TESTEXEC=ourunit


#VECTORS
#LSH
all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

test:
	$(CC) -o $(TESTEXEC) $(TESTSOURCE)

runtest:
	./$(TESTEXEC)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(TESTEXEC)
