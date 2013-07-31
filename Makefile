CC=clang++
CFLAGS= -std=c++11 -stdlib=libc++ -Wall
LDFLAGS= -c
SOURCES=tests/utils.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test

run: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@	
	
clean:
	rm -rf *o test