SOURCES = $(wildcard *.cpp)
OBJECTS = $(subst .cpp,.o,$(SOURCES))

all: $(OBJECTS)

%.o: %.cpp
	g++ -c $< -I ../include/

clean:
	rm -rf *.o all
