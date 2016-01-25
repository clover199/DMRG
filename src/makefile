SOURCES = $(wildcard *.cpp)
OBJECTS = $(subst .cpp,.o,$(SOURCES))

all: $(OBJECTS)

%.o: %.cpp
	g++ -c $< -I ../include/ -I ../

clean:
	rm -rf *.o all
