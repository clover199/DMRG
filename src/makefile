SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES, .cpp=.o)

all: $(OBJECTS)

%.o: %.cpp
	$(CC) $(FLAGS) -c $< -I ../include/ -I ../

clean:
	rm -rf *.o all
