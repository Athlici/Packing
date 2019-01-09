CXXFLAGS= -std=c++1z -Ofast

release: main.o tinyxml2.o
	g++ $(CXXFLAGS) -lipopt -o packer tinyxml2.o main.o

debug: main.o tinyxml2.o
	g++ $(CXXFLAGS) -ggdb -lipopt -o packer tinyxml2.o main.o

examplegen: tinyxml2.o
	g++ -o examplegen tinyxml2.o examplegen.cpp

main.o: *.cpp
	g++ $(CXXFLAGS) -c main.cpp

tinyxml2.o: tinyxml2.cpp
	g++ $(CXXFLAGS) -c tinyxml2.cpp

clean:
	rm -f *.o packer
