CXXFLAGS= -Ofast -std=c++1z

release:
	g++ $(CXXFLAGS) -lipopt -o packer main.cpp

clean:
	rm -f *.o packer
