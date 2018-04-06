CXXFLAGS= -Ofast -std=c++1z
CXXLINKFLAGS =  -Wl,--rpath -Wl,/home/empirelord/Dokumente/Devel/Packing/Ipopt-3.12.8/build/lib
INCL = `PKG_CONFIG_PATH=/home/empirelord/Dokumente/Devel/Packing/Ipopt-3.12.8/build/lib64/pkgconfig:/home/empirelord/Dokumente/Devel/Packing/Ipopt-3.12.8/build/lib/pkgconfig:/home/empirelord/Dokumente/Devel/Packing/Ipopt-3.12.8/build/share/pkgconfig: pkg-config --cflags ipopt`
LIBS = `PKG_CONFIG_PATH=/home/empirelord/Dokumente/Devel/Packing/Ipopt-3.12.8/build/lib64/pkgconfig:/home/empirelord/Dokumente/Devel/Packing/Ipopt-3.12.8/build/lib/pkgconfig:/home/empirelord/Dokumente/Devel/Packing/Ipopt-3.12.8/build/share/pkgconfig: pkg-config --libs ipopt`

release:
	g++ $(CXXFLAGS) $(INCL) $(CXXLINKFLAGS) $(LIBS) -o packer main.cpp
#	$(MAKE) clean

clean:
	rm -f *.o cubesolver
