CC=g++

all:
	g++  -std=c++0x -D__GXX_EXPERIMENTAL_CXX0X__ -O3 -Wall -c -fmessage-length=0 -o "src/WENO.o" "./src/WENO.cpp"	
	g++  -std=c++0x -D__GXX_EXPERIMENTAL_CXX0X__ -O3 -Wall -c -fmessage-length=0 -o "src/cse-proj.o" "./src/cse-proj.cpp"
	g++ -o "cse-proj" ./src/WENO.o ./src/cse-proj.o -lm
	rm -rf ./src/*.o ./src/*.d
clean:
	rm -rf cse-proj
