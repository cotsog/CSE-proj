all:
	${CXX}  -std=c++0x -D__GXX_EXPERIMENTAL_CXX0X__ -O3 -Wall -c -fmessage-length=0 -o "src/WENO.o" "./src/WENO.cpp"	
	${CXX}  -std=c++0x -D__GXX_EXPERIMENTAL_CXX0X__ -O3 -Wall -c -fmessage-length=0 -o "src/cse-proj.o" "./src/cse-proj.cpp"
	${CXX} -o "cse-proj" ./src/WENO.o ./src/cse-proj.o -lm
	rm -rf ./src/*.o ./src/*.d
clean:
	rm -rf cse-proj
