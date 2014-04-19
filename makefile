all:
	# g++  -std=c++0x -fopenmp -D__GXX_EXPERIMENTAL_CXX0X__ -O3 -Wall -c -fmessage-length=0 -o "src/WENO.o" "./src/WENO.cpp"	
	# g++  -std=c++0x -D__GXX_EXPERIMENTAL_CXX0X__ -O3 -Wall -c -fmessage-length=0 -o "src/cse-proj.o" "./src/cse-proj.cpp"
	# g++ -o "cse-proj" ./src/WENO.o ./src/cse-proj.o -lm
	${CXX} -std=c++0x -fopenmp -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/Interpolation.d" -MT"src/Interpolation.d" -o "src/Interpolation.o" "./src/Interpolation.cpp"
	${CXX} -std=c++0x -fopenmp -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/WENO.d" -MT"src/WENO.d" -o "src/WENO.o" "./src/WENO.cpp"
	${CXX} -std=c++0x -fopenmp -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"src/main.d" -MT"src/main.d" -o "src/main.o" "./src/main.cpp"
	${CXX} -fopenmp  -o "WENO"  ./src/Interpolation.o ./src/WENO.o ./src/main.o   
clean:
	rm -rf  ./src/Interpolation.o ./src/WENO.o ./src/main.o  ./src/Interpolation.d ./src/WENO.d ./src/main.d  WENO
