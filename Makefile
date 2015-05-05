# make file to include stuff 
# in MY_CPP 
#Inc= -I/home/walter/Desktop/4Simulation/My_CPP
CC=g++
OBJs=InitialGeometry.o  DumpMesh.o

default: DumpMesh
DumpMesh.o: DumpMesh.C
	$(CC) -c $< -o $@ 
InitialGeometry.o: InitialGeometry.cpp
	$(CC) -c $< -o $@ 
DumpMesh: $(OBJs)
	$(CC)  $^ -o $@
test: test.C
	$(CC)  $^ -o $@ $(Inc)
.PHONEY:clean

clean:
	rm -f DumpMesh *.o
	
