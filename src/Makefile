#CC=clang++
CC=g++
CXXOPTS=-std=c++17
INCLUDE=
LIBS=
PERF=-O3 -mavx2
#PERF=-ggdb

CXXOPTS+=-DUSE_VC
INCLUDE+=-I/opt/Vc/include
LIBS+=-L/opt/Vc/lib -lVc
PERF+=-mfma
CXXOPTS+=-fopenmp

all : Omega3D.bin

Omega3D.bin : main.cpp
	$(CC) $(CXXOPTS) $(PERF) $(INCLUDE) -o $@ $< $(LIBS)

clean : 
	rm -f Omega3D.bin
