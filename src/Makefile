#CC=clang++
CC=g++
CPPOPTS=-std=c++17
PERF=-O3 -mavx2
#PERF=-ggdb

all : Omega3D.bin

Omega3D.bin : main.cpp
	$(CC) $(CPPOPTS) $(PERF) -o $@ $<

clean : 
	rm -f Omega3D.bin
