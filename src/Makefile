#CC=clang++
CC=g++
CPPOPTS=-std=c++17
PERF=-O3 -mavx2
#PERF=-ggdb

all : omega3d

omega3d : main.cpp
	$(CC) $(CPPOPTS) $(PERF) -o $@ $<

clean : 
	rm -f omega3d
