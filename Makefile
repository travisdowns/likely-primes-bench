CPPFLAGS=-g -O3 -march=native -std=c++11 -Wall -Wno-unused

all: likely-primes

likely-primes: main.o process.o tables512.o
	g++ $^ -o likely-primes
 
process.o: process.asm Makefile
	nasm -w+all -f elf64 process.asm
	
main.o : main.cpp tables.hpp tables256.hpp Makefile

clean:
	rm *.o likely-primes
