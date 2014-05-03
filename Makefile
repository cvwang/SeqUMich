PATH := /usr/um/gcc-4.7.0/bin:$(PATH)
LD_LIBRARY_PATH := /usr/um/gcc-4.7.0/lib64
LD_RUN_PATH := /usr/um/gcc-4.7.0/lib64

FLAGS = -Wall -pedantic -std=c++11

all: main.o hash_table.o
	g++ $(FLAGS) main.o hash_table.o -o main

main.o: main.cpp
	g++ $(FLAGS) -c main.cpp

hash_table.o: hash_table.cpp
	g++ $(FLAGS) -c hash_table.cpp

clean: 
	rm -f *.o main
