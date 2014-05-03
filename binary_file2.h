#ifndef _binary_file2_h
#define _binary_file2_h
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <istream>
#include <vector>
#include <assert.h>

using namespace std;

unsigned char* read2char(string read, int max_read_length);
char* binary_filename;

string convert(unsigned char *kmers){
	int size = (int)kmers[0];
	string result = "";
	for(int i = 0; i < size; ++i){
		int char_index = i / 4 + 1;
		int char_shift = 3 - i % 4;
		int DNA_val = (kmers[char_index] >> (char_shift * 2)) & 0x3;
        result.push_back(DNAVALS[DNA_val]);
	}
	return result;
}

void set_binary_file(const char *bin_name){
  binary_filename = (char*) bin_name;
}

void make_binary_file(const char *filename, int max_read_length, bool fastq){
  ifstream infile(filename);
  ofstream outfile;
  outfile.open(binary_filename, ios::out | ios::binary);
  if(fastq)
	infile.ignore(1000,'\n');
  int number_chars = (int) ceil(max_read_length/4.0);
  while(!infile.eof()){
	bool skip = false;
	string read;
	infile >> read;
	int size = read.size();
	if(fastq)
	  for(int i=0;i<4;i++)
		infile.ignore(1000,'\n');
	for(int i=0;i<size;i++)
	  if(read[i] != 'A' && read[i] != 'C' && read[i] != 'G' && read[i] != 'T' && read[i] != 'a' && read[i] != 'c' && read[i] != 'g' && read[i] != 't')
		skip = true;
	if(!skip){
	  unsigned char* data = read2char(read, number_chars);
        outfile.write((char*)data, number_chars + 1);
	  delete data;
	}
  }
  infile.close();
  outfile.close();
}

unsigned char* read2char(string read, int number_chars){
  unsigned char* return_char = new unsigned char[number_chars + 1];
    return_char[0] = read.size();
  unsigned char* read_nucleotides = new unsigned char[number_chars*4];
  for(int i=0;i<read.size();++i){
    read_nucleotides[i] = read[i];
  }
  for(int i=read.size();i<number_chars*4;++i){
    read_nucleotides[i] = 'A';
  }
  for(int i=0;i<number_chars;++i){
    int char_value = 0;
    for(int j=0;j<4;++j){
      char_value += (mapper::DNAmap.at(read_nucleotides[4*i+j]) << 2*(3-j));
    }
    return_char[i + 1] = (unsigned char)(char_value);
  }
  return return_char;
}

unsigned char* find_read(int read, int max_read_length){
  int number_chars =(int) ceil(max_read_length/4.0);
  char * memblock = new char[number_chars + 1];
  ifstream infile(binary_filename, ios::out | ios::binary);
  infile.seekg((number_chars + 1) * read);
  infile.read(memblock, (streampos) (number_chars + 1));
  infile.close();
  return (unsigned char*)memblock;
}
#endif
