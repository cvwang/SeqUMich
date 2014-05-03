#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <cassert>
#include <string>
#include <map>
#include <cstdlib>
#include <cstdio>
#include "hash_table.h"
using namespace std;

// takes a nucleotide (char) input, and returns the value of the nucleotide
// A = 0, C = 1, G = 2, T = 3
bool hasher(long kmer){
	return (kmer % 7 == 5 || kmer % 17 == 7 || kmer % 19 == 13 || kmer % 53 == 31 || kmer % 71 == 47);
}
int main(){
    string input_filename = "SMALL_target_string.txt";
    string reads_filename = "SMALL_reads.txt";
    ifstream reads_input(reads_filename.c_str());
    ifstream input_target_file(input_filename.c_str());
    unordered_map<string, int> kmer_totals;
    
    int kmer_size = 16;
	int counter = 0;
	string target = "";
    while(!reads_input.eof()){
        string line;
        getline(reads_input, line);
        for(int i = 0; i < line.length() - kmer_size + 1; ++i){
            string s = line.substr(i, kmer_size);
            if(hasher(DNA2num(s))){
                kmer_totals[s]++;
            }
        }
    }
    counter = 0;
    
	while(input_target_file.peek() != '\n' && !input_target_file.eof()){
        if(input_target_file.eof() == true || input_target_file.peek() == '\n') break;
		// cout << "here1" << endl;
		if(target.length() >= kmer_size){
            string sKmer = target.substr(counter, kmer_size);
            cout << "sKmer: " << sKmer << endl;
            for(auto &x: kmer_totals){
                if(x.first == sKmer){
                    ++counter;
                }
            }
		}
        
		// counter++;
		// if(counter == 40)
		// 	break;
	}
    cout << "string " << target << endl;
    cout << "found: " << counter << endl;
}