#include "hash_table.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cassert>
#include <algorithm>


using namespace std;

const map<char, int> mapper::DNAmap = mapper::create_map();

Hash_table::Hash_table(int prefix_size_in, int kmer_size_in)
	:prefix_size(prefix_size_in), kmer_size(kmer_size_in)
{
	prefix.resize(1 << (prefix_size * 2), 0);
	suffix.resize(0);
}
Hash_table::Hash_table(const char* filename, int prefix_size_in, int kmer_size_in, bool fastq)
	:prefix_size(prefix_size_in), kmer_size(kmer_size_in)
{
	prefix.resize(1 << (prefix_size * 2), 0);
	suffix.resize(0);
	populate_from_file(filename, fastq);
}
long DNA2num(string kmer){
	long kmer_num = 0;
	for(unsigned int i = 0; i < kmer.size(); ++i){
		long shift = (kmer.size() - 1 - i);
		// crash if input is invalid
		assert(kmer[i] != 'T' || kmer[i] != 'A' || kmer[i] != 'C' || kmer[i] != 'G');
		kmer_num += ((unsigned long) mapper::DNAmap.at(kmer[i])) << (shift * 2);
	}
	return kmer_num;
}

bool Hash_table::hash(long kmer){
	return (kmer % 7 == 5 || kmer % 17 == 7 || kmer % 19 == 13 || kmer % 53 == 31 || kmer % 71 == 47);
}

void Hash_table::populate_from_file(const char* filename, bool fastq){
// add file not found exception
	ifstream infile(filename);
	if(!infile.fail()){
		long kmer_num;
		string kmer;
		string line;
		map <long, kmer_info> temp_map;
		int lineNum = 0;
		if(fastq)
			getline(infile, line);
		while(getline(infile, line)){
			bool skip = false;
			for(auto &x: line){
				if(x != 'A' && x != 'C' && x != 'G' && x != 'T' && x != 'a' && x != 'c' && x != 'g' && x != 't'){
					skip = true;
					break;
				}
			}
			if(!skip){
				for(unsigned int i = 0; i < line.length() - kmer_size + 1; ++i){
				
					kmer = line.substr(i, kmer_size);
				
					kmer_num = DNA2num(kmer);
					if(hash(kmer_num)){
						temp_map[kmer_num].kmer = kmer_num;
						temp_map[kmer_num].read_info.push_back(Read_Info(lineNum, i));
					}
				}
				++lineNum;
			}
			if(fastq){
				getline(infile, line);
				getline(infile, line);
				getline(infile, line);
			}
			
		}

		for(auto &x: temp_map){
			suffix.push_back(x.second);
		}
//		sort(suffix.begin(), suffix.end(), Kmer_comp());
		int num_kmers = 0;
		int num_individual = 0;
		for(unsigned long i = 0; i < suffix.size(); ++i){
			num_kmers += suffix[i].read_info.size();
			++num_individual;

			int prefix_index = DNA2num(num2DNA(suffix[i].kmer,kmer_size).substr(0, prefix_size));
			
			++prefix[prefix_index];
		}
		//cumulative sum of prefix
		int sum = 0;
		for(unsigned int i = 0; i < prefix.size(); ++i){
			sum += prefix[i];
			prefix[i] = sum;
		}
	}
	else{
		cout << "error: file not found" << endl;
		exit(1);
	}
}

string num2DNA(unsigned long num, int kmer_size){
	string DNA = "";
	unsigned long DNAindex;
	int count = 0;
	while(count < kmer_size){
		DNAindex = num & 0x3;
		DNA = DNAVALS[DNAindex] + DNA;
		num = num >> 2;
		++count;
	}
	return DNA;
}

int Hash_table::find_kmer_index(long kmer, int begin, int end){
	// binary search through hash table, find the kmer using kmer_value
	long value = kmer;
	int lower_bound = begin;
	int upper_bound = end;
	// basic binary search
	while(lower_bound <= upper_bound){
		long temp_comp_value = suffix[(upper_bound + lower_bound)/2].kmer;
		if(temp_comp_value > value){
			upper_bound = (upper_bound + lower_bound) / 2 - 1;
		}
		else if(temp_comp_value < value){
			lower_bound = (lower_bound + upper_bound) / 2 + 1;
		}
		// found a result
		else{
			// return index as a positive value if found (incremented by 1 to prevent it from
			// having a value of 0, which would be problematic
			return ((upper_bound + lower_bound)/2);
		}
	}
	return -1;
}
// returns kmer of -1 if the kmer is not found
kmer_info Hash_table::get_kmer_info(string kmer){
	assert(kmer.length() == (unsigned int) kmer_size);
	if(hash(DNA2num(kmer))){
		int prefix_index = DNA2num(kmer.substr(0, prefix_size));
		int endIndex = prefix[prefix_index];
		int beginIndex = 0;
		if(prefix_index != 0)
			beginIndex += prefix[prefix_index - 1];
		if(beginIndex != endIndex){
			int suffixIndex = find_kmer_index(DNA2num(kmer), beginIndex, endIndex);
			if(suffixIndex >= 0){
				return suffix[suffixIndex];
			}
		}
	}
	kmer_info temp;
	temp.kmer = -1;
	return temp;
}

bool Kmer_comp::operator()(kmer_info a, kmer_info b){
	return a.kmer < b.kmer;
}
