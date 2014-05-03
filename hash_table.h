#ifndef HASH_TABLE_H
#define HASH_TABLE_H
#include <vector>
#include <iostream>
#include <map>

using namespace std;

// opposite of the nucleotide mapper
const char DNAVALS[] = {'A', 'C', 'G', 'T'};

// takes a nucleotide (char) input, and returns the value of the nucleotide
// A = 0, C = 1, G = 2, T = 3
struct mapper{
	static std::map<char, int> create_map(){
		std::map<char, int> m;
		m.insert(pair<char, int>('A', 0));
		m.insert(pair<char, int>('C', 1));
		m.insert(pair<char, int>('G', 2));
		m.insert(pair<char, int>('T', 3));
		m.insert(pair<char, int>('a', 0));
		m.insert(pair<char, int>('c', 1));
		m.insert(pair<char, int>('g', 2));
		m.insert(pair<char, int>('t', 3));
		return m;
	}
	static const std::map<char, int> DNAmap;
};

// converts the kmer into an integer
long DNA2num(string kmer);

// convers an integer into its corresponding kmer
string num2DNA(unsigned long num, int kmer_size);

struct Read_Info{
	Read_Info(int read_in, int read_loc_in)
		:read(read_in), readLoc(read_loc_in){}
	int read;
	int readLoc;
};
// stored in the suffix hash table, contains the kmer (as an int)
// and the locations where the kmer was found
struct kmer_info{
	long kmer;
	vector<Read_Info> read_info;
};

class Hash_table{
  private:

  public:
	// stores the prefix size and kmer size
	int prefix_size;
	int kmer_size;

	vector<int> prefix;
	vector<kmer_info> suffix;

	// only numbers that pass the hash function are placed into the hash
	// table
	bool hash(long kmer);

	// finds the index of the kmer
	int find_kmer_index(long kmer, int begin, int end);

	// constructor: default prefix size is 5, default kmer size is 20
	Hash_table(int prefix_size, int kmer_size);

	Hash_table(const char* filename, int prefix_size_in, int kmer_size_in, bool fastq);

	// inputs the kmers from a file
	void populate_from_file(const char* filename, bool fastq);

	// get information from the table about a particular kmer
	kmer_info get_kmer_info(string kmer);
};

struct Kmer_comp{
	bool operator()(kmer_info a, kmer_info b);
};
#endif /* HASH_TABLE_H ENDIF */
