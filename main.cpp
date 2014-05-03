#include <iostream>
#include "hash_table.h"
#include "lu_wang.h"
#include "binary_file2.h"
#include <cstring>
#include <string>
#include <fstream>
#include <vector>
#include "getopt.h"
#include <ctime>
const int MATCH = 1;
const int MISMATCH = -1;
const int DEL = -1;
const int INS = -1;
const double PERCENT = 0;

using namespace std;
string convert(unsigned char *kmers);
unsigned char* process_read(unsigned char* read, int read_loc, int kmer_size, int read_buffer);
void printHelp();
bool checkNumeric(char *c_str);
string to_upper(string x);

int main(int argc, char* argv[]){
	struct option longOpts[] = {
		{"debug", no_argument, NULL, 'd'},
		{"target_string_filename", required_argument, NULL, 't'},
		{"read_filename", required_argument, NULL, 'r'},
		{"output_binary_filename", required_argument, NULL, 'o'},
		{"help", no_argument, NULL, 'h'},
		{"kmer_size", required_argument, NULL, 'k'},
		{"prefix_size", required_argument, NULL, 'p'},
		{"max_read_size", required_argument, NULL, 'm'},
		{"buffer", required_argument, NULL, 'b'},
        {"input_binary_filename", required_argument, NULL, 'i'}
	};
	int opt = 0, index = 0;
	string input_read_filename = "";
	string input_target_filename = "";
	string output_read_filename = "";
	int max_read_length = 150;
	int prefix_length = 4;
	int read_buffer = 10;
	int kmer_length = 16;
	int required_args = 0;
	bool debug = false;
    bool make_binary = true;
	while((opt = getopt_long (argc, argv, "dt:r:o:hk:p:m:b:i:", longOpts, &index)) != -1){
		switch(opt){
			case 'd':
				debug = true;
				break;
			case 't':
				if(input_target_filename == "" && strcmp(optarg, "")){
					input_target_filename = optarg;
					++required_args;
				}
				else{
					cout << "no argument for input target filename or input target filename specified twice. run with -h tag to display help\n";
					exit(1);
				}
				break;
			case 'r':
				if(input_read_filename == "" && strcmp(optarg, "")){
					input_read_filename = optarg;
					++required_args;
				}
				else{
					cout << "no argument for input read filename or input read filename specified twice. run with -h tag to display help\n";
					exit(1);
				}
				break;
			case 'o':
				if(output_read_filename == "" && strcmp(optarg, "")){
					output_read_filename = optarg;
					++required_args;
				}
				else{
					cout << "no argument for output binary filename or output binary filename specified twice. run with -h tag to display help\n";
					exit(1);
				}
				break;
			case 'h':
				printHelp();
				exit(0);
			case 'k':
				if(checkNumeric(optarg)){
					kmer_length = atoi(optarg);
				}
				else {
					cout << "kmer size was not numeric!\n";
					exit(0);
				}
				break;
			case 'p':
				if(checkNumeric(optarg)){
					prefix_length = atoi(optarg);
				}
				else{
					cout << "prefix size was not numeric!\n";
					exit(0);
				}
				break;
			case 'm':
				if(checkNumeric(optarg)){
					max_read_length = atoi(optarg);
					++required_args;
				}
				else{
					cout << "max read length is not numeric!\n";
					exit(0);
				}
				break;
			case 'b':
				if(checkNumeric(optarg)){
					read_buffer = atoi(optarg);
				}
				else {
					cout << "buffer size is not numeric!\n";
					exit(0);
				}
				break;
            case 'i':
                make_binary = false;
                if(output_read_filename == "" && strcmp(optarg, "")){
					output_read_filename = optarg;
					++required_args;
				}
				else{
					cout << "no argument for output binary filename or output binary filename specified twice. run with -h tag to display help\n";
					exit(1);
				}
				break;
			default:
				cout << "unknown tag: " << opt << "\n";
				printHelp();
				exit(0);
		}
	}
	
	if(required_args < 4){
		cout << "not all required arguments were met. here are all the required arguments\n";
		printHelp();
	}
	bool fastq = input_read_filename.substr(input_read_filename.size() - 6) == ".fastq";
//	cout << fastq << endl;
	auto start = clock();
	set_binary_file(output_read_filename.c_str());
	// if binary file is already created....we don't have to do this
    if(make_binary)
        make_binary_file(input_read_filename.c_str(), max_read_length, fastq);

    double total_bin_time = (clock() - start) / (double)CLOCKS_PER_SEC;


    start = clock();

	Hash_table kmer_table(input_read_filename.c_str(), prefix_length, kmer_length, fastq);

	double total_index_time = (clock() - start) / (double)CLOCKS_PER_SEC;
	// // take each kmer from input target & find them the hash
	ifstream input_target_file(input_target_filename.c_str());

	// get all kmers from the input target string and get their information
	// CHARLES: YOU ARE MISSING THE FIRST KMER HERE
	start = clock();
	int print_counter = 0;
	int target_number = 0;
	string target_string = "";

	getline(input_target_file, target_string);

	while(!input_target_file.eof()){
		target_string = "";
		int counter = 0;
		// while(input_target_file.peek() != '\n'){
		// 	target_string.push_back(input_target_file.get());
		// 	if(target_string.length() >= kmer_length){
		// 		string sKmer = target_string.substr(counter, kmer_length);
		// 		// only if it passes the hash function!
		// 		kmer_info temp = kmer_table.get_kmer_info(sKmer);
		// 		if(temp.kmer != -1){
		// 		// cout << num2DNA(temp.kmer, kmer_length) << endl;
		// 			all_kmer_info.push_back(temp);
		// 			target_loc.push_back(counter);
		// 		}
		// 		++counter;
		// 	}
		// }
		// input_target_file.get();

		//concatenate the strings of the .fa file
		//stringstream target_stream;
		string waste_string = "";
		while(input_target_file.peek() != '>' && !input_target_file.eof()){
			string temp = "";
			input_target_file >> temp;
			target_string = target_string + temp;
			getline(input_target_file, waste_string);
		}
		
		getline(input_target_file, waste_string);

		// cout << target_string << endl;
		// 
		// cout << "here" << endl;
		// exit(1);


		vector<kmer_info> all_kmer_info;
		vector<int> target_loc;
		for(int i=0; i<target_string.size()-kmer_length+1; i++){
			string sKmer = target_string.substr(i, kmer_length);
			// only if it passes the hash function!
			kmer_info temp = kmer_table.get_kmer_info(sKmer);
			if(temp.kmer != -1){
			// cout << num2DNA(temp.kmer, kmer_length) << endl;
				all_kmer_info.push_back(temp);
				target_loc.push_back(i);
			}
		}

		// cout << "done with finding kmers, found " << all_kmer_info.size() << endl;
		// collect all scores for each found kmer & read pair

		counter = 0;
		vector<int> lu_wang_scores;
		int target_length = target_string.size();
		for(auto &kmerInfo: all_kmer_info){
			for(auto &readInfo: kmerInfo.read_info){			
				string l = to_upper(convert(find_read(readInfo.read, max_read_length)));
				int begin_read = ((target_loc[counter] - readInfo.readLoc - read_buffer) < 0) ? 0 : (target_loc[counter] - readInfo.readLoc - read_buffer);
				int end_read = ((target_loc[counter] + l.size() + read_buffer - readInfo.readLoc) > target_length) ? target_length : (target_loc[counter] + l.size() + read_buffer - readInfo.readLoc);
				string k = to_upper(target_string.substr(begin_read, end_read - begin_read));
				lu_wang_scores.push_back(lu_wang(k, l,  PERCENT, DEL, INS, MISMATCH, MATCH)); //k is vertical col. l is horizontal row
				if(lu_wang_scores[lu_wang_scores.size() - 1] >= PERCENT * l.size()){
					cout << "Comparison #" << print_counter << " kmerInfo.kmer " << num2DNA(kmerInfo.kmer,kmer_length) << " readInfo.read " << readInfo.read << " readInfo.readLoc " << readInfo.readLoc << " target_number " << target_number << " target_loc " << target_loc[counter] << '\n';
					cout << "col: " << k << '\n' << "row: " << l << '\n';
					cout << "LWScore: " << lu_wang_scores[lu_wang_scores.size() - 1] << '\n';
					++print_counter;
				}
			}
			++counter;
		}
		// figure out how to process the scores and give back a numerical answer!

		++target_number;
	}
	double total_align_time = (clock() - start) / (double)CLOCKS_PER_SEC;
	cout << "binary file generation time: " << total_bin_time << '\n';
	cout << "hash table indexing time: " << total_index_time << '\n';
	cout << "alignment time: " << total_align_time << '\n';
}

unsigned char* process_read(unsigned char* read, int read_loc, int kmer_size, int read_buffer){
	int num_nucleotides = (int)read[0];

	int begin_read = ((read_loc - read_buffer) < 0) ? 0 : (read_loc - read_buffer);
	int end_read = ((read_loc + kmer_size + read_buffer) > num_nucleotides) ? num_nucleotides : (read_loc + kmer_size + read_buffer);
	unsigned char *output_read = new unsigned char [end_read - begin_read + 1];
	int output_read_index = 0;
	int shift = 0;
	output_read[0] = end_read - begin_read;
	for(int i = begin_read; i < end_read; ++i){
		// fine tuning for char_shift
		int char_index = i / 4 + 1;
		int char_shift = 3 - i % 4;
        
		int DNA_val = ((read[char_index] >> (char_shift * 2)) & 0x3);
		int output_shift = (6 - ((shift) % 8));
		if(output_shift == 6){
			++output_read_index;
			output_read[output_read_index] = 0;
		}

		output_read[output_read_index] +=  DNA_val << output_shift;
		shift += 2;
	}

	return output_read;
}

void printHelp(){
	cout << "This program takes in the read filename, which can be either fastq format or pure reads format";
	cout << " and searches these reads for the string(s) specified in <target_read_filename>.\n";
	cout << "During this process, the reads are placed into a binary file specified by binary_read_filename tag, which also requires the user enter the maximum read size.\n";
	cout << " The user can also specify the kmer size and prefix size, otherwise they are set to a default 20 and 5 respectively\n";
	cout << "required arguments: \n";
	cout << "-t --target_string_filename <target_string_filename> enter the filename of the input target string\n";
	cout << "-r --read_filename <read_filename>: enter the filename of the input reads\n";
	cout << "-o --output_binary_filename <binary_filename>: enter the filename of the binary file output\n";
    cout << "OR\n";
    cout << "-i --input_binary_filename <binary_filename>: enter the filename of the binary file input so that the program does not need to re-write a binary file\n";
	cout << "-m --max_read_size <max read size>: enter the maximum number of reads in your read file\n";
	cout << "optional arguments: \n";
	cout << "-k --kmer_size <kmer size>: enter a number that specifies the kmer size (default: 20)\n";
	cout << "-p --prefix_size <prefix size>: enter a number that specifies the prefix size (default: 5)\n";
	cout << "-b --buffer <buffer size>: enter a number that specifies how large the buffer (on each side) should be when comparing the target string and read (default: 10)\n";
	cout << "-h --help displays help message\n";
	cout << "-d --debug displays debugging output\n";
}

bool checkNumeric(char *c_str){
	bool numeric = true;
	string temp = c_str;
	for(auto &x: temp){
		if(!isdigit(x)){
			numeric = false;
			break;
		}
	}
	return numeric;
}

string to_upper(string x){
    for(auto &i: x)
        i = toupper(i);
    return x;
}