#ifndef LU_WANG_H
#define LU_WANG_H
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
using namespace std;

typedef vector<vector<int> > TwodVec;
int score_read(string row, unsigned char* col, double match_percent, int del_penalty, int ins_penalty, int mismatch_penalty, int match_score){
	// create a two dimensional vector that match with the number
	int row_length = row.size();
	int col_length = col[0];
	vector<vector<int> > T(row_length + 1, vector<int>(col_length + 1));
	for(unsigned int i = 0; i < T.size(); ++i){
		T[i][0] = 0;
	}
	for(unsigned int j = 0; j < T[0].size(); ++j){
		T[0][j] = 0;
	}
	int max_value = -1;//, max_i = -1, max_j = -1;
	//cout << "row_length: " << row_length << " col length " << col_length << endl;
	for(unsigned int i = 0; i < row_length; ++i){
		for(unsigned int j = 0; j < col_length; ++j){
			// double complement = 1 - match_percent;
			//if(j <= i + complement * col_length && j <= i - (row_length - match_percent * col_length)){
				// calculate row[i] and col[i]
                int cur_row = mapper::DNAmap.at(row[i]);//(row >> ((row_length - i - 1) * 2)) & 0x3;
				// char* col needs to be found
				int char_index = j / 4 + 1;
				int char_shift = 3 - j % 4;
				int cur_col = ((col[char_index] >> (char_shift * 2)) & 0x3);
				int match = (int)(T[i][j]) + (cur_row == cur_col ? match_score : mismatch_penalty);
				int deletion = T[i][j + 1] + del_penalty;
				int insertion = T[i + 1][j] + ins_penalty;
				int value = max(max(match, deletion), insertion);
				T[i + 1][j + 1] = value;
				if(max_value < value && j == col_length - 1){
					max_value = value;
				}
			//}
		}
	}
	//cout << endl;
	//for(unsigned int i = 0; i < T.size(); ++i){
	//	for(unsigned int j = 0; j < T[i].size(); ++j){
	//		cout << T[i][j] << "\t";
	//	}
	//	cout << endl;
	//}

	//cout << max_value << " " << max_i << " " << max_j << endl;
	return max_value;
}
// row string is 100% globally aligned
// col string is locally aligned based on the match_percent
// all settings can  be configured with the penalty/match scores
// change to work with longs!

int lu_wang(string row, string col, double match_percent, int del_penalty, int ins_penalty, int mismatch_penalty, int match_score){
	//  begin creating needleman-wunsch table
	TwodVec T (row.length() + 1, vector<int>(col.length() + 1, -999999));

	// write the first column as all 0
	for(unsigned int i = 0; i < T.size(); ++i){
//        T[i][0] = i * ins_penalty;
		T[i][0] = 0;
	}

	// write first row as 0, del_penalty, 2 * del_penalty, 3 * del_penalty....

	for(unsigned int j = 0; j < T[0].size(); ++j){
//        T[0][j] = j * del_penalty;
		T[0][j] = ins_penalty*j;
	}

	// calculate hybrid local/global alignment
	// local alignment for ROW, global alignent for COL
	for(unsigned int i = 0; i < row.length(); ++i){
		for(unsigned int j = 0; j < col.length(); ++j){
			// allow only
			double complement = 1 - match_percent;
			if (j <= i + complement * col.length() && j >= i - (row.size() - match_percent * col.length())) {
//                ++T[i + 1][j + 1];
				int match = (int)(T[i][j]) + (row[i] == col[j] ? match_score : mismatch_penalty);   //diag
				int deletion = T[i][j + 1] + del_penalty;                             //up
				int insertion = T[i + 1][j] + ins_penalty;                            //left
				int value = max(max(match, deletion), insertion);
				//            int value = max(match, insertion);
				//            value = max(value, deletion);
				T[i + 1][j + 1] = value;
//				if (max_value < value && j == col.length() - 1){
//					max_value = value;
//				}
			}
		}
	}
	//  end creating smith-waterman table
    int max_value;
    for(unsigned int i = 1; i < row.length() + 1; ++i){
        if(T[i][col.length()] > max_value || i == 1)
            max_value = T[i][col.length()];
    }
	// print table
//	for(unsigned int i = 0; i < T.size(); ++i){
//		for(unsigned int j = 0; j < T[i].size(); ++j){
//			cout << T[i][j];// << "\t";
//		}
//		cout << endl;
//	}

//    cout << LWscore << endl;
	// cout << "LWscore: " << max_value << endl;
	return max_value;
}
#endif /* LU_WANG ENDIF */
