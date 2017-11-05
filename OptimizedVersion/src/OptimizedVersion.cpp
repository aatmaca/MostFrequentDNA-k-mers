//============================================================================
// Name        : OptimizedVersion.cpp
// Author      : Abdullah Atmaca
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================

#include <ctime>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <deque>
#include <algorithm>
using namespace std;

// I divide input file into partitions and find top kMers for each partition as candidates for this partition.
// Then, I combine these candidates of each partition and calculate exact frequencies for them.
// Finally, I sort exact frequencies of candidates and get the results for top kMers.

// I assume that we have at least 8 GB Ram capacity.
// I want to use memory in 3 parts:
// One is for finding candidates, one is for combining candidates and the last part for remaining operations.

void updateMyMap(unordered_map<string, int> &umap, const string &line,
		const int &kmersize);
void updateCandidateMap(unordered_map<string, int> &umap, const string &line,
		const int &kmersize);

int main(int argc, char *argv[]) {
	string filename;
	int kmersize, topcount;
	if (argc == 4) {
		filename = argv[1];
		kmersize = atoi(argv[2]);
		topcount = atoi(argv[3]);
	} else if (argc == 7) {
		filename = argv[2];
		kmersize = atoi(argv[4]);
		topcount = atoi(argv[6]);
	} else {
//		filename = "ERR047698.filt.fastq";
//		filename = "ERR068396.filt.fastq";
//		filename = "ERR055763_1.filt.fastq";
//		kmersize = 25;
//		topcount = 50;
		cout << "Argument list is not valid. I am waiting 4 or 7 arguments." << endl;
		return -1;
	}

	if (kmersize < 1) {
		cout << "kmersize is not valid." << endl;
		return -1;
	}
	if (topcount < 1) {
		cout << "topcount is not valid." << endl;
		return -1;
	}

	ifstream preReader(filename);

	if (preReader.good() == false) {
		cout << filename << " not found." << endl;
		return -1;
	}

	cout << "Processing file..." << endl;

	string pre;
	getline(preReader, pre);
	getline(preReader, pre);
	unsigned int kMerPerLine = pre.length() - kmersize + 1;

	cout << "Number of kmer per line...:" << kMerPerLine << endl;

	// we initialize total_number_of_lines as 2 because we already read first two lines.
	unsigned int total_number_of_lines = 2;
	while (getline(preReader, pre)) {
		total_number_of_lines++;
	}
	cout << "Total number of lines in input file:" << total_number_of_lines << endl;

	unsigned int total_number_of_DNA_lines = total_number_of_lines/4;
	cout << "Total number of lines containing DNA sequence:" << total_number_of_DNA_lines << endl;
	preReader.close();

	// Value given below is an experimental value which have produced best results in my tests.
	// This value can be updated according to inputs and memory capacity,
	// but since inputs and memory capacity are not given to us, I try to assign most suitable value.
	unsigned int proposedMapSize = 8000000; // 8 Million elements per map.

	unsigned int number_of_lines_per_map = proposedMapSize / kMerPerLine;
	unsigned int no_of_file_partitions = (total_number_of_DNA_lines
			/ number_of_lines_per_map) + 1;
	// I will use quarter of proposed map size because of memory limitations.
	unsigned int no_of_candidates_per_partition = (proposedMapSize/4)
			/ no_of_file_partitions;

	// will be used for calculation overall completion percentage.
	unsigned int one_percent_of_total_number_of_lines = total_number_of_DNA_lines
			/ 100;

	cout << "Maximum number of lines per map: " << number_of_lines_per_map
			<< endl;
	cout << "Number of file partitions: " << no_of_file_partitions << endl;
	cout << "Maximum number of candidates per partition: "
			<< no_of_candidates_per_partition << endl;

	// top kmer candidates of each partition will be combined in this map
	unordered_map<string, int> top_kmer_candidates;

	// umap will be used for finding top kmers of each partition.
	unordered_map<string, int> umap;
	ifstream reader(filename);
	string line;
	time_t now = time(0);
	unsigned int lineCounter = 0;
	for (unsigned int i = 0; i < no_of_file_partitions; ++i) {
		cout << "File Partition Number: " << i + 1 << "/"
				<< no_of_file_partitions << ", Elapsed time: "
				<< (time(0) - now) / 60
				<< " min, Overall completion percentage: "
				<< lineCounter / one_percent_of_total_number_of_lines
				<< "%" << endl;
		while (!reader.eof()) {
			getline(reader, line);
			getline(reader, line);

			updateMyMap(umap, line, kmersize);

			getline(reader, line);
			getline(reader, line);
			lineCounter++;

			if ((umap.size() > proposedMapSize) || reader.eof()) {

				deque<pair<int, string>> items;

				for (auto& it : umap) {
					items.push_back( { it.second, it.first });
				}

				umap.clear();

				if (items.size() > no_of_candidates_per_partition) {
					// extract the highest elements
					nth_element(items.begin(),     // beginning of range
							items.end() - no_of_candidates_per_partition, // element that should be sorted correctly
							items.end());      // end of range
				}

				unsigned int internalCounter = 0;
				for (auto it = items.rbegin(); it != items.rend(); ++it) {
					if (internalCounter < no_of_candidates_per_partition) {
						// we need only kmer info at this point. We will calculate exact frequencies later.
						top_kmer_candidates.insert( { (*it).second, 0 });
					} else {
						break;
					}
					internalCounter++;
				}

				items.clear();
				break;
			}
		}

	}
	reader.close();

	/**
	 * This part is for calculation of exact frequencies.
	 */
	ifstream finalReader(filename);
	cout << "Calculating exact frequencies... " << endl;
	string str;
	while (getline(finalReader, str)) {
		getline(finalReader, str);

		updateCandidateMap(top_kmer_candidates, str, kmersize);

		getline(finalReader, str);
		getline(finalReader, str);
	}

	finalReader.close();

	cout << "Sorting according to exact frequencies... " << endl;
	deque<pair<int, string>> resultDeque;

	for (auto& it : top_kmer_candidates) {
		resultDeque.push_back( { it.second, it.first });
	}
	top_kmer_candidates.clear();

	sort(resultDeque.begin(), resultDeque.end());

	cout << endl;
	cout << topcount << " most frequent k-mers: " << endl;

	//(including all equal frequent kmers for last frequency value)
	int minValueForTopKmers = resultDeque[resultDeque.size() - topcount].first;
	for (auto it = resultDeque.rbegin(); it != resultDeque.rend(); ++it) {
		if ((*it).first >= minValueForTopKmers) {
			cout << (*it).second << " :: " << (*it).first << endl;
		} else {
			break;
		}
	}
	resultDeque.clear();

	return 0;
}

/**
 * This method will be used for insert and update operations
 */
void updateMyMap(unordered_map<string, int> &umap, const string &line,
		const int &kmersize) {

	int kMerCounter = line.length() - kmersize + 1;

	// Declare an iterator to unordered_map
	unordered_map<string, int>::iterator it;

	// Declare a variable for kMer
	string kMer;
	for (int i = 0; i < kMerCounter; i++) {
		kMer = line.substr(i, kmersize);

		it = umap.find(kMer);

		if (it == umap.end()) {
			// Element Not Found
			umap.insert( { kMer, 1 });
		} else {
			// Element Found
			it->second = it->second + 1;
		}
	}
}

/**
 * This method will be used for calculation of exact frequencies and only have update operation.
 */
void updateCandidateMap(unordered_map<string, int> &umap, const string &line,
		const int &kmersize) {

	int kMerCounter = line.length() - kmersize + 1;

	// Declare an iterator to unordered_map
	unordered_map<string, int>::iterator it;

	// Declare a variable for kMer
	string kMer;
	for (int i = 0; i < kMerCounter; i++) {
		kMer = line.substr(i, kmersize);

		it = umap.find(kMer);

		if (it == umap.end()) {
			// Element Not Found
			// Do nothing...
		} else {
			// Element Found
			it->second = it->second + 1;
		}
	}
}
