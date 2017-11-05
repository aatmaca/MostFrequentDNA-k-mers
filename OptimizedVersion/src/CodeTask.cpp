//============================================================================
// Name        : CodeTask.cpp
// Author      : Abdullah Atmaca
// Version     :
// Copyright   : Your copyright notice
// Description :
//============================================================================

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <deque>
#include <algorithm>
using namespace std;

void updateMyMap(unordered_map<string, int> &umap, const string &line,
		const int &kmersize);

int main(int argc, char *argv[]) {
	string filename;
	//filename = "ERR047698.filt.fastq";
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
		filename = "ERR055763_1.filt.fastq";
		filename = "ERR047698.filt.fastq";
		filename = "ERR068396.filt.fastq";
		kmersize = 45;
		topcount = 50;
		//return -1;
	}

	if (kmersize < 1) {
		cout << "kmersize is not valid." << endl;
		return -1;
	}
	if (topcount < 1) {
		cout << "topcount is not valid." << endl;
		return -1;
	}

	ifstream screen_reader(filename);

	if (screen_reader.good() == false) {
		cout << filename << " not found." << endl;
		return -1;
	}

	cout << "Processing file..." << endl;

	string line;
	unordered_map<string, int> umap;
	while (getline(screen_reader, line)) {
		getline(screen_reader, line);

		updateMyMap(umap, line, kmersize);

		getline(screen_reader, line);
		getline(screen_reader, line);
	}
	screen_reader.close();

	deque<pair<int, string>> items;

	for (auto& it : umap) {
		items.push_back( { it.second, it.first });
	}
	umap.clear();
	sort(items.begin(), items.end());

	int minValueForTopKmers = items[items.size() - topcount].first;
	for (auto it = items.rbegin(); it != items.rend(); ++it) {
		if ((*it).first >= minValueForTopKmers) {
			cout << (*it).second << " :: " << (*it).first << endl;
		} else {
			break;
		}
	}

	return 0;
}

void updateMyMap(unordered_map<string, int> &umap, const string &line,
		const int &kmersize) {

	int lineMax = line.length() - kmersize + 1;

	// Declare an iterator to unordered_map
	unordered_map<string, int>::iterator it;

	// Declare a variable for kMer
	string kMer;
	for (int i = 0; i < lineMax; i++) {
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
