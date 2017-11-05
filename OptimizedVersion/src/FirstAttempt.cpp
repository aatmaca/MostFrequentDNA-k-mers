//============================================================================
// Name        : aaa.cpp
// Author      : Atmaca
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#define __USE_LARGEFILE64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSETS_BITS 64

#if __WIN32__
#define stat64 _stat64
#define mystruct __stat64
#else
#define mystruct stat64
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <bitset>
#include <sstream>

#include <sys/stat.h>
using namespace std;


long GetFileSize(string filename);
unsigned int FileRead( istream & is, vector <char> & buff );
unsigned int CountLines( const vector <char> & buff, int sz );

int main(int argc, char * argv[]) {
	string filename = "ERR055763_1.filt.fastq";
//	cout << GetFileSize(filename) << endl;

	//cout << binToHex(kMerToBinary("CCGGG")) << endl;


	time_t now = time(0);
	    if ( argc == 1  ) {
	        cout << "lines\n";
	        ifstream ifs( filename );
	        int n = 0;
	        string s;
	        while( getline( ifs, s ) ) {
	            n++;
	        }
	        cout << n << endl;
	    }
	    else {
	        cout << "buffer\n";
	        const int SZ = 1024 * 1024;
	        std::vector <char> buff( SZ );
	        ifstream ifs( filename );
	        int n = 0;
	        while( int cc = FileRead( ifs, buff ) ) {
	            n += CountLines( buff, cc );
	        }
	        cout << n << endl;
	    }
	    cout << time(0) - now << endl;







	return 0;
}
unsigned int FileRead( istream & is, vector <char> & buff ) {
    is.read( &buff[0], buff.size() );
    return is.gcount();
}

unsigned int CountLines( const vector <char> & buff, int sz ) {
    int newlines = 0;
    const char * p = &buff[0];
    for ( int i = 0; i < sz; i++ ) {
        if ( p[i] == '\n' ) {
            newlines++;
        }
    }
    return newlines;
}
long GetFileSize(string filename) {
	struct mystruct buf[26];
	int result = stat64(filename.c_str(), buf);
	return result == 0 ? buf->st_size / (1024 * 1024) : -1;
}

string charToBinaryString(char &c) {
	switch (c) {
	case 'A':
		return "00";
	case 'T':
		return "01";
	case 'G':
		return "10";

	default: //if (c == 'C')
		return "11";
	}
}

string kMerToBinary(string kMer) {

	string binary_str = "";
	for (unsigned int i = 0; i < kMer.size(); i++) {
		binary_str += charToBinaryString(kMer.at(i));
	}

	return binary_str;
}

long binToHex(string binary_str) {
	bitset<4> set(binary_str);
	return set.to_ulong();
}

string binToHexString(string binary_str, const int kmersize) {

	bitset<4> set(binary_str);
	cout << hex << set.to_ulong() << endl;

	stringstream res;
	res << hex << set.to_ulong();
	return res.str();
}
