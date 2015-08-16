/*
 *  This file is part of ONTotalRecaller.
 *  Copyright (c) by
 *  Nicola Prezza 		<nicolapr@gmail.com>
 *  Bud Mishra 			<mishra@nyu.edu>
 *  Giuseppe Narzisi	<gnarzisi@nygenome.org>
 *
 *  ONTotalRecaller is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.

 *  ONTotalRecaller is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details (<http://www.gnu.org/licenses/>).
 */

/*
 * definitions.hpp
 *
 *  Created on: Jun 12, 2015
 *      Author: nicola
 */

#ifndef DEFINITIONS_HPP_
#define DEFINITIONS_HPP_

#include <string>
#include <sdsl/suffix_arrays.hpp>
#include <fast5.hpp>
#include <math.h>

using namespace std;
using namespace sdsl;

namespace ontrc{

typedef uint64_t ulint;
typedef int64_t lint;

typedef unsigned char uchar;

//as FM index we choose a wavelet tree over a succinct bitvector. Since we are
//dealing with DNA, it wouldn't make much sense to use a more compressed bitvector,
//e.g. a rrr_vector: it would just be slower.
typedef csa_wt<wt_huff<bit_vector>, 512, 1024> fm_index_t;

enum strand {TEMPLATE, COMPLEMENT};

double neg_inf = -std::numeric_limits<double>::infinity();
double inf = std::numeric_limits<double>::infinity();

//we use this value to label an undefined score (i.e. score to be computed in the DP matrix). 1 can never be the value of a log-likelihood
double UNDEFINED = 1;

static constexpr double log_sigma = log(4.0);	//log of alphabet cardinality
static constexpr double log_1_2pi = -log(2*3.141592653589793); // log(1/(2pi))

bool is_defined(double x){return x!=UNDEFINED;}

enum base {	A,C,G,T, base_end };

double max(double x, double y, double z){

	double maxxy = (x>y?x:y);
	return (maxxy>z?maxxy:z);

}

struct event{

	double mean;
	double length;
	double stdv;
	double start_time;
	ulint number_of_samples;

};

struct normal_d{

	double mean;
	double stdev;

	double sd_mean;
	double sd_stdev;

};

base complement(base b){

	if(b==A) return T;
	if(b==C) return G;
	if(b==G) return C;
	return A;

}

char complement_char(char b){

	if(b=='A') return 'T';
	if(b=='C') return 'G';
	if(b=='G') return 'C';
	if(b=='T') return 'A';
	return b;

}

base char_to_base(uchar c){

	if(c=='A') return A;
	if(c=='C') return C;
	if(c=='G') return G;
	return T;

}

uchar base_to_char(base b){

	if(b==A) return 'A';
	if(b==C) return 'C';
	if(b==G) return 'G';
	return 'T';

}



string reverse_complement(string seq){

	string rc;

	for(ulint i=0;i<seq.size();++i)
		rc += complement_char(seq[seq.size()-i-1]);

	return rc;

}

string reverse(string seq){

	string r;

	for(ulint i=0;i<seq.size();++i)
		r += seq[seq.size()-i-1];

	return r;

}

bool contains_separator(string s){

	for(auto c:s)
		if(c=='$')
			return true;

	return false;

}

double square(double x){return x*x;}

struct counters{

	//number of outgoing transitions of each kind
	uint64_t vv_cnt=0;
	uint64_t ee1_cnt=0;
	uint64_t ee2_cnt=0;

	//number of total outgoing transitions exiting states e and v
	uint64_t tot_e_cnt=0;
	uint64_t tot_v_cnt=0;

	//actual probabilities
	double ve=0;
	double vv=0;
	double ev=0;
	double ee1=0;
	double ee2=0;

};

string format(char c){

	switch(c){

		case 'a' : case 'A' : return "A"; break;
		case 'c' : case 'C' : return "C"; break;
		case 'g' : case 'G' : return "G"; break;
		case 't' : case 'T' : return "T"; break;
		case 'n' : case 'N' : return "N"; break;

		default: return "";break;

	}

	return "";

}

/*
 * read a fasta file and return the DNA string
 */
string read_fasta_file(string path){

	std::ifstream file(path);

	string sequence = string();

	ulint contig_nr=0;

	while (file.good()){

		string line = string();
		getline(file, line);

		if(line[0] != '>'){

			for(auto c:line) sequence = sequence.append(format(c));

		}else{

			if(contig_nr>0)
				sequence = sequence.append("$");//put separator between contigs

			contig_nr++;

		}

	}

	return sequence;

}

/*
 * replace Ns with random bases
 */
void remove_Ns(string &s){

	srand(time(NULL));

	for(ulint i=0;i<s.size();++i){

		if(s[i] == 'N'){

			uint8_t rb = rand()%4;

			switch(rb){

				case 0: s[i] = 'A'; break;
				case 1: s[i] = 'C'; break;
				case 2: s[i] = 'G'; break;
				default: s[i] = 'T';

			}

		}

	}

}

}

#endif /* DEFINITIONS_HPP_ */
