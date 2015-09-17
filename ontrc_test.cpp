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
 * ontrc_test.cpp
 *
 *  Created on: Sep 17, 2015
 *      Author: nicola
 */

#include <string>
#include <iostream>

#include <definitions.hpp>
#include <fmi_trie.hpp>

#include <nanopore_signal.hpp>
#include <on_hmm.hpp>
#include <total_recaller.hpp>

using namespace ontrc;

void help(){
	cout << "ontrc_test : load genome, mutate it with SNPs and index it. Then, align the reads on the mutated genome" << endl;
	cout << "             and call variations. Computes also precision of calls." << endl << endl;
	cout << "Usage: ontrc_test [options] <param_file> <fasta> <fast5 files>" << endl;
	cout <<	"   <param_file>     HMM parameters file" << endl;
	cout <<	"   <fasta>          fasta file" << endl;
	cout <<	"   <fast5 files>    regexp matching all input fast5 files to be aligned." << endl;
	cout <<	"   -r <arg>         Mutation rate in [0,1]. SNPs are inserted uniformly." << endl;
	cout <<	"   -p <arg>         align at most <arg> events (template strand) for each fast5 file" << endl;
	exit(0);
}

string param_file;
string fasta_file;
vector<string> fast5_files;
double mutation_rate=0;
ulint prefix_temp=0;

void parse_parameters(char** argv, int argc, int &idx){

	if(idx >= argc) help();

	if(string(argv[idx]).compare("-r")==0){

		idx++;
		if(idx >= argc) help();

		mutation_rate = atof(argv[idx]);

		cout << "Mutation rate is " << mutation_rate << endl;

		idx++;
		return;

	}

	if(string(argv[idx]).compare("-p")==0){

		idx++;
		if(idx >= argc) help();

		prefix_temp = atoi(argv[idx]);

		cout << "Aligning first " << prefix_temp << " events for each fast5 file" << endl;

		idx++;
		return;

	}

	//if no -r or -p encountered, then this is the param file
	param_file = string(argv[idx]);
	idx++;
	if(idx >= argc) help();

	//fasta file
	fasta_file = string(argv[idx]);
	idx++;
	if(idx >= argc) help();

	//now fast5 files
	while(idx<argc) fast5_files.push_back(argv[idx++]);

}

void mutate(string & s){

	srand(time(NULL));

	for(ulint i = 0;i<s.size();++i){

		if((double) rand() / RAND_MAX < mutation_rate){

			cout << "mutation in pos " << i << " from " << s[i] << " to " << flush;

			s[i] = base_to_char(base((char_to_base((uchar)s[i]) + ((rand()%3) + 1))%4));

			cout << s[i] << endl;

		}

	}

}

int main(int argc, char** argv){

	if(argc < 4) help();

	int idx = 1;
	while(idx < argc) parse_parameters(argv, argc, idx);

	cout << "Loading reference from file " << fasta_file << " ... " << flush;
	string forward_strand = read_fasta_file(fasta_file);
	cout << "done." <<endl;

	//copy original reference
	string original_reference = string(forward_strand);

	//mutate fwd strand

	mutate(forward_strand);

	//build reference text: forward strand concatenated to its reverse complement, plus a separator in the middle
	string reference = forward_strand + '$' + ontrc::reverse_complement(forward_strand);
	//string reference = forward_strand + '$';

	ontrc::remove_Ns(reference);

	cout << "Building FM index ... " << flush;
	fmi_trie fmi(reference);
	cout << "done."  << endl;

	//read parameters and build TRC object
	on_hmm::parameters_t param;
	param.load_from_file(param_file);

	on_total_recaller trc(&fmi,param);

	for(auto f5 : fast5_files){

		cout << "processing " << f5 << " ..." << endl;

		fast5::File f_p(f5);

		nanopore_signal ns_t(f_p,TEMPLATE);
		nanopore_signal ns_c(f_p,COMPLEMENT);

		//cout << ns_t.size() << " template events and " << ns_c.size() << " complement events. " << endl;

		ulint prefix_comp = (double) prefix_temp * ((double)ns_c.size()/ns_t.size());

		if(prefix_temp==0){

			prefix_temp = ns_t.size();
			prefix_comp = ns_c.size();

		}

		//cout << "prefix comp is " << prefix_comp << endl;

		//cut events

		//template
		vector<ontrc::event> subsample1;
		for(ulint i=0;i<prefix_temp;++i)
			subsample1.push_back( ns_t.get_event(i) );
		ns_t.replace_events(subsample1);

		//complement
		vector<ontrc::event> subsample2;
		for(ulint i=0;i<prefix_comp;++i)
			subsample2.push_back( ns_c.get_event(i) );
		ns_c.replace_events(subsample2);

		//align signal on reference

		auto calls = trc.call(ns_t,1, true,true);

		//if found alignment
		if(calls[0].first.size()>0){

			strand str = calls[0].second.s;






			on_hmm hmm_t(&ns_t,200,param);
			on_hmm hmm_c(&ns_c,200,param);

		}


	}


}

















