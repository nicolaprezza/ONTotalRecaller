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
 *
 *  Created on: Jun 30, 2015
 *      Author: nicola
 *
 *  given FM index, fast5 file, and parameters file, perform base-call of the fast5 signal
 *
 */

#include <string>
#include <iostream>

#include <definitions.hpp>
#include <nanopore_signal.hpp>
#include <total_recaller.hpp>

using namespace ontrc;

void help(){
	cout << "ontrc_call : perform base calling of a ONT signal, using the reference as bayesian prior." << endl << endl;
	cout << "Usage: ontrc_call [options]" << endl;
	cout <<	"   --idx <fmi.sdsl.fmi>      FM index file created with ontrc_build. [Required]" << endl;
	cout <<	"   --param <parameters>      File containing parameters of the HMM, created with ontrc_param. [Required]" << endl;
	cout <<	"   --input <in.fast5>        FAST5 file containing the ONT signal [Required]" << endl;
	cout <<	"   --output <out.fa>         Output file. Sequence will be written here in fasta format. [Required]" << endl;
	cout <<	"   --w <value>               Reference weight in the range [0,1]. Default: 1." << endl;
	cout <<	"   --start <value>           Call events starting from number <value>. Default: 0." << endl;
	cout <<	"   --end <value>             Last called event is number <value>. If 0, last event is the last in the original sequence. Default: 0." << endl;
	cout <<	"   --verbose                 Verbose output. Default: false." << endl;
	exit(0);
}

string fmi_file;
string param_file;
string fast5_file;
string fasta_file;
double W;
ulint start_;
ulint end_;
bool verbose;

void parse_parameters(char** argv, int argc, int &idx){

	if(idx >= argc) help();

	if(string(argv[idx]).compare("--idx")==0){

		idx++;
		if(idx >= argc) help();

		fmi_file = string(argv[idx]);

		idx++;
		return;

	}

	if(string(argv[idx]).compare("--param")==0){

		idx++;
		if(idx >= argc) help();

		param_file = string(argv[idx]);
		idx++;
		return;

	}

	if(string(argv[idx]).compare("--input")==0){

		idx++;
		if(idx >= argc) help();

		fast5_file = string(argv[idx]);
		idx++;
		return;

	}

	if(string(argv[idx]).compare("--output")==0){

		idx++;
		if(idx >= argc) help();

		fasta_file = string(argv[idx]);
		idx++;
		return;

	}

	if(string(argv[idx]).compare("--w")==0){

		idx++;
		if(idx >= argc) help();

		W = atof(argv[idx]);
		idx++;
		return;

	}

	if(string(argv[idx]).compare("--start")==0){

		idx++;
		if(idx >= argc) help();

		start_ = atoi(argv[idx]);
		idx++;
		return;

	}

	if(string(argv[idx]).compare("--end")==0){

		idx++;
		if(idx >= argc) help();

		end_ = atoi(argv[idx]);
		idx++;
		return;

	}

	if(string(argv[idx]).compare("--verbose")==0){

		idx++;
		verbose=true;
		return;

	}

	cout << "Unrecognized '" << argv[idx] << "' option." << endl;
	help();

}

void check_param(){

	if( fmi_file.compare("")==0 or
		param_file.compare("")==0 or
		fast5_file.compare("")==0 or
		fasta_file.compare("")==0 )
		help();

}

int main(int argc, char** argv){

	// DEFAULTS
	W = 1;
	start_=0;
	end_=0;
	verbose=false;
	//END DEFAULTS

	//PARSE PARAMETERS
	int idx = 1;
	while(idx < argc) parse_parameters(argv, argc, idx);
	//END PARSE PARAMETERS

	//CHECK MANDATORY PARAMETERS
	check_param();
	//END CHECK MANDATORY PARAMETERS

	cout << "Loading FM index form file " << fmi_file << " ... " << flush;
	fmi_trie corpus;

	//Detect if suffix .sdsl.fmi is present. If it is, cut it.
	ulint p = fmi_file.find(corpus.get_extension());
	if(p != string::npos) fmi_file = fmi_file.substr(0,p);

	corpus.load_from_disk(fmi_file);
	cout << "done." << endl;

	//load HMM parameters
	cout << "Loading HMM parameters from file " << param_file << " ... " << flush;
	on_total_recaller::parameters_t param;
	param.load_from_file(param_file);
	cout << "done." << endl;

	//build ONTRC object
	on_total_recaller trc(&corpus,param);

	//load ON signal
	cout << "Loading FAST5 file from " << fast5_file << " ... " << flush;
    fast5::File f_p(fast5_file);
    nanopore_signal ns_t(f_p,TEMPLATE);
    nanopore_signal ns_c(f_p,COMPLEMENT);

    end_ = (end_==0 or end_ >= ns_t.size()?ns_t.size()-1:end_);
    start_ = (start_ >= end_ ? 0: start_);

    vector<ontrc::event> subsample;//TODO test
    for(ulint i=start_;i<=end_;++i)
    	subsample.push_back( ns_t.get_event(i) );
    ns_t.replace_events(subsample);


	cout << "done." << endl;

	cout << "W is " << W << endl;

    //call sequence
	cout << "Calling sequence (" << end_-start_ << " ONT events) ... " << flush;

	auto calls = trc.call(ns_t,W,verbose);

	cout << "done. Found a match of length " << calls[0].size() << " on the reference. " << endl;

	ofstream ofs(fasta_file);

	cout << "Saving sequences ..."<<flush;

	if(calls[0].size()>0){
		ofs << ">ontrc_reference_match_with_snps\n";
		ofs << calls[0] << "\n";
	}

	ofs << ">ontrc_reference_unbiased\n";
	ofs << calls[1] << "\n";

	cout << " done." << endl;

	ofs.close();


	cout << "Called sequence has been saved in " << fasta_file << endl;

}
