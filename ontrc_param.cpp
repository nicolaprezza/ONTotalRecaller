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
 * ontrc_param.cpp
 *
 *  Created on: Jun 30, 2015
 *      Author: nicola
 *
 *  infer HMM parameters from a set of FAST5 files
 *
 */

#include <string>
#include <iostream>

#include <definitions.hpp>
#include <nanopore_signal.hpp>
#include <on_hmm.hpp>

using namespace ontrc;

typedef on_hmm::parameters_t hmm_parameters;

void help(){
	cout << "ontrc_param : compute HMM parameters from a set of FAST5 files" << endl << endl;
	cout << "Usage: ontrc_param <output_file> <FAST5 file 1>...<FAST5 file n>" << endl;
	cout <<	"   <output_file>     Output parameters file. If file already exists, then parameters are updated with the new data. Otherwise, a new file is created." << endl;
	cout <<	"   <FAST5 file 1>...<FAST5 file n>   input FAST5 files." << endl;
	exit(0);
}

/*
 * extract k value from fast5 file
 */
uint8_t get_k_value(string fast5_filename){

	fast5::File* f_p;
	f_p = new fast5::File(fast5_filename);

	assert(f_p->is_open());

	if(f_p->have_model(0)){

		auto v = f_p->get_model(0);
		uint8_t k_val = string(v[0].kmer).size();
		delete f_p;

		return k_val;

	}

	//Default: k = 5
	delete f_p;
	return 5;

}

void update_parameters(hmm_parameters& parameters, string fast5_filename){

	fast5::File* f_p;
	f_p = new fast5::File(fast5_filename);

	// Check that it opened successfully
	assert(f_p->is_open());

	//if base calls are not available, return.
	if(not f_p->have_basecalled_2D()){

		cout << " File does not have base calls. Skipping file." << endl;

		delete f_p;
		return;

	}

	//template strand
	if (f_p->have_model(0) and f_p->have_events(0)){

		if(f_p->have_basecalled_2D()){

			parameters.tot_template_events += f_p->get_events(0).size();
			parameters.tot_template_called_bases += f_p->basecalled_2D().size();

		}

		nanopore_signal ns_t(*f_p,TEMPLATE);

		std::vector< fast5::Event_Alignment_Entry > temp_events;//only events that have a template_index>=0

		{
			auto v = f_p->get_event_alignments();
			for(auto e:v)
				if(e.template_index>=0)
					temp_events.push_back(e);

		}

		cout << " processing " << temp_events.size() << " events on template strand ..." << flush;

		for (uint64_t i = 0;i<temp_events.size()-1;++i){

			auto km_now = kmer(string(temp_events[i].kmer));//current called kmer
			auto km_next = kmer(string(temp_events[i+1].kmer));//next called kmer

			auto chain = km_now.chain(km_next);// chain of kmers (length > 2 if there are skipped kmers)

			assert(chain.size()>0);

			parameters.temp[km_now].tot_e_cnt++;

			if(chain.size()==2){

				//E -> E transition of type 2 (next kmer)
				parameters.temp[km_now].ee2_cnt++;

			}else if(chain.size()==1){

				//E -> E transition of type 1 (we stay in the same kmer)
				parameters.temp[km_now].ee1_cnt++;

			}else{

				//chain.size() > 2: there are skipped kmers

				//don't need to increment ev_cnt (which does not exist), since all
				//counters must sum up to tot_e_cnt

				//for each skipped kmer except the last, increment V->V counte
				//Note: first and last kmers are emitting!
				for(ulint j=1;j<chain.size()-2;++j){

					parameters.temp[chain[j]].tot_v_cnt++;
					parameters.temp[chain[j]].vv_cnt++;

				}

				//for the last skipped kmer, increment V->E counter
				//(actually, just tot_v_cnt counter since ve_cnt is not
				//Explicitly stored)

				parameters.temp[chain[chain.size()-2]].tot_v_cnt++;

			}

		}

		//now compute probabilities for each kmer, plus events/called bases rate

		parameters.bases_per_template_event = (double)parameters.tot_template_called_bases/(double)parameters.tot_template_events;

		for(auto& cnt : parameters.temp){

			if(cnt.tot_e_cnt > 0){

				cnt.ee1 = (double)cnt.ee1_cnt/(double)cnt.tot_e_cnt;
				cnt.ee2 = (double)cnt.ee2_cnt/(double)cnt.tot_e_cnt;

				assert(cnt.ee1+cnt.ee2<=1);
				cnt.ev = 1 - (cnt.ee1+cnt.ee2);

			}

			if(cnt.tot_v_cnt > 0){

				cnt.vv = (double)cnt.vv_cnt/(double)cnt.tot_v_cnt;

				assert(cnt.vv <= 1);
				cnt.ve = 1 - cnt.vv;

			}

		}

		cout << " done." << endl;

	}else{

		cout << " File does not have template events or template model. Skipping template." << endl;

	}

	//complement strand
	if (f_p->have_model(1) and f_p->have_events(1)){

		if(f_p->have_basecalled_2D()){

			parameters.tot_complement_events += f_p->get_events(1).size();
			parameters.tot_complement_called_bases += f_p->basecalled_2D().size();

		}

		nanopore_signal ns_c(*f_p,COMPLEMENT);

		std::vector< fast5::Event_Alignment_Entry > comp_events;//only events that have a complement index>=0

		{

			auto v = f_p->get_event_alignments();

			for(ulint i = 0;i<v.size();++i){

				//read events backwards
				auto e = v[v.size()-i-1];

				if(e.complement_index>=0)
					comp_events.push_back(e);

			}

		}

		cout << " processing " << comp_events.size() << " events on complement strand ..." << flush;

		for (uint64_t i = 0;i<comp_events.size()-1;++i){

				auto km_now = kmer(string(comp_events[i].kmer)).rev_compl();//current called kmer
				auto km_next = kmer(string(comp_events[i+1].kmer)).rev_compl();//next called kmer

				auto chain = km_now.chain(km_next);// chain of kmers (length > 2 if there are skipped kmers)

				assert(chain.size()>0);

			parameters.comp[km_now].tot_e_cnt++;

			if(chain.size()==2){

				//E -> E transition of type 2 (next kmer)
				parameters.comp[km_now].ee2_cnt++;

			}else if(chain.size()==1){

				//E -> E transition of type 1 (we stay in the same kmer)
				parameters.comp[km_now].ee1_cnt++;

			}else{

				//chain.size() > 2: there are skipped kmers

				//don't need to increment ev_cnt (which does not exist), since all
				//counters must sum up to tot_e_cnt

				//for each skipped kmer except the last, increment V->V counte
				//Note: first and last kmers are emitting!
				for(ulint j=1;j<chain.size()-2;++j){

					parameters.comp[chain[j]].tot_v_cnt++;
					parameters.comp[chain[j]].vv_cnt++;

				}

				//for the last skipped kmer, increment V->E counter
				//(actually, just tot_v_cnt counter since ve_cnt is not
				//Explicitly stored)

				parameters.comp[chain[chain.size()-2]].tot_v_cnt++;

			}

		 }

		 //now compute probabilities for each kmer, plus events/called bases rate

		 parameters.bases_per_complement_event = (double)parameters.tot_complement_called_bases/(double)parameters.tot_complement_events;

		 for(auto& cnt : parameters.comp){

			if(cnt.tot_e_cnt > 0){

				cnt.ee1 = (double)cnt.ee1_cnt/(double)cnt.tot_e_cnt;
				cnt.ee2 = (double)cnt.ee2_cnt/(double)cnt.tot_e_cnt;

				assert(cnt.ee1+cnt.ee2<=1);
				cnt.ev = 1 - (cnt.ee1+cnt.ee2);

			}

			if(cnt.tot_v_cnt > 0){

				cnt.vv = (double)cnt.vv_cnt/(double)cnt.tot_v_cnt;

				assert(cnt.vv <= 1);
				cnt.ve = 1 - cnt.vv;

			}

		}

		cout << " done." << endl;

	}

	delete f_p;

}

int main(int argc, char** argv){

	uint64_t threshold = 50;

	if(argc<3) help();

	string param_filename(argv[1]);

	vector<string> fast5_filenames;
	for(uint i=2;i<argc;++i)
		fast5_filenames.push_back(argv[i]);

	bool param_exists=false;
	std::ifstream p_f(param_filename);
	param_exists = p_f.good();
	p_f.close();

	uint8_t k_val = 5;

	try{

		k_val = get_k_value(fast5_filenames[0]);

	}catch(const hdf5_tools::Exception& e) {}

	auto parameters = hmm_parameters(k_val);

	cout << "k is " << (ulint)parameters.k << endl;
	cout << "number of kmers : " << parameters.temp.size() << endl;

	//if parameters file exists, load it. Otherwise, all counters are initialized to 0
	if(param_exists){

		cout << "Parameters file exists. Updating frequencies. Note: existing parameters file will be overwritten." << endl;
		parameters.load_from_file(param_filename);

		std::remove(param_filename.c_str());//file will be overwritten

	}else{

		cout << "Parameters file does not exist. A new file will be created." << endl;

	}

	//loop through all fast5 files and update parameters
	for(auto in : fast5_filenames){

		try{

			cout << "Processing file " << in << endl;
			update_parameters(parameters,in);

		}catch(const hdf5_tools::Exception& e) {

			//if malformed file, just skip it.
			cout << " Error while reading file. File will be skipped." << endl;

		}

	}

	cout << "done." << endl;

	//for all k-mers that still have tot counters<threshold, use as probability the average of probabilities of kmers with enough observations (>threshold)
	parameters.normalize(threshold);

	cout << endl << parameters.bases_per_template_event << " average bases per template event " << endl;
	cout << parameters.bases_per_complement_event << " average bases per complement event " << endl;

	//for(auto& cnt : parameters.temp) cout << cnt.e << " " << cnt.b << endl;

	parameters.save_to_file(param_filename);

}
