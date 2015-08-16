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
 * nanopore_signal.hpp
 *
 *  Created on: Jun 12, 2015
 *      Author: nicola
 *
 *  brief: this class represents a series of nanopore events, each characterized by a duration, a mean, and a standard deviation. The class represents also the
 *  pore model (i.e. set of gaussian distributions) associated with the signal.
 *
 *	sample usage:
 *
 *	string filename(argv[1]);
 *  fast5::File f_p(filename);
 *
 *  nanopore_signal ns_t(f_p,TEMPLATE);		//template signal
 *  nanopore_signal ns_c(f_p,COMPLEMENT);	//complement signal (for 2D reads)
 *
 *
 */

#ifndef NANOPORE_SIGNAL_HPP_
#define NANOPORE_SIGNAL_HPP_

#include <fast5.hpp>
#include <tuple>
#include <vector>
#include <definitions.hpp>
#include <kmer.hpp>

using namespace std;

namespace ontrc{

class nanopore_signal{

public:

	/*
	 * constructor
	 *
	 * \param f_p the fast5 file
	 * \param s is this the TEMPLATE or COMPLEMENT strand?
	 * \param auto_correct automatically correct normal distributions and events using ONT correction parameters?
	 *
	 */
	nanopore_signal(fast5::File &f_p, strand s = TEMPLATE, bool auto_correct = true){

		str = s;

		int str_int = (s==TEMPLATE?0:1);

		auto params = f_p.get_model_parameters(str_int);

    	/*cout << "\ndrift = " << params.drift << endl;
    	cout << "scale = " << params.scale << endl;
    	cout << "shift = " << params.shift << endl;
    	cout << "var = " << params.var << endl;
    	cout << "scale_sd = " << params.scale_sd << endl;
    	cout << "var_sd = " << params.var_sd << endl;*/

		//If samplig rate is present in the fast5 file, use it. Otherwise, use the default 3000Hz.
		sampling_rate = 3000;
		if(f_p.have_sampling_rate())
			sampling_rate = f_p.get_sampling_rate();

		if( f_p.have_events(str_int) ){

			auto v = f_p.get_events(str_int);

            for (const auto& e : v)
            	events.push_back({e.mean,e.length,e.stdv,e.start,ulint(e.length*sampling_rate)});

		}

		if(f_p.have_model(str_int)){

			//4^k normal distributions
			model = vector<normal_d>(number_of_kmers());

            auto v = f_p.get_model(str_int);

            //get normal distributions

            for (const auto& e : v)
            	model[ kmer(string(e.kmer)) ] = {e.level_mean, e.level_stdv, e.sd_mean, e.sd_stdv};


            //automatically correct distributions and events (can be done only if f_p.have_model(str_int) )
            if(auto_correct){

                auto params = f_p.get_model_parameters(str_int);

                /*cout << "params.drift = " << params.drift << endl;
            	cout << "params.scale = " << params.scale << endl;
            	cout << "params.scale_sd = " << params.scale_sd << endl;
            	cout << "params.shift = " << params.shift << endl;
            	cout << "params.var = " << params.var << endl;
            	cout << "params.var_sd = " << params.var_sd << endl;*/

            	//correct the events
                //we disabled this correction because it messes up
                //everything... without drift correction, the caller
                //works just fine
            	/*if( f_p.have_events(str_int) ){

            		for(ulint i=0;i<events.size();++i){

            			double start = events[i].start_time;

            			events[i].mean -= (start*params.drift);

            		}

            	}*/

            	//correct the normal distributions

            	for(kmer km(k);km<number_of_kmers();km++){

            		//These transformations are provided by ONT

            		//mean and stdev of the signal mean
            		model[km].mean = model[km].mean*params.scale + params.shift;
            		model[km].stdev = model[km].stdev*params.var;

            		//mean and stdev of the signal stdev
                    model[km].sd_mean = model[km].sd_mean * params.scale_sd;
                    model[km].sd_stdev = model[km].sd_stdev * sqrt(pow(params.scale_sd, 3.0) / params.var_sd);

            	}

            }

		}

	}

	/*
	 * get event number i
	 */
	event operator[](ulint i){

		return get_event(i);

	}

	event get_event(ulint i){

		assert(i<events.size());
		return events[i];

	}

	bool has_events(){return events.size()>0;}

	bool has_model(){return model.size()>0;};

	/*
	 * number of events
	 */
	ulint size(){ return events.size(); }

	/*
	 * 4^k
	 */
	ulint number_of_kmers(){return ulint(1)<<(2*k);}

	uint8_t get_k_value(){return k;}

	/*
	 * get parameters of the normal distribution associated with kmer km
	 */
	normal_d get_distribution(kmer km){ return model[km]; }

	double stdev(kmer km){ return get_distribution(km).stdev; }
	double mean(kmer km){ return get_distribution(km).mean; }

	double sd_stdev(kmer km){ return get_distribution(km).sd_mean; }
	double sd_mean(kmer km){ return get_distribution(km).sd_stdev; }

	void replace_events(vector<event> events){
		this->events = events;
	}

	/*
	 * P(b|km) = log-likelihood of observing event b given kmer km
	 * accepts also a scaling parameter for the variance.
	 */
	double log_P(event e, kmer km, double scale_var=1){

		double var = scale_var*square(stdev(km)); //variance associated with the kmer
		double mu = mean(km); //mean associated with the kmer

		double e_var = square(e.stdv);
		double e_mean = e.mean;

		//use only mean
		//double logp = 0.5*log_1_2pi - (square( e_mean - mu ))/(2*var);

		//use mean and stdev
		double logp = log_1_2pi -
						(1/(2*var)) *
						( square( e_mean - mu ) +
						  square( e_var - var )/(2*var)
						);

		return logp;

	}

	strand get_strand(){return str;}

private:

	uint8_t k = 5;
	vector<event> events;
	strand str;

	//4^k normal distributions (one per kmer)
	vector<normal_d> model;

	//frequency (Hz) of sampling
	double sampling_rate = 0;

};

}

#endif /* NANOPORE_SIGNAL_HPP_ */
