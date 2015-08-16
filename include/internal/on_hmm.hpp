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
 * on_hmm.hpp
 *
 *  Created on: Jun 24, 2015
 *      Author: nicola
 *
 *  This class implements the Hidden Markov model used to call Oxford Nanopore events. The HMM is
 *  implemented as a Dynamic Programming algorithm. The DP matrix is a vector of (complete) trees; each entry DP[node, i] represents the
 *  log-likelihood of being in node 'node' at time i, i.e. the log-likelihood that the DNA sequence seq(node) has generated the
 *  nanopore signal beta_1, ..., beta_i
 *
 *  The functions accessed by TotalRecaller (and which must be implemented by this class in order to be used by TRC) are:
 *
 *	- get_root() : root node of the tree
 *	- get_child(node_t n, base b) : children of node n that is reached descending the edge labeled b
 *	- log_likelihood(node_t* x, ulint i) : log-likelihood that the DNA sequence seq(x) has generated the
 *    nanopore signal beta_0, ..., beta_i (note: event indexes start from 0)
 *  - get_exhaustive_depth() : inclusive depth at which ALL tree nodes should be evaluated. Example: for ONT, this is k+1 because
 *      nodes below height k do not have a defined score and we want to take into account at least 1 signal emission in order to
 *      take a decision on which path to take.
 *
 *	The class must moreover define:
 *
 *	- the nodes type node_t
 *	- the signal type signal_t
 *	- a struct/class named precomputed_tables which must contain all information (e.g. pre-computed arrays)
 *	  shared between multiple istances of this class. This is needed in order to avoid building multiple
 *	  times heavy tables
 *	- a static method get_precomputed_tables() that builds and returns all precomputed arrays
 *	- a type parameters_t, containing all parameters of the underlying model (in this case, a HMM)
 *
 *	Note that the memory management must be done OUTSIDE this class: methods get_root() and get_child(..) return pointers which must be freed by the user
 *
 */

#ifndef ON_HMM_HPP_
#define ON_HMM_HPP_

#include <definitions.hpp>
#include <nanopore_signal.hpp>
#include <kmer.hpp>

using namespace std;

namespace ontrc{

class on_hmm{

private:

	class node_internal{

		using node_ptr = node_internal*;

	public:

		//empty constructor
		//node_internal(){}

		/*
		 * Root constructor: builds the root node.
		 */
		node_internal(on_hmm * tree){

			this->tree = tree;

			//initialize new kmer
			//km = kmer(tree->get_k_value());

			depth_ = 0;

			//allocate memory for the 4 children
			children_ = vector<node_ptr>(4,NULL);

			//root is its own parent
			parent_ = this;

			b = base_end;

		}

		/*
		 * Child constructor: build descendent of a node
		 */
		node_internal(node_ptr x, base b){

			this->tree = x->tree;

			//left-shift parent's kmer and add base b to the right
			//km = x->get_kmer();
			//km += b;

			this->b = b;

			depth_ = x->depth() + 1;

			parent_ = x;

			children_ = vector<node_ptr>(4,NULL);

			//assign this node as b-th children of the parent node x
			x->children_[b] = this;

			//compute size of E_ and V_ arrays

			//only nodes after depth k have vectors V_ and E_
			if(depth()>=tree->get_k_value()){

				uint half_band = floor(tree->bandwidth/2);//half the size of the band

				//this is the sequence index in the DP matrix. i is the event index
				uint seq_idx = depth()-tree->get_k_value();

				//compute central diagonal position on the events axis
				uint event_diagonal = (uint)((double)seq_idx/tree->bases_per_event);

				//compute (inclusive) extremes of the band on the events axis

				band_left = (half_band>event_diagonal?0:event_diagonal-half_band);
				band_right = half_band+event_diagonal;

				//cout << half_band << " " << band_left << "  " << band_right << endl;

				uint vec_size = (band_right-band_left)+1;

				V_ = vector<double>(vec_size, UNDEFINED);
				E_ = vector<double>(vec_size, UNDEFINED);

			}

			assert(parent_->children_.size()==4);
			assert(depth()==x->depth()+1);

		}

		void erase_scores(){

			ulint vec_size = V_.size();

			//nothinf to do if E_ and V_ are empty (i.e. node at depth < k)
			if(vec_size>0){

				V_ = vector<double>(vec_size, UNDEFINED);
				E_ = vector<double>(vec_size, UNDEFINED);

			}

		}

		bool is_leaf(){

			return  (not has_child(A)) and
					(not has_child(C)) and
					(not has_child(G)) and
					(not has_child(T));

		}

		node_ptr get_parent(){

			assert(tree!=NULL);
			assert(parent_->children_.size()==4);
			assert(is_root() or depth()==parent_->depth()+1);

			return parent_;

		}

		/*
		 * debugging function
		 */
		void check_tree(){

			assert(children_.size()==4);

			assert(tree!=NULL);
			assert(parent_!=NULL);
			assert(is_root() or depth()==parent_->depth()+1);
			assert(parent_->children_.size()==4);

			if(not is_root())
				get_parent()->check_tree();

		}

		node_ptr get_child(base b){

			assert(tree!=NULL);
			assert(has_child(b));

			return children_[b];

		}

		ulint depth(){

			assert(tree!=NULL);
			return depth_;

		}

		bool is_root(){

			assert(tree!=NULL);
			return depth_==0;

		};

		base get_base(){return b;}
		void set_base(base b1){b = b1;}

		//these functions access arrays E_ and V_ and compute the band automatically
		//input: event number, >=0
		//output: a reference to the corresponding array element. Consequently, these functions
		//are l-values and can be assigned a value
		double& E(ulint i){

			assert(tree!=NULL);

			if(depth()<tree->get_k_value()) return neg_inf;

			//if outside band, return -inf
			if(i<band_left or i>band_right) return neg_inf;

			//else, return ref to value stored in array E_
			return E_[i-band_left];

		}

		double& V(ulint i){

			assert(tree!=NULL);

			if(depth()<tree->get_k_value()) return neg_inf;

			//if outside band, return -inf
			if(i<band_left or i>band_right) return neg_inf;

			//else, return ref to value stored in array V_
			return V_[i-band_left];

		}

		/*
		 * kmer is well-defined only if depth>=k!
		 */
		kmer get_kmer(){

			assert(tree!=NULL);
			uchar k = tree->get_k_value();

			if(depth_<k) return kmer(k);

			return get_kmer(k);

		}

		/*
		 * does this node have child b?
		 */
		bool has_child(base b){

			assert(tree!=NULL);
			return children_[b] != NULL;

		}

		ulint get_band_left(){
			return band_left;
		}

		ulint get_band_right(){
			return band_right;
		}

		//inclusive bounds of the banded DP (on events axis)
		uint band_left=0;
		uint band_right=0;

	private:

		/*
		 * kmer is well-defined only if depth>=k!
		 */
		kmer get_kmer(uchar i){

			assert(tree!=NULL);
			uchar k = tree->get_k_value();
			assert(depth_>= i);

			if(i==0) return kmer(k);

			//if i>0 this node cannot be the root!
			assert(not is_root());

			kmer km = get_parent()->get_kmer(i-1);
			km += get_base();

			return km;

		}

		//tree topology. Nodes are indexes in array 'nodes' of parent class
		node_ptr parent_ = NULL; //parent node
		vector<node_ptr> children_;

		//node's depth
		ulint depth_ = 0;

		//DP entries
		//ideally, here we would have m entries for each vector.
		//in practice, to speed-up computation here we will use
		//banded Dynamic Programming: outside the band, we assign
		//score = neg_inf
		vector<double> E_;
		vector<double> V_;

		//kmer associated with this node
		//kmer km;

		//base associated with this node
		base b;

		on_hmm * tree = NULL;

	};

	class hmm_parameters{

	public:

		hmm_parameters(){}

		hmm_parameters(uint8_t k){

			this->k = k;

			number_of_kmers = 1<<(2*k);//4^k

			temp = vector<counters>(number_of_kmers);
			comp = vector<counters>(number_of_kmers);

		}

		uint64_t save_to_file(string filename){

			/*
			 * file format: k value followed by a series of 64-bits counters and probabilities following the schema:
			 *
			 * k value
			 * e_1_template	b_1_template tot_1_template e b
			 * ...
			 * e_{4^k}_template	b_{4^k}_template b_{4^k}_template e b
			 * e_1_complement	b_1_complement tot_1_complement e b
			 * ...
			 * e_{4^k}_complement	b_{4^k}_complement b_{4^k}_complement e b
			 *
			 * (no tabs or newlines; file is binary). Counters reflect number of times a transition on the HMM has been taken in the fast5 file
			 *
			 */

			std::ofstream out(filename);

			uint64_t w_bytes = 0;

			out.write((char*)&tot_template_events,sizeof(tot_template_events));
			out.write((char*)&tot_complement_events,sizeof(tot_complement_events));
			out.write((char*)&tot_template_called_bases,sizeof(tot_template_called_bases));
			out.write((char*)&tot_complement_called_bases,sizeof(tot_complement_called_bases));
			out.write((char*)&bases_per_template_event,sizeof(bases_per_template_event));
			out.write((char*)&bases_per_complement_event,sizeof(bases_per_complement_event));

			w_bytes += sizeof(uint64_t)*4 + sizeof(double)*2;

			out.write((char*)&k,sizeof(k));
			w_bytes += sizeof(k);

			out.write((char*)temp.data(),sizeof(counters)*temp.size());
			out.write((char*)comp.data(),sizeof(counters)*comp.size());

			w_bytes += sizeof(counters)*(temp.size()+comp.size());

			out.close();

			return w_bytes;

		}

		void load_from_file(string filename){

			std::ifstream in(filename);

			in.read((char*)&tot_template_events,sizeof(tot_template_events));
			in.read((char*)&tot_complement_events,sizeof(tot_complement_events));
			in.read((char*)&tot_template_called_bases,sizeof(tot_template_called_bases));
			in.read((char*)&tot_complement_called_bases,sizeof(tot_complement_called_bases));
			in.read((char*)&bases_per_template_event,sizeof(bases_per_template_event));
			in.read((char*)&bases_per_complement_event,sizeof(bases_per_complement_event));

			in.read((char*)&k,sizeof(k));
			number_of_kmers = 1<<(2*k);//4^k

			uint64_t nr_of_kmers = 1<<(2*k);//4^k

			temp = vector<counters>(nr_of_kmers);
			comp = vector<counters>(nr_of_kmers);

			in.read((char*)temp.data(),sizeof(counters)*temp.size());
			in.read((char*)comp.data(),sizeof(counters)*comp.size());

			in.close();

		}

		void normalize(uint64_t threshold){

			//means
			double ev = 0;
			double ee1 = 0;
			double vv = 0;

			double ee2 = 0;
			double ve = 0;

			uint64_t good_kmers_e = 0;
			uint64_t good_kmers_v = 0;


			for(auto& cnt : temp){

				if(cnt.tot_e_cnt >= threshold){

					good_kmers_e++;
					ev += cnt.ev;
					ee1 += cnt.ee1;

				}

				if(cnt.tot_v_cnt >= threshold){

					good_kmers_v++;
					vv += cnt.vv;

				}

			}

			ev /= good_kmers_e;//this is the average e
			ee1 /= good_kmers_e;//this is the average b
			vv /= good_kmers_v;

			cout << "Kmers with at least " << threshold << " observations for e edges/template: " << good_kmers_e << endl;
			cout << "Kmers with at least " << threshold << " observations for v edges/template: " << good_kmers_v << endl;

			assert(ev+ee1<=1);
			assert(vv<=1);
			ee2 = 1 - (ev+ee1);
			ve = 1 - vv;

			cout << "average transition probabilities on template strand: " << endl;
			cout << " ev = " << ev << endl;
			cout << " ee1 = " << ee1 << endl;
			cout << " ee2 = " << ee2 << endl;
			cout << " ve = " << ve << endl;
			cout << " vv = " << vv << endl;

			for(auto& cnt : temp){

				if(cnt.tot_v_cnt < threshold or cnt.vv==0 or cnt.ve==0){
					cnt.ve = ve;
					cnt.vv = vv;
				}

				if(cnt.tot_e_cnt < threshold or cnt.ev==0 or cnt.ee1==0 or cnt.ee2==0){
					cnt.ev = ev;
					cnt.ee1 = ee1;
					cnt.ee2 = ee2;
				}


			}

			//now for complement

			ev = 0;
			ee1 = 0;
			vv = 0;

			ee2 = 0;
			ve = 0;

			good_kmers_e = 0;
			good_kmers_v = 0;

			for(auto& cnt : comp){

				if(cnt.tot_e_cnt >= threshold){

					good_kmers_e++;
					ev += cnt.ev;
					ee1 += cnt.ee1;

				}

				if(cnt.tot_v_cnt >= threshold){

					good_kmers_v++;
					vv += cnt.vv;

				}

			}

			ev /= good_kmers_e;//this is the average e
			ee1 /= good_kmers_e;//this is the average b
			vv /= good_kmers_v;

			cout << "Kmers with at least " << threshold << " observations for e edges/complement: " << good_kmers_e << endl;
			cout << "Kmers with at least " << threshold << " observations for v edges/complement: " << good_kmers_v << endl;

			assert(ev+ee1<=1);
			assert(vv<=1);
			ee2 = 1 - (ev+ee1);
			ve = 1 - vv;

			cout << "\naverage transition probabilities on complement strand: " << endl;
			cout << " ev = " << ev << endl;
			cout << " ee1 = " << ee1 << endl;
			cout << " ee2 = " << ee2 << endl;
			cout << " ve = " << ve << endl;
			cout << " vv = " << vv << endl;

			for(auto& cnt : comp){

				if(cnt.tot_v_cnt < threshold or cnt.vv==0 or cnt.ve==0){
					cnt.ve = ve;
					cnt.vv = vv;
				}

				if(cnt.tot_e_cnt < threshold or cnt.ev==0 or cnt.ee1==0 or cnt.ee2==0){
					cnt.ev = ev;
					cnt.ee1 = ee1;
					cnt.ee2 = ee2;
				}

			}

		}

		//counters to keep track of average number of bases per event seen in the fast5 files
		//this is needed in order to estimate length of sequence to be called
		uint64_t tot_template_events=0;
		uint64_t tot_complement_events=0;
		uint64_t tot_template_called_bases=0;
		uint64_t tot_complement_called_bases=0;
		double bases_per_template_event=0;
		double bases_per_complement_event=0;

		uint64_t number_of_kmers=0;
		uint64_t k=0;
		vector<counters> temp;//template strand
		vector<counters> comp;//complement strand

	};


public:

	//the nodes type
	typedef node_internal* node_t;
	typedef nanopore_signal signal_t;
	typedef hmm_parameters parameters_t;

	struct precomputed_tables{

		vector<double> log_1_e_;

	};

	/*
	 * constructor 1: precomputed log tables are computed inside the constructor
	 *
	 * \param nps nanopore signal
	 * \param n max depth: max DNA sequence length
	 * \param bandwidth: 	width of the DP matrix band. Is the maximum number of events
	 * 						that can be associated to a single sequence position in the DP matrix. If
	 * 						>= m, then matrix is full.
	 * \param parameters parameters of the HMM.
	 *
	 *
	 */
	on_hmm(nanopore_signal * nps, uint bandwidth, parameters_t& parameters){

		pt = get_precomputed_tables();
		build(nps, bandwidth, parameters, pt);

	}

	/*
	 * constructor 2: pass precomputed log tables from the extern
	 *
	 * \param nps nanopore signal
	 * \param n max depth: max DNA sequence length
	 * \param bandwidth: 	width of the DP matrix band. Is the maximum number of events
	 * 						that can be associated to a single sequence position in the DP matrix. If
	 * 						>= m, then matrix is full.
	 * \param parameters parameters of the HMM.
	 * \param precomputed_1e double array with precomputed values for funcion log(1+exp(x)). To build it, call
	 * 								pre_vec = get_precomputed_1e_vector(), and pass to this function pre_vec.data(),
	 * 								being careful not to destroy pre_vec before using this class! (this prevents
	 * 								the use of new and delete)
	 *
	 *
	 */
	on_hmm(nanopore_signal * nps, uint bandwidth, parameters_t& parameters, precomputed_tables& precomputed_t){

		build(nps, bandwidth, parameters, precomputed_t);

	}

	/*
	 * get identifier of this scoring strategy
	 */
	string get_ID(){
		return string("ONT");
	}

	/*
	 * get the root.
	 */
	node_t get_root() {

		return new node_internal(this);

	}

	/*
	 * create and return child node of n reached by descending the edge labeled with b
	 */
	node_t get_child(node_t x, base b){

		//if this child already exists, return it.
		if(x->has_child(b))
			return x->get_child(b);

		//otherwise, create new child.
		return new node_internal(x,b);

	}

	bool has_child(node_t x, base b){

		return x->has_child(b);

	}

	/*
	 * is node x the root?
	 */
	bool is_root(node_t x){

		return x->is_root();

	}

	/*
	 * By calling this function, the class computes recursively the required DP matrix entries (if they were not already computed)
	 * returned value is V(x,i)
	 *
	 * i >= 0
	 *
	 */
	double log_likelihood(node_t x, lint i){

		assert(i < nps->size());
		assert(i >= 0);

		double llh = log_x_p_y(V_recursive(x,i),E_recursive(x,i));

		assert(llh<=0);

		//sum of V and E probabilities
		return llh;

	}

	/*
	 * compute the log-likelihood of input sequence on signal beta_1,...,beta_{i-1}
	 */
	double log_likelihood(string& seq, ulint i=0){

		//default: evaluate sequence on the last signal
		if(i==0) i = nps->size()-1;

		auto root = get_root();
		return build_linear_tree(root, i, seq, 0);

	}

	/*
	 * build and retrieve precomputed arrays
	 */
	static precomputed_tables get_precomputed_tables(){

		precomputed_tables pt;

		pt.log_1_e_ = vector<double>(log_1e_size);

		for(ulint i=0;i<log_1e_size;++i){

			double x = -(double)i*step;
			pt.log_1_e_[i] = log(1+exp(x));

		}

		return pt;

	}

	uint8_t get_k_value(){ return k; }

	/*
	 * all nodes below this tree depth (included) should be assessed before to start pruning.
	 * In this case (ONT), we have to reach at least depth k since shallower nodes do not
	 * have a defined score (which is always -infinity). We choose k+1, so to take into
	 * account at least one signal emission.
	 */
	uint8_t get_exhaustive_depth(){
		return k+1;
	}

private:

	double build_linear_tree(node_t parent, ulint i, string& seq, ulint pos){

		if(pos<seq.size()){

			auto x = get_child(parent,char_to_base(seq[pos]));
			return build_linear_tree(x, i, seq, pos+1);

		}else{

			double llh = log_likelihood(parent, i-1);
			delete_linear_tree(parent);

			return llh;

		}

		return 0;

	}

	void delete_linear_tree(node_t x){

		if(x->is_root()){

			delete x;

		}else{

			auto p = x->get_parent();
			delete x;
			delete_linear_tree(p);

		}

	}

	void build(nanopore_signal * nps, uint bandwidth, parameters_t& parameters, precomputed_tables& precomputed_t){

		log_1e_array = precomputed_t.log_1_e_.data();

		this->nps = nps;
		this->strand_ = nps->get_strand();

		if(strand_==TEMPLATE)
			bases_per_event = parameters.bases_per_template_event;
		else
			bases_per_event = parameters.bases_per_complement_event;

		init_transition_prob(parameters);

		k = parameters.k;

		assert(k>0);

		assert(bandwidth>0);

		//min between bandwidth and m
		this->bandwidth = bandwidth;//(bandwidth>number_of_events?number_of_events:bandwidth);

	}

	/*
	 * recursive functions to compute E and V
	 * indexes start from 0 and can be negative (base cases)
	 */
	double E_recursive(node_t x, lint i){

		if(x->depth()<k or i<0) return neg_inf;

		if(is_defined(E(x,i))) return E(x,i);

		double p_ee1 = nps->log_P(nps->get_event(i),x->get_kmer());
		double p_ee2 = nps->log_P(nps->get_event(i),x->get_kmer());

		//base case
		if(x->depth()==k and i==0){

			E(x,i) = -k*log_sigma + p_ee2;

			assert(E(x,i)<=0);
			return E(x,i);

		}

		//else compute recursively E(x,i)

		//Forward algorithm
		double logC = log_p_ve(get_parent(x)->get_kmer()) + p_ee2 + V_recursive(get_parent(x),i);
		double logD = log_p_ee2(get_parent(x)->get_kmer()) + p_ee2 + E_recursive(get_parent(x),i-1);
		double logE = log_p_ee1(x->get_kmer()) + p_ee1 + E_recursive(x,i-1);
		double logF = log_x_p_y(logC,logD);
		E(x,i) = -log_sigma + log_x_p_y(logF,logE);

		//Viterbi algorithm
		/*double logA = log_p_ve(get_parent(x)->get_kmer()) + p_ee2 + V_recursive(get_parent(x),i);
		double logB = log_p_ee2(get_parent(x)->get_kmer()) + p_ee2 + E_recursive(get_parent(x),i-1);
		double logC = log_p_ee1(x->get_kmer()) + p_ee1 + E_recursive(x,i-1);
		E(x,i) = -log_sigma + max(logA,logB,logC);*/

		assert(E(x,i)<=0);//probability!

		return E(x,i);

	}

	double V_recursive(node_t x, lint i){

		//base case 1 : not reachable configurations
		if(x->depth()<k or i<0)	return neg_inf;

		//if the quantity is defined, return it.
		if(is_defined(V(x,i))) return V(x,i);

		//Forward algorithm
		double logA = log_p_vv(get_parent(x)->get_kmer()) + V_recursive(get_parent(x),i);
		double logB = log_p_ev(get_parent(x)->get_kmer()) + E_recursive(get_parent(x),i-1);
		V(x,i) = -log_sigma + log_x_p_y(logA,logB);

		//Viterbi algorithm
		/*double logA = log_p_vv(get_parent(x)->get_kmer()) + V_recursive(get_parent(x),i);
		double logB = log_p_ev(get_parent(x)->get_kmer()) + E_recursive(get_parent(x),i-1);
		double maxab = std::max(logA,logB);
		V(x,i) = -log_sigma + maxab;*/

		assert(V(x,i)<=0);//probability!

		return V(x,i);

	}

	node_t get_parent(node_t x){

		node_t p = x->get_parent();

		assert(x->is_root() or x->depth() == p->depth()+1);

		return p;

	}

	/*
	 * E and V entries of the DP matrix. These are already log-likelihoods (in the description document these quantities are called log E and log V)
	 *
	 */
	double& E(node_t x, lint i){

		assert(i>=0);
		assert(x->depth()>=k);
		assert(i<nps->size());

		return x->E(i);

	}

	double& V(node_t x, lint i){

		assert(i>=0);
		assert(x->depth()>=k);
		assert(i<nps->size());

		return x->V(i);

	}

	/*
	 * returns logx + log(1+exp(logy-logx)) = log(x+y)
	 */
	double log_x_p_y(double logx, double logy){

		//x and y are probabilities
		assert(logx<=0);
		assert(logy<=0);

		//particular cases
		if(logy==neg_inf) return logx; //in this case y=0
		if(logx==neg_inf) return logy; //in this case x=0

		//else: use function log_1e_eff
		double result = logx + log_1e_eff(logy-logx);

		//result is log of a probability
		assert(result <= 0);

		return result;

	}

	/*
	 * returns log(1+exp(x)). Works with any x
	 */
	double log_1e_eff(double x){

		assert(x<inf and x>neg_inf);

		double result = 0;

		if(x>=0){

			//problem: if x is large, exp(x) overflows!!
			//solution: if x>=0, we use the alternative formula log(1+exp(x)) = x + log( 1+e^-x )

			result = x + log_1e_eff_neg(-x);

		}else{

			//here x<0
			result = log_1e_eff_neg(x);

		}

		assert(result<inf and result>neg_inf);

		return result;

	}

	/*
	 * returns log(1+exp(x)), assumes x<=0
	 */
	double log_1e_eff_neg(double x){

		assert(x<=0);

		ulint idx = (ulint)(-x/step);

		//value not cached
		if(idx>=log_1e_size) return log(1+exp(x));

		//value cached
		return log_1e_array[idx];

	}

	void init_transition_prob(hmm_parameters& parameters){

		//pre-computed logarithms, template strand

		log_p_ve_ = vector<double>(parameters.number_of_kmers);
		log_p_vv_ = vector<double>(parameters.number_of_kmers);
		log_p_ev_ = vector<double>(parameters.number_of_kmers);
		log_p_ee1_ = vector<double>(parameters.number_of_kmers);
		log_p_ee2_ = vector<double>(parameters.number_of_kmers);

		uint64_t km = 0;

		auto param_vec = (strand_==TEMPLATE?parameters.temp:parameters.comp);

		for( auto p : param_vec ){

			log_p_ve_[km] = log(p.ve);
			log_p_vv_[km] = log(p.vv);
			log_p_ev_[km] = log(p.ev);
			log_p_ee1_[km] = log(p.ee1);
			log_p_ee2_[km] = log(p.ee2);

			assert(log_p_ve_[km] <= 0 and log_p_ve_[km]>neg_inf);
			assert(log_p_vv_[km] <= 0 and log_p_vv_[km]>neg_inf);
			assert(log_p_ev_[km] <= 0 and log_p_ev_[km]>neg_inf);
			assert(log_p_ee1_[km] <= 0 and log_p_ee1_[km]>neg_inf);
			assert(log_p_ee2_[km] <= 0 and log_p_ee2_[km]>neg_inf);

			km++;

		}

	}

	//to speed up computation, we pre-compute all logarithms of the probabilities

/*	double log_b(kmer km){ return log_b_[km]; } //log(b(km))
	double log_e(kmer km){ return log_e_[km]; } //log(e(km))
	double log_n(kmer km){ return log_n_[km]; } //log(1-(e(km)+b(km)))
*/


/*
	static constexpr double pve = 0.929312;
	static constexpr double pvv = 0.0706881;
	static constexpr double pev = 0.0594006;
	static constexpr double pee1 = 0.267851;
	static constexpr double pee2 = 0.672748;
*/

	/*static constexpr double pvv = 0.01;
	static constexpr double pve = 1-pvv;
	static constexpr double pev = 0.05;
	static constexpr double pee2 = 0.8;
	static constexpr double pee1 = 1-(pev+pee2);


	double log_p_ve(kmer km){ return log(pve); }
	double log_p_vv(kmer km){ return log(pvv); }
	double log_p_ev(kmer km){ return log(pev); }
	double log_p_ee1(kmer km){ return log(pee1); }
	double log_p_ee2(kmer km){ return log(pee2); }*/

	double log_p_ve(kmer km){ return log_p_ve_[km]; }
	double log_p_vv(kmer km){ return log_p_vv_[km]; }
	double log_p_ev(kmer km){ return log_p_ev_[km]; }
	double log_p_ee1(kmer km){ return log_p_ee1_[km]; }
	double log_p_ee2(kmer km){ return log_p_ee2_[km]; }

	//pre-computed logarithms, template strand
	vector<double> log_p_ve_;
	vector<double> log_p_vv_;
	vector<double> log_p_ev_;
	vector<double> log_p_ee1_;
	vector<double> log_p_ee2_;

	//precomputed values of log(1+exp(x)), for negative x
	static constexpr ulint log_1e_size = 1000000;//number of sampled values
	static constexpr double step = 0.001;//sample step

	//while assessing called sequence size, add this number of bases
	//to the max sequence length
	static constexpr double extra_bases = 100;

	precomputed_tables pt;
	double * log_1e_array;

	//the nanopore events together with the normal distributions
	nanopore_signal * nps;

	//number of kmers in the sequence. Equals max sequence length minus k + 1
	//uint number_of_kmers_in_sequence = 0;

	uint bandwidth;

	double bases_per_event = 0;

	//kmer length
	uint8_t k = 5;

	strand strand_;

};

}

#endif /* ON_HMM_HPP_ */
