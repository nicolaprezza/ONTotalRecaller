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
 * total_recaller.hpp
 *
 *  Created on: Jun 24, 2015
 *      Author: nicola
 *
 * This class implements TotalRecaller algorithm is such a way that is independent from the underlying
 * corpus implementation and scoring strategy. See classes fmi_trie and on_hmm for more details on the
 * template models accepted as parameters.
 *
 */

#ifndef TOTAL_RECALLER_HPP_
#define TOTAL_RECALLER_HPP_

#include <definitions.hpp>
#include <fmi_trie.hpp> //trie of corpus substrings implemented as an FM index
#include <on_hmm.hpp>	//scoring algorithm for Oxford nanopore signals

using namespace std;

namespace ontrc{

template<
	class corpus_index = fmi_trie,		//default: FM index
	class scoring_strategy = on_hmm		//default: oxford nanopore scoring scheme (HMM)
>
class total_recaller{

public:

	//this is the input signal type. TRC aligns a signal_t object against the corpus
	//using the scoring function defined in scoring_strategy
	typedef typename scoring_strategy::signal_t signal_t;

	//TRC navigates in parallel 2 trees: the one associated with the base caller and the corpus suffix trie.

	//this is a node of the complete tree associated with the underlying sequencing technology (e.g. Oxford Nanopore).
	//This tree represents the base caller.
	//note: this node is a pointer, and must be freed inside this class
	typedef typename scoring_strategy::node_t node_caller_t;

	//corpus suffix trie node.
	//note: this node is not implemented as a pointer, and thus does not have to be freed inside this class.
	typedef typename corpus_index::node_t node_corpus_t;

	//a TRC node
	class node_t{

		using node_ptr = node_t*;

	public:

		/*
		 * constructor.
		 */
		node_t( node_caller_t node_caller,
				node_corpus_t node_corpus,
				node_ptr parent){

			base b = node_caller->get_base();

			if(parent!=NULL){

				corpus_kmer = parent->corpus_kmer;
				assert(corpus_kmer.get_k()==kmer_length);

				corpus_kmer += b;
				assert(corpus_kmer.get_k()==kmer_length);

			}else{

				corpus_kmer = kmer(kmer_length);
				assert(corpus_kmer.get_k()==kmer_length);

			}

			this->node_caller = node_caller;
			this->node_corpus = node_corpus;
			this->parent = parent;

			has_child_ = vector<bool>(4,false);
			children = vector<node_ptr>(4,NULL);

		}

		bool has_child(base b){
			return has_child_[b];
		}

		/*
		 * assign c as child of this node on edge s
		 */
		void assign_child(node_ptr c){

 			//x must not already have child b
			assert(not has_child(c->get_corpus_base()));

			children[c->get_corpus_base()] = c;
			has_child_[c->get_corpus_base()] = true;

			nr_of_children_++;

		}

		/*
		 * free memory of child b and decrease children counter
		 */
		void delete_child(base b){

			assert(has_child(b));
			assert(nr_of_children()>0);

			assert(children[b]!=NULL);
			assert(children[b]->get_parent()==this);

			delete children[b];
			children[b] = NULL;

			has_child_[b] = false;
			nr_of_children_--;

		}

		node_ptr get_child(base b){

			assert(has_child(b));

			return children[b];

		}

		bool is_leaf(){return nr_of_children()==0;}

		bool is_root(){return parent==NULL;}

		int nr_of_children(){return nr_of_children_;}

		node_ptr get_parent(){return parent;}

		ulint depth(){

			assert(node_corpus.depth==node_caller->depth());

			return node_caller->depth();

		}

		kmer get_corpus_kmer(){

			return corpus_kmer;

		}

		base get_caller_base(){return node_caller->get_base();}
		base get_corpus_base(){return node_corpus.b;}

		void set_caller_base(base b){ node_caller->set_base(b);}


		node_caller_t node_caller;
		node_corpus_t node_corpus;

		node_ptr parent;

		typename corpus_index::node_start_id start_id = corpus_index::default_node_start_id;
		double score = neg_inf;

		vector<node_ptr> children;
		vector<bool> has_child_;

		int nr_of_children_=0;

		//last kmer of the caller's node
		kmer corpus_kmer;

	};

	//a pointer to a TRC tree node
	typedef node_t* node_ptr;

	//precomputed arrays used by the scoring strategy.
	typedef typename scoring_strategy::precomputed_tables precomputed_tables;

	//parameters of the scoring strategy
	typedef typename scoring_strategy::parameters_t parameters_t;

	//empty constructor
	total_recaller(){}

	/*
	 * Constructor: corpus is passed from outside since it can be any kind of data structure encoding a set of strings
	 */
	total_recaller(corpus_index * corpus, parameters_t parameters){

		this->corpus = corpus;
		this->parameters = parameters;
		pt = scoring_strategy::get_precomputed_tables();

	}

	/*
	 * align a signal against the corpus & perform base calling at the same time. Then, return all matches
	 * against the reference, in decreasing order of likelihood.
	 * This function implements the Branch&Bound algorithm.
	 *
	 * returns a vector of 2  strings: the match on reference and the called bases (reference unbiased)
	 *
	 */
	vector<pair<string, alignment>> call(signal_t& s, double W, bool verbose = false, bool free_mem = true){

		this->W = W;
		this->verbose=verbose;

		assert(corpus != NULL);

		//create scoring_strategy class
		scoring_strategy sc_st(&s, band_width, parameters,pt);

		//create roots of the 2 trees
		auto root_w0 = get_root(sc_st);
		root_w1 = get_root(sc_st);

		//fill the vector 'candidates' with all nodes up to the exhaustive depth
		vector<node_ptr> candidates_w0 = vector<node_ptr>(); // empty vector
		vector<node_ptr> candidates_w1 = vector<node_ptr>(); // empty vector

		fill_candidates_up_to_exhaustive_depth(sc_st, candidates_w0, root_w0, sc_st.get_exhaustive_depth());

		int step = 5;
		int last_perc = -step;
		int perc;
		if(verbose) cout << endl;

		//core of the B&B algorithm: iterate over the signals and, at each step, create and
		//select new candidate nodes
		for(ulint i=0;i<s.size();++i){

			//call weighted B&B
			branch_and_bound(sc_st,candidates_w1,i,W);

			//do B&B without reference.
			node_ptr best_node = branch_and_bound(sc_st,candidates_w0,i,0);

			kmer best_kmer = best_node->get_corpus_kmer();			//best node's last kmer
			ulint substring_length = best_node->depth();			//best node's depth

			//now locate best kmer
			//and add all candidates to the weighted B&B

			if(substring_length>=kmer_length){

				vector<string> paths = corpus->retrieve_substrings_by_suffix(best_kmer.to_string(), substring_length);

				//turns string paths into tree paths, and adds them in the tree
				vector<node_ptr> new_candidates = paths_to_nodes(paths,sc_st);

				for(auto x : new_candidates){

					x->score = get_bayesian_score(sc_st,x,i,W);

					candidates_w1.push_back(x);

				}

				if(insert_snps)
					for(auto x : candidates_w1)
						insert_one_best_snp(sc_st, x,i);

			}

			//insert one SNP in the best path

			if(verbose){

				perc = (100*i)/s.size();
				if(perc>=last_perc+step){

					cout << " " << perc << "% done" << endl;
					last_perc=perc;

				}

			}

		}

		//cluster nodes that have the same start position, keep the best.
		candidates_w1 = uniq(candidates_w1,1);

		//sort candidates by decreasing score
		auto comp = [](const node_ptr a, const node_ptr b) -> bool
						{

							return a->score > b->score;

						};

		std::sort(candidates_w1.begin(),candidates_w1.end(),comp);

		if(verbose){

			cout << "Position, length, score of candidates:" << endl;

			for(node_ptr c : candidates_w1){

				cout << c->start_id << " " << c->depth() << " " << c->score << endl;

			}

		}

		node_ptr best_node_reference = best_candidate(sc_st,candidates_w1,s.size()-1,W);
		node_ptr best_node_caller = best_candidate(sc_st,candidates_w0,s.size()-1,0);

		//execute SNP strategy only if enabled and if we use reference as prior
		/*if(insert_snps and W>0){

			//for each candidate x, climb (k-1) nodes of the tree from x
			//and substitute the caller base that maximizes x->score
			best_node_reference_with_snps = insert_best_snps(sc_st, best_node_reference,s.size()-1);

			if(verbose){
				cout << "Position, length, score of best candidate after calling SNPs:" << endl;
				cout << best_node_reference_with_snps->start_id << " " << best_node_reference_with_snps->depth() << " " << best_node_reference_with_snps->score << endl;
			}

		}*/

		vector<pair<string, alignment>> calls;

		vector<typename corpus_index::coordinate> occ;

		if(candidates_w1.size()>0)
			occ = corpus->occ(best_node_reference->node_corpus);

		alignment ali;
		if(occ.size()>0){

			auto coord = occ[0];
			ali = global_to_local(coord,corpus->corpus_size());

		}

		calls.push_back({path(best_node_reference),ali});
		calls.push_back({path(best_node_caller),{}});

		if(free_mem) free_memory();

		return calls;

	}


private:

	node_ptr get_ancestor(node_ptr x,ulint i){

		if(i==0 or x->is_root()) return x;
		return get_ancestor(x->get_parent(),i-1);

	}

	/*node_ptr insert_best_snps(scoring_strategy& sc_st, node_ptr x, ulint i){

		string ref_match = path(x);

		//build new linear tree with SNPs
		node_ptr current_node = get_root(sc_st);

		for(ulint j = 0; j<ref_match.size();++j){

			current_node = get_and_assign_child(sc_st, current_node, char_to_base((uchar)ref_match[j]));
			insert_one_best_snp(sc_st, current_node, i);

		}

		current_node->start_id = x->start_id;

		return current_node;

	}*/

	void insert_one_best_snp(scoring_strategy& sc_st, node_ptr x, ulint i){

		if(x==NULL) return;

		uchar k = sc_st.get_k_value();

		//for each candidate in increasing score order

		//we are going to mutate the (k-1)-th base climbing
		//from x

		if(x->depth() >= k){

			//get corpus kmer of x, and extract first base
			base orig_b = x->get_corpus_kmer()[0];

			auto scores = vector<double>(4,0);
			scores[orig_b] = x->score;

			//try all possible substitutions and re-compute score.
			for(base b = A; b!=base_end; b = base(b+1)){

				//do not try again base orig_b
				if(b != orig_b){

					//get (2k-2)-th ancestor, or the root if this node's heigth is < 2k-2
					node_ptr ancestor_2k = get_ancestor(x,2*k-2);

					//get (k-1)-th ancestor
					node_ptr ancestor_k = get_ancestor(x,k-1);

					//mutate base of the (k-1)-th ancestor
					ancestor_k->set_caller_base(b);

					//reset and re-compute all scores in the subtree rooted ancestor_2k:
					//in this way, we take into account the new mutation

					erase_and_re_compute_scores_in_subtree(sc_st, i, ancestor_2k);

					//update current score
					scores[b] = x->score;

				}

			}

			//now find the best base
			double best_score = scores[A];
			base best_base = A;

			for(base b = C; b!=base_end; b = base(b+1)){

				if(scores[b] > best_score){

					best_score = scores[b];
					best_base = b;

				}

			}

			//finally, apply the best mutation

			//get (2k-2)-th ancestor, or the root if this node's height is < 2k-2
			node_ptr ancestor_2k = get_ancestor(x,2*k-2);

			//get (k-1)-th ancestor
			node_ptr ancestor_k = get_ancestor(x,k-1);

			//mutate base of the (k-1)-th ancestor
			ancestor_k->set_caller_base(best_base);

			//reset and re-compute all scores in the subtree rooted ancestor_2k:
			//in this way, we take into account the new mutation

			erase_and_re_compute_scores_in_subtree(sc_st, i, ancestor_2k);

		}

	}

	void erase_and_re_compute_scores_in_subtree(scoring_strategy& sc_st, ulint i, node_ptr x){

		erase_scores_in_subtree(x);
		re_compute_scores_in_subtree(sc_st, i, x);

	}

	void re_compute_scores_in_subtree(scoring_strategy& sc_st, ulint i, node_ptr x){

		x->score = get_bayesian_score(sc_st, x, i, W);

		for(base b = A; b!=base_end; b = base(b+1)){

			if(x->has_child(b)){

				re_compute_scores_in_subtree(sc_st, i, x->get_child(b));

			}

		}

	}

	vector<node_ptr> get_leafs(node_ptr x){

		vector<node_ptr> result;

		//x is leaf: return just {x}
		if(x->is_leaf()){

			result.push_back(x);
			return result;

		}

		//x is not leaf

		for(base b = A; b!=base_end; b = base(b+1)){

			if(x->has_child(b)){

				vector<node_ptr> sub_leafs = get_leafs(x->get_child(b));
				result.push_back(sub_leafs);

			}

		}

		return result;

	}

	void erase_scores_in_subtree(node_ptr x){

		x->node_caller->erase_scores();

		for(base b = A; b!=base_end; b = base(b+1)){

			if(x->has_child(b)){

				erase_scores_in_subtree(x->get_child(b));

			}

		}

	}

	/*
	 * input: set of strings
	 * behavior: traverses the weighted tree (w1) from its
	 * 	root following the strings in the set, creating new
	 * 	nodes if necessary. At the end of each path, return
	 * 	the nodes corresponding to each path.
	 */
	vector<node_ptr> paths_to_nodes(vector<string> paths, scoring_strategy& sc_st){

		vector<node_ptr> nodes;

		assert(root_w1!=NULL);

		//for each path
		for(auto s : paths){

			node_ptr x = root_w1;

			//for each character
			for(auto c : s){

				base b = char_to_base(c);

				//traverse the path labeled b from x. If x already has
				//a b-child, return it. Otherwise, create a brand new.
				x = get_and_assign_child(sc_st, x, b);

			}

			nodes.push_back(x);

		}

		return nodes;

	}

	/*
	 * input: a leaf x
	 * behavior: this procedure implements a B&B strategy which explores
	 * a tree having seq(x) as backbone, and where SNPs are inserted in order to
	 * maximize the score of the final node.
	 *
	 * output: a leaf L such that the path from the root to L includes the best
	 * combination of SNPs
	 *
	 */
/*	node_ptr insert_SNPs(node_ptr x, scoring_strategy& sc_st){

		//create a new tree. This tree is explored inserting
		//SNPs and maximizing the final score
		auto root = get_root(sc_st);

		return NULL;

	}

	//try to insert a SNP in this node?
	bool try_SNP(node_ptr x){

		//return number_of_occurrences(x) == 1;
		return true;

	}*/

	/*
	 * return best candidate at signal index i
	 */
	node_ptr best_candidate(scoring_strategy& sc_st, vector<node_ptr>& candidates, ulint i, double W){

		if(candidates.size()==0) return NULL;

		node_ptr best = candidates[0];
		double best_score_s_e = get_bayesian_score(sc_st, best, i,W);
		double score_e_s = get_caller_score(sc_st, best, i);

		for(auto x : candidates){

			double score_s_e = get_bayesian_score(sc_st,x,i,W);

			if(score_s_e>best_score_s_e){

				best_score_s_e = score_s_e;
				best = x;
				score_e_s = get_caller_score(sc_st, x, i);

			}

		}

		return best;

	}

	/*
	 * return called bases path from root to this node as a string of called bases
	 */
	string path(node_ptr x){

		if(x==NULL) return string();

		if(x->is_root()) return string();

		//return called base
		//(not the base on reference, which is different in case of SNPs)
		return path(x->get_parent()) + (char)base_to_char(x->get_caller_base());

	}

	/*
	 * create children of nodes in vector 'candidates', re-evaluate children and
	 * parent nodes on the new signal (number i), and then keep only the best
	 * beamwidth nodes (re-writing vector candidates)
	 *
	 * \param i signal index
	 *
	 * \return best node, or NULL if candidates set becomes empty
	 */
	node_ptr branch_and_bound(scoring_strategy& sc_st, vector<node_ptr>& candidates, ulint i, double W){

		if(candidates.size()==0) return NULL;

		vector<node_ptr> uniq_candidates;

		{

			vector<node_ptr> next_candidates_vec;

			{

				set<node_ptr> next_candidates;//we use a set to remove duplicates (i.e. same pointer)

				//fill next_candidates with all next candidate nodes

				//for each candidate
				for(auto x : candidates){

					x->score = get_bayesian_score(sc_st,x,i,W);

					if(x->start_id==0 and W>0)
						x->start_id = corpus->get_start_id(x->node_corpus);

					//push the node itself: this is needed in order to take into account
					//oversegmentation (kmer emitting multiple signals)

					if(x->score>neg_inf){ //consider only nodes of non-null probability

						next_candidates.insert(x);

						auto cA = get_and_assign_child(sc_st,x,A);
						auto cC = get_and_assign_child(sc_st,x,C);
						auto cG = get_and_assign_child(sc_st,x,G);
						auto cT = get_and_assign_child(sc_st,x,T);

						cA->score = get_bayesian_score(sc_st,cA,i,W);
						cC->score = get_bayesian_score(sc_st,cC,i,W);
						cG->score = get_bayesian_score(sc_st,cG,i,W);
						cT->score = get_bayesian_score(sc_st,cT,i,W);

						if(W>0){

							cA->start_id = corpus->get_start_id(cA->node_corpus);
							cC->start_id = corpus->get_start_id(cC->node_corpus);
							cG->start_id = corpus->get_start_id(cG->node_corpus);
							cT->start_id = corpus->get_start_id(cT->node_corpus);

						}

						if(cA->score>neg_inf) next_candidates.insert(cA);
						if(cC->score>neg_inf) next_candidates.insert(cC);
						if(cG->score>neg_inf) next_candidates.insert(cG);
						if(cT->score>neg_inf) next_candidates.insert(cT);

					}

				}

				next_candidates_vec = vector<node_ptr>(next_candidates.begin(),next_candidates.end());

			}

			//now filter nodes and put filtered candidates in uniq_candidates
			uniq_candidates = vector<node_ptr>();

			if(W==0 or i<cluster_depth){

				uniq_candidates = vector<node_ptr>(next_candidates_vec.begin(),next_candidates_vec.end());

			}else{

				uniq_candidates = uniq(next_candidates_vec, sub_beam_width);

			}

		}

		//sort candidates according to scores
		std::sort(	uniq_candidates.begin(),
					uniq_candidates.end(),
					[](const node_ptr a, const node_ptr b) -> bool
						{
		    				return a->score > b->score;
						});

		//erase old candidates
		candidates = vector<node_ptr>();

		//now extract best beamwidth nodes

		ulint bw = beam_width;
		for(ulint j = 0;j < std::min(bw,uniq_candidates.size()) ;++j){

			candidates.push_back(uniq_candidates[j]);

		}

		return best_candidate(sc_st, candidates, i, W);

	}


	vector<node_ptr> uniq(vector<node_ptr>& old_candidates, ulint sub_beam_width){

		//2 nodes are considered as equal if they have
		//the same start ID: in this way we can cluster together
		//nodes starting in the same genomic position

		//nodes are sorted by decreasing start_ID and, in case
		//of same start_ID, decreasing score.

		auto comp = [](const node_ptr a, const node_ptr b) -> bool
						{

							if(a->start_id/block_size>b->start_id/block_size) return true;
							if(a->start_id/block_size<b->start_id/block_size) return false;

							return a->score > b->score;

						};

		std::sort(old_candidates.begin(),old_candidates.end(),comp);

		vector<node_ptr> uniq_candidates;

		if(old_candidates.size()>0)
			uniq_candidates.push_back(old_candidates[0]);

		ulint nr_clones = 1;//number of nodes with same start position inserted in uniq_candidates

		for(ulint j=1;j<old_candidates.size();++j){

			//look at the start position. Since nodes with same start ID are ordered
			//by decreasing score, the first sub_beam_width have the highest score

			auto start_id_candidates = old_candidates[j]->start_id/block_size;
			auto start_id_uniq = uniq_candidates[uniq_candidates.size()-1]->start_id/block_size;

			if(start_id_candidates != start_id_uniq){

				nr_clones=1;
				uniq_candidates.push_back(old_candidates[j]);

			}else if(nr_clones<sub_beam_width){

				nr_clones++;
				uniq_candidates.push_back(old_candidates[j]);

			}

		}

		return uniq_candidates;

	}

	/*
	 * if x has child b already, return it. otherwise, create new and return it.
	 */
	node_ptr get_and_assign_child(scoring_strategy& sc_st, node_ptr x, base b){

		return get_and_assign_child(sc_st, x, b, b);

	}

	/*
	 * if x has child b already, return it. otherwise, create new and return it.
	 */
	node_ptr get_and_assign_child(scoring_strategy& sc_st, node_ptr x, base b_corpus, base b_caller){

		if(not x->has_child(b_caller)){
			auto cA = new_child(sc_st,x,b_corpus, b_caller);
			assign_child(x,cA);
			return cA;
		}

		return x->children[b_caller];

	}

	/*
	 * starting from the leafs, free memory of all nodes
	 */
	void free_memory(){

		for(auto l : leaves){

			assert(l->is_leaf());
			free(l);

		}

		//delete allocated caller nodes. By doing it here, we
		//avoid freeing multiple times caller nodes (remember
		//that different node_t objects can share the same
		//caller object
		for(auto x_c : nodes_caller){

			delete x_c;

		}

	}

	//recursively free memory from node x
	void free(node_ptr& x){

		assert(x->is_leaf());

		if(x->is_root()){

			assert(x!=NULL);
			delete x;
			x=NULL;

		}else{

			node_ptr p = x->get_parent();
			auto b = x->get_corpus_base();

			assert(p->has_child(b));

			//if we don't delete a leaf, there will be memory leaks!
			assert(p->children[b]->is_leaf());
			assert(p->children[b] == x);
			assert(p->children[b]->get_parent() == p);

			p->delete_child(b);

			assert(not p->has_child(b));

			//if p has no more children, recursion
			if(p->is_leaf()){

				free(p);

			}

		}

	}

	void fill_candidates_up_to_exhaustive_depth(scoring_strategy& sc_st, vector<node_ptr>& candidates, node_ptr x, ulint exhaustive_depth){

		//create children. No SNPs for now (this is not a big assumption since we do this up do depth k=5)
		auto cA = get_and_assign_child(sc_st,x,A);
		assert(cA->depth()==x->depth()+1);

		auto cC = get_and_assign_child(sc_st,x,C);
		assert(cC->depth()==x->depth()+1);

		auto cG = get_and_assign_child(sc_st,x,G);
		assert(cG->depth()==x->depth()+1);

		auto cT = get_and_assign_child(sc_st,x,T);
		assert(cT->depth()==x->depth()+1);

		candidates.push_back(cA);
		candidates.push_back(cC);
		candidates.push_back(cG);
		candidates.push_back(cT);

		if(x->depth() < exhaustive_depth){

			//call recursively

			fill_candidates_up_to_exhaustive_depth(sc_st,candidates,cA,exhaustive_depth);
			fill_candidates_up_to_exhaustive_depth(sc_st,candidates,cC,exhaustive_depth);
			fill_candidates_up_to_exhaustive_depth(sc_st,candidates,cG,exhaustive_depth);
			fill_candidates_up_to_exhaustive_depth(sc_st,candidates,cT,exhaustive_depth);

		}

	}


	/*
	 * get root of the search tree
	 */
	node_ptr get_root(scoring_strategy& sc_st){

		//get_root() creates a new object
		auto r_caller = sc_st.get_root();
		auto r_corpus = corpus->get_root();

		nodes_caller.insert(r_caller);

		auto root = new node_t(r_caller, r_corpus, NULL);

		leaves.insert(root);

		return root;

	}

	/*
	 * new child constructor. NO SNP: caller and corpus base coincide
	 *
	 * build and get child node. Note that memory must be freed afterwards!
	 * this function does not assign the child to its parent!
	 */
	node_ptr new_child(scoring_strategy& sc_st, node_ptr x, base b){

		//call the SNP constructor with 2 equal bases
		return new_child(sc_st, x, b, b);

	}

	/*
	 * new child constructor. insert SNP: caller and corpus base are different!
	 *
	 * build and get child node. Note that memory must be freed afterwards!
	 * this function does not assign the child to its parent!
	 */
	node_ptr new_child(scoring_strategy& sc_st, node_ptr x, base b_corpus, base b_caller){

		//get_child() creates a new object
		//create node caller
		auto x_caller = sc_st.get_child(x->node_caller,b_caller);

		//create node corpus
		auto x_corpus = corpus->get_child(x->node_corpus,b_corpus);

		assert(x_corpus.depth == x->node_corpus.depth+1);

		nodes_caller.insert(x_caller);

		node_ptr child = new node_t(x_caller, x_corpus,x);

		assert(child->depth() == x->depth()+1);

		//do not assign child to parent

		return child;

	}

	/*
	 * create the edge x->c and update leafs vector
	 */
	void assign_child(node_ptr x, node_ptr c){

		bool was_leaf = x->is_leaf();
		x->assign_child(c);

		//since we just created it, child node is a leaf
		leaves.insert(c);

		//if x was leaf, now it is not. remove x from vector leaves
		if(was_leaf){

			assert(leaves.find(x) != leaves.end()); //x was a leaf

			ulint size = leaves.size();
			leaves.erase(x);
			assert(leaves.size()==size-1);

		}

	}

	ulint number_of_occurrences(node_ptr x){

		return corpus->count(x->node_corpus);

	}

	/*
	 * get log P(node|signal) <proportional to> log P(signal|node) + W * log P(node)
	 *
	 * we use a weight W in order to determine the importance of the reference in the base calling process
	 *
	 */
	double get_bayesian_score(scoring_strategy& sc_st, node_ptr x, ulint i, double W){

		//if W==0, multiplication 0*-infinity produces a -nan in debug mode
		if(W==0) return sc_st.log_likelihood(x->node_caller,i);

		return sc_st.log_likelihood(x->node_caller,i) + W*corpus->log_frequency(x->node_corpus);

	}

	/*
	 * get log P(node|signal) <proportional to> log P(signal|node) + W * log P(node)
	 *
	 * we use a weight W in order to determine the importance of the reference in the base calling process
	 *
	 */
	double get_caller_score(scoring_strategy& sc_st, node_ptr x, ulint i){

		return sc_st.log_likelihood(x->node_caller,i);

	}

	ulint get_depth(node_ptr x){

		return corpus->get_depth(x->node_corpus);

	}

	//all substrings that match with the corpus
	corpus_index * corpus = NULL;

	precomputed_tables pt;

	parameters_t parameters;

	double W = 1;

	//this is the length of kmers used to
	//find matches with the reference.
	static constexpr uchar kmer_length = 12;

	static constexpr ulint band_width = 200; //width of banded DP matrix

	static constexpr ulint beam_width = 100; //beam width of the branch&bound algorithm
	//start clustering nodes from this event depth. clustering is performed according to start
	//position of the nodes. If multiple nodes have the same start position, we keep
	//only the node with highest score among them
	static constexpr ulint cluster_depth = 40;
	//if more than one node have the same start position, keep
	//the best sub_beam_width ones (highest scores) with that
	//start position
	static constexpr ulint sub_beam_width = beam_width;

	//reference is divided in blocks of block_size bases. if 2 sequences fall inside the
	//same block, use the one with the highest score
	static constexpr ulint block_size = 1;

	//we keep leafs here in order to be able to free memory in the end
	//we use a set in order to efficiently retrieve and delete nodes that
	//are no more a leaf because we created a child
	set<node_ptr> leaves;

	//we store here all caller nodes (which must be implemented as pointers) in order to be able
	//to free memory after the alignment
	set<node_caller_t> nodes_caller;

	bool verbose;

	//try to insert SNPs in order to maximize score?
	bool insert_snps = false;

	//zero-weight caller (non reference-biased):

	node_ptr root_w1 = NULL;

};

typedef total_recaller<fmi_trie, on_hmm> on_total_recaller;

}

#endif /* TOTAL_RECALLER_HPP_ */
