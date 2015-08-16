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
 * fmi_trie.hpp
 *
 *  Created on: Jun 24, 2015
 *      Author: nicola
 *
 *  Brief:
 *
 *  TotalRecaller accepts as corpus index any implementation of a trie of sequences; the trie represents all strings that match with the
 *  corpus, and could be, e.g. a suffix tree, a FM index, a stringome, a LZ-index, a generalized CSA, ecc...
 *
 *  A trie implementation must provide these functions:
 *
 *		- root node: get_root()
 *  	- descend to a children node (edges are labeled with DNA bases): get_child(node_t n, base b)
 *  	- number of leafs in the subtree rooted in a node (i.e. # of occurrences of the string associated with the node): count(node_t n)
 *  	- occurrences of the string associated with the node occ(node_t n)
 *  	- size of the corpus: corpus_size()
 *
 * 	Moreover, a trie must define the types
 *
 * 		- node_t : type of a trie's nodes
 * 		- coordinate : type of a coordinate in the corpus (e.g. an integer text position or a path in a graph)
 *
 *  This class implements the trie using an FM index over the reversed text
 *
 */

#ifndef FMI_TRIE_HPP_
#define FMI_TRIE_HPP_

#include <definitions.hpp>

using namespace sdsl;
using namespace std;

namespace ontrc{

class fmi_trie{

public:

	//in this implementation of the corpus trie, a node is simply an interval [first, second] on the BWT, plus the node's depth
	struct node_t{

		ulint first;
		ulint second;
		ulint depth;
		base b;

	};

	//start position identifier for a node. In this case, start position on the genome
	typedef ulint node_start_id;

	static const ulint default_node_start_id = 0;

	//type of a substring coordinate in the corpus. In this case is simply a pair <begin, length>
	//more complex structures (e.g. a labeled graph) could be a more complex type
	typedef pair<ulint,ulint> coordinate;

	/*
	 * Empty constructor: if called, fm index needs to be loaded from file
	 */
	fmi_trie(){}

	/*
	 * string constructor: build FM index of input text (not a file path!)
	 */
	fmi_trie(string &input_text){

		//reverse string and build FM index on the reverse

		string rev;

		for(ulint i=0;i<input_text.size();++i)
			rev += input_text[input_text.size()-i-1];

		construct_im(fmi, rev.c_str(), 1);

	}

	/*
	 * root node: get full BWT range. Ranges are inclusive!
	 */
	node_t get_root(){

		assert(not empty());

		//inclusive range
		return {0,bwt_size()-1,0,base_end};

	}

	/*
	 * \param n current node
	 * \param c character
	 * \return node associated with string s(n)c
	 */
	node_t get_child(node_t& x, base b){

		//LF function

		auto c = base_to_char(b);

		ulint l = x.first;
		ulint r = x.second;

		if(l > r) return {1,0,x.depth+1,b};

		assert(r < fmi.size());

	    ulint c_begin = fmi.C[fmi.char2comp[c]];

	    ulint c_before_l = fmi.bwt.rank(l, c); // count c in bwt[0..l-1]
	    ulint c_before_r = fmi.bwt.rank(r+1, c); // count c in bwt[0..r]

	    assert(c_before_r >= c_before_l);

	    ulint n_occ = c_before_r-c_before_l;

	    if(n_occ==0) return {1,0,x.depth+1,b};

	    ulint l_start = c_begin+c_before_l;

	    return {l_start, l_start+n_occ-1, x.depth+1,b};

	}

	/*
	 * return all occurrences of the string s(n), i.e. the string associated
	 * with node n of the trie
	 *
	 * \param n the trie node
	 * \return occurrences of s(n) in the corpus
	 *
	 */
	vector<coordinate> occ(node_t& x){

		vector<coordinate> result;

		for(ulint i=x.first;i<=x.second;++i){

			assert(corpus_size() >= fmi[i]);
			assert(corpus_size() - fmi[i] >= x.depth);

			ulint o = (corpus_size() - fmi[i]) - x.depth;

			result.push_back({o,x.depth});

		}

		return result;

	}

	/*
	 * return number of leafs reachable from node n, i.e. number of occurrences of the current string s(n) in the corpus
	 */
	ulint count(node_t& x){

		if(x.second<x.first) return 0;
		return (x.second - x.first) + 1;

	}

	/*
	 * return logarithm of empirical frequency of a string (represented by a node) in the corpus
	 *
	 * alternatively, this is the log-likelihood of extracting randomly sequence s(x) from the
	 * corpus, given that the extraction is uniform
	 *
	 */
	double log_frequency(node_t& x){

		double nr_occ = count(x); //number of occurrences of the string

		if(x.depth==0) return neg_inf;//root

		return log(nr_occ/(corpus_size()-x.depth+1));

	}

	/*
	 * save the structure to the specified path.
	 * \param path_prefix prefix of the index files. suffix ".sdsl.fmi" is automatically added
	 */
	void save_to_disk(string path_prefix){

		assert(not empty());

		string path = string(path_prefix).append(get_extension());
		store_to_file(fmi, path);

	}

	/*
	 * load the structure from the specified path.
	 * \param path_prefix prefix of the index files. suffix ".sdsl.fmi" is automatically added
	 */
	void load_from_disk(string path_prefix){

		string path = string(path_prefix).append(get_extension());
		load_from_file(fmi, path);
		assert(not empty());

	}

	/*
	 * return size (=number of characters) of the corpus
	 */
	ulint corpus_size(){

		if(bwt_size()==0) return 0;	//index has not been constructed yet
		return bwt_size()-1;

	}

	/*
	 * return file extension used
	 */
	string get_extension(){
		return string(".sdsl.fmi");
	}

	bool empty(){ return corpus_size()==0; }

	node_start_id get_start_id(node_t& x){

		if(x.second<x.first) return 0;

		assert(corpus_size() >= fmi[x.first]);
		assert(corpus_size() - fmi[x.first] >= x.depth);

		ulint o = (corpus_size() - fmi[x.first]) - x.depth;

		return o;

	}

	ulint get_depth(node_t& x){return x.depth;}


	/*
	 * searches and returns all length-len substrings in the corpus which are suffixed
	 * by suf
	 */
	vector<string> retrieve_substrings_by_suffix(string suf, ulint len){

		assert(suf.size()<=len);

		vector<string> result;

		string rev_suf = reverse(suf);//because we index the reverse sequence

		auto locations = sdsl::locate(fmi, rev_suf.begin(), rev_suf.begin()+rev_suf.size());

		for(auto l : locations){

			if(l+len<corpus_size()){

				string s = sdsl::extract(fmi, l, l+len-1);

				if(not contains_separator(s))
					result.push_back(reverse(s));

			}

		}

		return result;

	}

private:

	ulint bwt_size(){

		return fmi.size();

	}

	fm_index_t fmi;	//the underlying FM index (sdsl)

};

}


#endif /* FMI_TRIE_HPP_ */
