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
 * kmer.hpp
 *
 *  Created on: Jun 12, 2015
 *      Author: nicola
 *
 *  a kmer, with intuitive shift operations and automatic cast to numeric value
 *
 */

#ifndef KMER_HPP_
#define KMER_HPP_

#include <definitions.hpp>
#include <string>

using namespace std;

namespace ontrc{

class kmer{

public:

	/*
	 * empty constructor: initialize kmer as AAAAA
	 */
	kmer(uint8_t k = 5){

		assert(k<64);

		this->k = k;
		MASK = (ulint(1)<<(2*k))-1;

	}

	/*
	 * constructor from string: convert array of char to a kmer object
	 */
	kmer(string km_c){

		uint8_t k = km_c.size();

		this->k = k;
		MASK = (ulint(1)<<(2*k))-1;

		val = 0;

		for(ulint i=0;i<k;++i)
			val = val*4 + char_to_base(km_c[i]);

	}

	/*
	 * get integer value associated with the kmer
	 */
	operator int() const { return val; }

	/*
	 * convert k-mer to string, for visualization purposes
	 */
	string to_string(){

		string k_str;

		for(uint8_t i=0;i<k;++i)
			k_str += base_to_char(base_at(i));

		return k_str;

	}

	/*
	 * get base at position i (0 being the leftmost base)
	 */
	base operator[](uint8_t i){

		assert(i<k);
		return base_at(i);

	}

	/*
	 * add a base to the right, and shift other bases
	 */
	kmer &operator+=(base b){

		val = ((val<<2) + ulint(b)) & MASK;
		return *this;

	}

	/*
	 * increment value by 1: useful to iterate over all possible kmers
	 */
	kmer &operator++(int){

		val++;
		return *this;

	}

	/*
	 * increment value by 1: useful to iterate over all possible kmers
	 */
	bool operator==(kmer km){

		return this->val == km.val;

	}

	uint8_t get_k(){return k;}

	/*
	 * returns minimum right-shift of input kmer such that this kmer and input kmer overlap
	 * returns k if there is no overlap
	 *
	 * ex:
	 *
	 * ACTGG.overlap(ACTGG) returns 0
	 * ACTGG.overlap(CTGGC) returns 1
	 * ACTGG.overlap(TGGCT) returns 2
	 *
	 */
	uint8_t overlaps(kmer right){

		bool ov = *this==right;

		kmer copy = *this;
		uint8_t len = 0;

		while(not ov){

			len++;

			right.val = right.val >> 2;
			copy.val = ( (copy.val << (2*len)) & copy.MASK ) >> (2*len);

			 ov = copy==right;

		}

		return len;

	}

	kmer rev_compl(){

		kmer km(k);

		for(uint8_t i = 0;i<k;++i)
			km += complement(base_at(k-i-1));

		return km;

	}

	kmer get_complement(){

		kmer km(k);

		for(uint8_t i = 0;i<k;++i)
			km += complement(base_at(i));

		return km;

	}

	/*
	 * returns the chain of kmers starting from this kmer and ending in last
	 */
	vector<kmer> chain(kmer last){

		auto d = overlaps(last);
		vector<kmer> chain_;

		kmer start = *this;

		chain_.push_back(start);

		uint8_t idx_r = k - d;
		while(chain_.size()<d){

			start += last[idx_r];
			chain_.push_back(start);
			idx_r++;

		}

		if(d>0) chain_.push_back(last);

		return chain_;

	}

private:

	base base_at(uint8_t i){

		assert(i<k);
		return (base)((val >> (2*(k-i-1)) ) & ulint(3));

	}

	uint8_t k = 5;
	ulint val = 0;
	ulint MASK = 0;

};

}

#endif /* KMER_HPP_ */
