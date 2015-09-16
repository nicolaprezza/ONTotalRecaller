/*
 * total_recaller_snps.hpp
 *
 *  Created on: Sep 16, 2015
 *      Author: nicola
 */

#ifndef TOTAL_RECALLER_SNPS_HPP_
#define TOTAL_RECALLER_SNPS_HPP_

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
 * This class implements TotalRecaller SNP calling algorithm
 *
 */


#include <definitions.hpp>
#include <fmi_trie.hpp> //trie of corpus substrings implemented as an FM index
#include <on_hmm.hpp>	//scoring algorithm for Oxford nanopore signals

using namespace std;

namespace ontrc{

/*
 * input:
 * 	- an alignment and its corresponding DNA sequence
 * 	- a set of variable positions on the reference (forward strand). Must be ordered by increasing coordinate!
 * 	- HMMs for forward and complement strands
 *
 * output: the input DNA sequence is modified to accept only SNPs that
 * increase the HMM score.
 *
 */
template<
	class scoring_strategy = on_hmm		//default: oxford nanopore scoring scheme (HMM)
>
static pair<string, alignment> call_snps(pair<string, alignment> input, set<ulint> snp_positions, scoring_strategy &fwd_hmm, scoring_strategy &rev_hmm){

	string forward_call;
	string rev_call;
	strand s = input.second.s;

	//inclusive start/end of alignment
	ulint start = input.second.start;
	ulint length = input.second.length;
	ulint end = input.second.start + length -1;

	if(s == TEMPLATE){

		forward_call = input.first;
		rev_call = reverse_complement(forward_call);

	}else{

		rev_call = input.first;
		forward_call = reverse_complement(rev_call);

	}

	//extract SNPS that intersect the alignment position
	for(	set<ulint>::iterator it = snp_positions.lower_bound(start);
			it != snp_positions.end() and *it <= end;
			++it	){

		ulint i = *it - start;

		vector<double> scores_fwd(4,0);
		vector<double> scores_rev(4,0);

		for(base b = A; b != base_end; ++b){

			forward_call[i] = base_to_char(b);
			rev_call[length - i - 1] = base_to_char(complement(b));

			scores_fwd[b] = fwd_hmm.log_likelihood(forward_call);
			scores_rev[b] = rev_hmm.log_likelihood(rev_call);

		}

		//maximize score
		base best_snp = A;
		double max_score = combine_scores(scores_fwd[A],scores_rev[A]);

		for(base b = C;b != base_end;++b){

			double best_score = combine_scores(scores_fwd[best_snp],scores_rev[best_snp]);
			double this_score = combine_scores(scores_fwd[b],scores_rev[b]);

			if( this_score >  best_score) best_snp = b;

		}

		forward_call[i] = base_to_char(best_snp);
		rev_call[length - i - 1] = base_to_char(complement(best_snp));

	}

	if(s==TEMPLATE)
		return {forward_call, input.second};

	return {rev_call, input.second};

}

double combine_scores(double a, double b){

	return a+b;

}

}

#endif /* TOTAL_RECALLER_SNPS_HPP_ */
