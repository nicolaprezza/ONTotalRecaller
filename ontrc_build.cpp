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
 *  build FM index for ontrc
 *
 */

#include <string>
#include <iostream>

#include <definitions.hpp>
#include <fmi_trie.hpp>

using namespace ontrc;

void help(){
	cout << "ontrc_build : build FM index used by ONTRC." << endl << endl;
	cout << "Usage: ontrc_build <fasta> <output prefix>" << endl;
	cout <<	"   <fasta>          fasta file" << endl;
	cout <<	"   <output prefix>  name of the index. Suffix '.sdsl.fmi' is automatically added." << endl;
	exit(0);
}


int main(int argc, char** argv){

	if(argc!=3) help();

	string fasta_file(argv[1]);
	string out_prefix(argv[2]);

	//load reference genome
	cout << "Loading reference from file " << fasta_file << " ... " << flush;
	string forward_strand = read_fasta_file(fasta_file);
	cout << "done." <<endl;

	//build reference text: forward strand concatenated to its reverse complement, plus a separator in the middle
	string reference = forward_strand + '$' + ontrc::reverse_complement(forward_strand);
	//string reference = forward_strand + '$';

	ontrc::remove_Ns(reference);

	cout << "Building FM index ... " << flush;
	fmi_trie fmi(reference);
	cout << "done."  << endl;

	cout << "Saving FM index ... " << flush;
	fmi.save_to_disk(out_prefix);
	cout << "done."  << endl;

}

















