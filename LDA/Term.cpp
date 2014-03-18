/*
 * Term.cpp
 *
 *  Created on: Mar 14, 2014
 *      Author: Varad Pathak, Arul Samuel
 */

#include "Term.h"
#include "SCVB0.h"
#include <string>
#include <vector>

Term::Term(int id, string vocabWord) {
	wordId = id;
	word = vocabWord;
	prob = new vector<double>();
}

