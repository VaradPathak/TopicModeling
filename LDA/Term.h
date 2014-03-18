/*
 * Term.h
 *
 *  Created on: Mar 14, 2014
 *      Author: Varad Pathak, Arul Samuel
 */

#ifndef Term_H_
#define Term_H_

#include <omp.h>
#include <vector>
#include <string>

using namespace std;

class Term {
public:
	int wordId;
	string word;
	vector<double> *prob;//Probability of K topics
	Term(int wordId, string vocabWord);
};

#endif /* Term_H_ */
