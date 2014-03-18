/*
 * LDA.h
 *
 *  Created on: Mar 7, 2014
 *      Author: Varad Pathak, Arul Samuel
 */

#ifndef LDA_H_
#define LDA_H_

#include <omp.h>
#include <map>
#include <string>
#include <vector>
#include "Term.h"
#include "SCVB0.h"

using namespace std;

class LDA {
public:
	int numOfDoc;
	int numOfTerms;
	int numOfWordsInCorpus;
	int iterations;
	int numOfTopics;
	string fileName;
	vector<Term*>* termVector;

	LDA(string fileName, int iter, int topics);
	int main(int argv, char *argc[]);
	SCVB0* parseDataFile(int nProcessors);
	void printResults(SCVB0* scvb0);
	void executeSCVB0(SCVB0* scvb0);
	void normalize(SCVB0* scvb0);
	double normalizeAndPerplexity(SCVB0* scvb0);
	virtual ~LDA();
};

#endif /* LDA_H_ */
