/*
 * SCVB0.h
 *
 *  Created on: Mar 7, 2014
 *      Author: Varad Pathak, Arul Samuel
 */

#ifndef SCVB0_H_
#define SCVB0_H_

#include <vector>
#include "Document.h"
#include "MiniBatch.h"
#include <iostream>

using namespace std;

class SCVB0 {
public:
	int iterations;
	int K; //Number of Topics
	int W; //no of terms in vocab
	int D; //Total no of docs in corpus
	int C; //Total no of words in corpus
	int numOfBurnInPasses;

	double **nPhi;
	double **nTheta;
	double *nz;

	double alpha;
	double eta;

	int s;
	int tau;
	double kappa;

	double rhoPhi;
	double rhoTheta;

	int rhoPhi_t;
	int rhoTheta_t;

	vector<MiniBatch*> *miniBatches;

	SCVB0(int iter, int numberOfTopics, int vocabSize, int numOfDocs,
			int corpusSize);
	virtual ~SCVB0();
	void run(MiniBatch* miniBatch);
};

#endif /* SCVB0_H_ */
