/*
 * MiniBatch.h
 *
 *  Created on: Mar 13, 2014
 *      Author: vspathak
 */

#ifndef MINIBATCH_H_
#define MINIBATCH_H_

#include <vector>
#include "Document.h"

using namespace std;

class MiniBatch {
public:
	int M;
	vector<Document> *docVector;
	MiniBatch();
	MiniBatch(int m, vector<Document> *documentVector);
	virtual ~MiniBatch();
};

#endif /* MINIBATCH_H_ */
