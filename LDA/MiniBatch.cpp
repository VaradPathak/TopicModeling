/*
 * MiniBatch.cpp
 *
 *  Created on: Mar 13, 2014
 *      Author: vspathak
 */

#include "MiniBatch.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

MiniBatch::MiniBatch(int m, vector<Document> *documentVector) {
	M = m;
	docVector = documentVector;
}

MiniBatch::MiniBatch() {
	M = 0;
	docVector = new vector<Document>();
}

MiniBatch::~MiniBatch() {
	docVector->clear();
}

