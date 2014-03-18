/*
 * Document.cpp
 *
 *  Created on: Mar 7, 2014
 *      Author: Varad Pathak, Arul Samuel
 */

#include "Document.h"

Document::Document() {
	docId = 0;
	Cj = 0;
}

Document::Document(int id, map<int, int> termMap) {
	docId = id;
	Cj = 0;
	termDict = termMap;
}

Document::~Document() {
	termDict.clear();
}

