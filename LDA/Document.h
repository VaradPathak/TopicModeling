/*
 * Document.h
 *
 *  Created on: Mar 7, 2014
 *      Author: Varad Pathak, Arul Samuel
 */

#ifndef DOCUMENT_H_
#define DOCUMENT_H_

#include <map>
#include <string>

using namespace std;

//typedef std::map<int,int> intMap;

class Document {
public:
	int docId;
	int Cj;
	map<int, int> termDict;
	Document();
	Document(int id, map<int, int> termMap);
	virtual ~Document();
};

#endif /* DOCUMENT_H_ */
