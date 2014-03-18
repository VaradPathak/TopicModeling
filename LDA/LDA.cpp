/*
 * LDA.cpp
 *
 *  Created on: Mar 7, 2014
 *      Author: Varad Pathak, Arul Samuel
 */

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "LDA.h"
#include "Document.h"
#include "SCVB0.h"
#include "MiniBatch.h"
#include "Term.h"
#include <omp.h>
#include <fstream>
#include <algorithm>
#include <math.h>
#define SUBMISSION true

int currTopic = 0;
struct myclass {
	bool operator()(Term* i, Term* j) {
		return ((*i->prob)[currTopic] > (*j->prob)[currTopic]);
	}
} myobject;

LDA::LDA(string file, int iter, int topics) {
	numOfDoc = 0;
	numOfTerms = 0;
	numOfWordsInCorpus = 0;
	fileName = file;
	iterations = iter;
	numOfTopics = topics;

	termVector = new vector<Term*>();
}

LDA::~LDA() {
	termVector->clear();
}

/* REQUIRED WHILE CALCULATING PERPLEXITY */
/*
double LDA::normalizeAndPerplexity(SCVB0* scvb0) {
	double sumProb = 0.0, perplexity = 0.0;
	for (int k = 0; k < scvb0->K; k++) {
		double k_total = 0;
		for (std::vector<Term*>::iterator it = termVector->begin();
				it != termVector->end(); it++) {
			Term* term = *it;
			k_total += scvb0->nPhi[term->wordId][k];
		}
		for (std::vector<Term*>::iterator it = termVector->begin();
				it != termVector->end(); it++) {
			Term* term = *it;
			double temp = scvb0->nPhi[term->wordId][k] / k_total;
			scvb0->nPhi[term->wordId][k] = temp;
			term->prob->push_back(temp);
		}
	}
	int d = 1;
	#pragma omp parallel for shared(d)
	for (d = 1; d < scvb0->D + 1; ++d) {

		double k_total = 0;
		for (int k = 0; k < scvb0->K; k++) {
			k_total += scvb0->nTheta[d][k];
		}

		for (int k = 0; k < scvb0->K; k++) {
			double temp = scvb0->nTheta[d][k] / k_total;
			scvb0->nTheta[d][k] = temp;
			for (std::vector<Term*>::iterator it = termVector->begin();
					it != termVector->end(); it++) {
				Term* term = *it;
				sumProb += log(temp * ((*term->prob)[k]));
			}
		}
	}
	cout << "sumProb: " << sumProb << " C: " << scvb0->C << endl;
	perplexity = exp(-sumProb / scvb0->C);
	return perplexity;
}
*/
void LDA::normalize(SCVB0* scvb0) {
	for (int k = 0; k < scvb0->K; k++) {
		double k_total = 0;
		for (std::vector<Term*>::iterator it = termVector->begin();
				it != termVector->end(); it++) {
			Term* term = *it;
			k_total += scvb0->nPhi[term->wordId][k];
		}
		for (std::vector<Term*>::iterator it = termVector->begin();
				it != termVector->end(); it++) {
			Term* term = *it;
			double temp = scvb0->nPhi[term->wordId][k] / k_total;
			scvb0->nPhi[term->wordId][k] = temp;
			term->prob->push_back(temp);
		}
	}
	int d = 1;
	#pragma omp parallel for shared(d)
	for (d = 1; d < scvb0->D + 1; ++d) {

		double k_total = 0;
		for (int k = 0; k < scvb0->K; k++) {
			k_total += scvb0->nTheta[d][k];
		}

		for (int k = 0; k < scvb0->K; k++) {
			double temp = scvb0->nTheta[d][k] / k_total;
			scvb0->nTheta[d][k] = temp;
		}
	}
}

void LDA::printResults(SCVB0* scvb0) {
	cout << "Writing results to file" << endl;
	int p = 0;
	#pragma omp parallel for shared(p)
	for (p = 0; p < 2; p++) {
		if (p == 0) {
			ofstream doctopicFile;
			doctopicFile.open("doctopic.txt");
			// Each document with its topic allocation
			for (int d = 1; d < scvb0->D + 1; ++d) {
				for (int k = 0; k < scvb0->K - 1; k++) {
					doctopicFile << scvb0->nTheta[d][k] << ",";
				}
				doctopicFile << scvb0->nTheta[d][scvb0->K - 1] << endl;
			}
			doctopicFile.close();
		} else {
			#if SUBMISSION
			ofstream topicFile;
			topicFile.open("topic.txt");
			int counter = 1;
			for (int k = 0; k < numOfTopics; k++) {
				currTopic = k;
				std::sort(termVector->begin(), termVector->end(), myobject);
				counter = 1;
				for (std::vector<Term*>::iterator it = termVector->begin();
						it != termVector->end(); it++) {
					Term* term = *it;
					if (counter < 100) {
						topicFile << term->wordId << ":" << (*term->prob)[k]
						<< ",";
					} else if (counter == 100) {
						topicFile << term->wordId << ":" << (*term->prob)[k];
					} else {
						break;
					}
					counter++;
				}
				topicFile << endl;
			}
			topicFile.close();
			#else
			ofstream topicFile;
			topicFile.open("topic.txt");
			int counter = 1;
			for (int k = 0; k < numOfTopics; k++) {
				currTopic = k;
				std::sort(termVector->begin(), termVector->end(), myobject);
				topicFile << "Topic " << k + 1 << " : " << endl;
				counter = 1;
				for (std::vector<Term*>::iterator it = termVector->begin();
						it != termVector->end(); it++) {
					if (counter < 21) {
						Term* term = *it;
						topicFile << term->word << " " << (*term->prob)[k]
								<< endl;
					} else {
						break;
					}
					counter++;
				}
				topicFile << endl;
			}
			topicFile.close();
			#endif
		}
	}
}

SCVB0 * LDA::parseDataFile(int nProcessors) {
	ifstream inputfile(fileName.c_str());
	inputfile >> numOfDoc;
	inputfile >> numOfTerms;
	inputfile >> numOfWordsInCorpus;
	cout << numOfDoc << " " << numOfTerms << " " << numOfWordsInCorpus << endl;
	cout << "--------------------" << endl;

	inputfile.close();

	SCVB0 *scvb0 = new SCVB0(iterations, numOfTopics, numOfTerms, numOfDoc,
			numOfWordsInCorpus);
	int megaBatchSize = (int) numOfDoc / nProcessors;

	int miniBatchSize = 100;
	int megaBatchStartDoc = 1;
	int megaBatchId = 0;
	int totalWork = (nProcessors * (nProcessors + 1) / 2);
	vector<int>* batchSizes = new vector<int>();
	for (megaBatchId = 0; megaBatchId < nProcessors; megaBatchId++) {
		batchSizes->push_back(
				(nProcessors - megaBatchId) * numOfDoc / totalWork);
	}
	int BatchCounter = nProcessors;

#pragma omp parallel for shared(BatchCounter, scvb0, batchSizes)
	for (megaBatchStartDoc = 1; megaBatchStartDoc <= numOfDoc;
			megaBatchStartDoc += (*batchSizes)[nProcessors - BatchCounter]) {
		BatchCounter--;
		int docId, wordId, freq;
		ifstream infile(fileName.c_str());
		//Ignoring first three lines of the file as they are already read
		int x, y, z;
		infile >> x;
		infile >> y;
		infile >> z;
		int eof = 0;
		if (!(infile >> docId >> wordId >> freq)) {
			eof = 1;
		}
		int megaBatchEndDoc = (megaBatchStartDoc + megaBatchSize);
		while (docId != megaBatchStartDoc) {
			infile >> docId >> wordId >> freq;
		}
		for (int a = megaBatchStartDoc; a <= megaBatchEndDoc; a +=
				miniBatchSize) {
			int docId, wordId, freq;
			int miniBatchStartDoc = a;
			int miniBatchEndDoc = miniBatchStartDoc + miniBatchSize;
			MiniBatch* miniBatch = new MiniBatch();
			miniBatch->M = 0;
			while (!eof && (docId < miniBatchEndDoc)) {
				map<int, int>* termMap = new map<int, int>();
				int oldDocId = docId;
				int Cj = 0;
				while (docId == oldDocId) {
					(*termMap)[wordId] = freq;
					Cj += freq;
					if (!(infile >> docId >> wordId >> freq)) {
						eof = 1;
						break;
					}
				}
				Document* newDoc = new Document(oldDocId, *termMap);
				newDoc->Cj = Cj;
				miniBatch->M += Cj;
				miniBatch->docVector->push_back(*newDoc);
			}
			scvb0->miniBatches->push_back(miniBatch);
		}
		infile.close();
	}
	//Reading the Vocabulary file
	ifstream myVocabFile;
	myVocabFile.open("vocab.kos.txt");
	int wordId = 1;
	string word;
	int eof = 0;
	while (!eof) {
		if (!(myVocabFile >> word)) {
			eof = 1;
			break;
		}
		Term* term = new Term(wordId, word);
		termVector->push_back(term);
		wordId++;
	}
	return scvb0;
}
LDA *parseCommandLine(int argv, char *argc[]) {
	argv--, argc++;
	if (argv < 3) {
		cerr<< "Malformed command. Correct format: ./fastLDA docword.txt iterations NumOfTopics"<< endl;
		exit(1);
	}
	LDA *lda = new LDA(string(argc[0]), atoi(argc[1]), atoi(argc[2]));
	return lda;
}

int main(int argv, char *argc[]) {

	//Setting the number of threads
	int nProcessors = omp_get_max_threads();
	omp_set_num_threads(nProcessors);
	
	double tStart = omp_get_wtime();
	LDA *lda = parseCommandLine(argv, argc);
	double tStart1 = omp_get_wtime();
	SCVB0 *scvb0 = lda->parseDataFile(nProcessors);
	double tEnd1 = omp_get_wtime();
	printf("Time taken to read files: %.2fs\n", (double) (tEnd1 - tStart1));

	for (int itr = 0; itr < lda->iterations; ++itr) {
		cout << "Iteration: " << itr + 1 << endl;
		int m = 0;
		#pragma omp parallel for shared(m)
		for (m = 0; m < (int) scvb0->miniBatches->size(); m++) {
			scvb0->rhoPhi_t = 1;
			scvb0->rhoTheta_t = 1;

			scvb0->run((*scvb0->miniBatches)[m]);
		}

		lda->normalize(scvb0);
		/*Call the normalizeAndPerplexity function for calculating the perplexity and comment the normalize function 
		double perplexity = lda->normalizeAndPerplexity(scvb0);
		cout << "Perplexity after iteration: " << itr + 1 << " is "
				<< perplexity << endl;
		*/
	}

	lda->printResults(scvb0);

	delete scvb0;
	delete lda;
	double tEnd = omp_get_wtime();
	cout << "Done " << endl;
	printf("Time taken: %.2fs\n", (double) (tEnd - tStart));
	return 0;
}
