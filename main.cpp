#include "ME.cpp"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
using namespace std;

void printBucket(Bucket b) {
	cout << "(" << b.getSize() << ", " << b.getTimestamp() << ") ";
}

void printInfo(Info info) {
	cout << info.getSize() << "; ";
	for (list<Bucket>::iterator it = info.getBuckets()->begin(); it != info.getBuckets()->end(); ++it) {
		printBucket(*it);
	}
	cout << endl;
}

void printHistogram(ME_Alg* alg) {
	for (list<Info>::iterator it = alg->histogram->begin(); it != alg->histogram->end(); ++it) {
		printInfo(*it);
	}
	cout << endl;
}

int correctness_test(int argc, char* argv[]) {
	if (argc < 3) {
        cerr << "Usage: " << argv[0] << " W e" << endl;
        return 1;
    }

	int W = atoi(argv[1]);
    double e = atof(argv[2]);

	ME_Alg* alg = new ME_Alg(W, e);

	ifstream myfile;
	myfile.open("my.dump");
	char output[100];
	clock_t time = clock();
	if (myfile.is_open()) {
		while (!myfile.eof()) {
		    myfile >> output;
		    int n = atoi(output);
	    	alg->insert(n);
		}
	}
	time = clock() - time;
	myfile.close();

    int as = alg->query();

    cout << "sum: " << as << endl;
	cout << "took " << time/(double)CLOCKS_PER_SEC << " seconds" << endl;

	return 0;
}

void printArray (int* arr, int n) {
	for (int i = n-1; i >=0; i--)
		cout << arr[i] << " ";
	cout << endl;
}
void printArray (unsigned int* arr, unsigned int n) {
	for (int i = 0; i < n; i++)
		cout << arr[i] << " ";
	cout << endl;
}

int main(int argc, char* argv[]) {
	return correctness_test(argc, argv);
}

