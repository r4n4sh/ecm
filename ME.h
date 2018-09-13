#ifndef __ME_H_INCLUDED__
#define __ME_H_INCLUDED__
#include <list>
using namespace std;

class Bucket {
	unsigned int size;
	unsigned int timestamp;
public:
	Bucket();
	Bucket(unsigned int timestamp);
	Bucket(unsigned int size, unsigned int timestamp);
	unsigned int getSize();
	unsigned int getTimestamp();
	void setTimestamp(unsigned int tstmp);
	void decrementSize(unsigned int s);
};

class Info {
public:
	unsigned int size;
	list<Bucket>* buckets;
public:
	Info();
	Info(Bucket b);
	unsigned int getSize();
	void setSize(unsigned int s);
	void incrementSize();
	void decrementSize();
	void decrement2Size();
	list<Bucket>* getBuckets();
	//~Info();
};

class ME_Alg {
public:
	list<Info>* histogram;
	unsigned int N;
	unsigned long long time;
	unsigned int total;
	unsigned int k;
	unsigned int k2; // k/2
	unsigned int B;
	unsigned int* buff;
	unsigned int index;
	unsigned int clearTime;
public:
	ME_Alg(unsigned int W, double e);
	ME_Alg();
	void init(unsigned int W, double e);
	void insert(int v);
	void insert_time(int t);
	unsigned int query();
	unsigned int query_time(int w);
	unsigned int slow_query(int w);
	void printHistogram();
	int* l_CanonicalRepresentation(int S);
	int usedSpace();
	//~ME_Alg();
public:
	void histogramInsert(Bucket b);
};

#endif // __ME_H_INCLUDED__
