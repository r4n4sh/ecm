#include "ME.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

Bucket::Bucket() {
	size = timestamp = 1;
}

Bucket::Bucket(unsigned int timestamp) {
	size = 1;
	this->timestamp = timestamp;
}

Bucket::Bucket(unsigned int size, unsigned int timestamp) {
	this->size = size;
	this->timestamp = timestamp;
}
unsigned int Bucket::getSize() {
	return size;
}
unsigned int Bucket::getTimestamp() {
	return timestamp;
}

void Bucket::setTimestamp(unsigned int tstmp) {
	timestamp = tstmp;
}

void Bucket::decrementSize(unsigned int s) {
	size -= s;
}

Info::Info() {
	size = 0;
	buckets = new list<Bucket>();
}

Info::Info(Bucket b) {
	buckets = new list<Bucket>();
	buckets->push_back(b);
	size = 1;
}

unsigned int Info::getSize() {
	return size;
}

void Info::setSize(unsigned int s) {
	size = s;
}

void Info::incrementSize() {
	size++;
}

void Info::decrementSize() {
	size--;
}

void Info::decrement2Size() {
	size -= 2;
}

list<Bucket>* Info::getBuckets() {
	return buckets;
}

/*Info::~Info() {
	delete buckets;
}*/

int* ME_Alg::l_CanonicalRepresentation(int S) {
	int l = ceil(k/2.0);
	int j = floor(log2(S/((double) l) + 1));
	int Sprime = S - (pow(2, j) - 1) * l;
	int m = 0;
	if (Sprime >= pow(2,j))
		m = floor(Sprime/pow(2,j));
	int Sroof = Sprime - m * pow(2, j);
	int* can = new int[j + 1];
	for (int i = 0; i < j; i++)
		can[i] = l;
	can[j] = m;
	int i = 0;
	int n = Sroof;
	while(n != 0) {
		can[i] = n % 2 + l;
		n /= 2;
		i++;
	}
	return can;
}

ME_Alg::ME_Alg(unsigned int W, double e) {
	histogram = new list<Info>();
	this->N = W;
	time = 1;
	total = 0;
	k = ceil(1.0/double(e));
	k2 = ceil(double(k)/2.0) + 2;
	printf("k: %d k2:%d\n",k, k2 );
	B = log2(N);/***/
	buff = new unsigned int[B];
	for (int i = 0; i < B; i++)
		buff[i] = 0;
	index = 0;
	clearTime = time;
}

ME_Alg::ME_Alg() {
	histogram = new list<Info>();
	time = 1;
	total = 0;
	index = 0;
	clearTime = time;
}

void ME_Alg::init(unsigned int W, double e) {
	histogram = new list<Info>();
	this->N = W;
	time = 1;
	total = 0;
	k = ceil(1.0/double(e));
	k2 = ceil(double(k)/2.0) + 2;
	B = log2(N);/***/
	buff = new unsigned int[B];
	for (int i = 0; i < B; i++)
		buff[i] = 0;
	index = 0;
	clearTime = time;
}

void ME_Alg::insert(int v) {
	time++;

	if (index < B) {
		buff[index] = v;
		index++;
	}
	else {
		while (!histogram->empty() && time - histogram->back().getBuckets()->back().getTimestamp() >= N) {
			total -= histogram->back().getBuckets()->back().getSize();
			if (histogram->back().getSize() > 1) {
				histogram->back().getBuckets()->pop_back();
				histogram->back().decrementSize();
			}
			else
				histogram->pop_back();
		}


		int buffSum = 0;
		for (int i = 0; i < B; i++)
			buffSum += buff[i];
		int S = total + buffSum;
		int* K = l_CanonicalRepresentation(S);
		int j = floor(log2(S/((double) ceil(k/2.0)) + 1)) + 1;
		list<Info>* oldHistogram = histogram;
		histogram = new list<Info>();
		for (int i = 0; i < j; i++) {
			Info info = Info();
			for (int t = 0; t < K[i]; t++)
				info.getBuckets()->push_front(Bucket(pow(2,i), 0));
			info.setSize(K[i]);
			histogram->push_back(info);
		}

		index = B - 1;
		list<Info>::iterator it;
		list<Bucket>::iterator it2;
		for (it = histogram->begin(); it != histogram->end(); ++it) {
			it2 = (*it).getBuckets()->begin();
			for (; it2 != (*it).getBuckets()->end(); ++it2) {
				if (index < 0)
					break;
				if (buff[index] == 0)
					index--;
				if (index == 0)
					break;

				if ((*it2).getSize() <= buff[index]) {
					(*it2).setTimestamp(clearTime + index);
					buff[index] -= (*it2).getSize();
				}
				else {
					int diff = (*it2).getSize() - buff[index];
					(*it2).setTimestamp(clearTime + index);
					buff[index] = 0;
					index--;
					if (index >= 0)
						buff[index] -= diff;
				}
			}
			if (it2 != (*it).getBuckets()->end())
				break;
		}
		for (list<Info>::iterator it3 = oldHistogram->begin(); it3 != oldHistogram->end(); ++it3) {
			list<Bucket>::iterator it4;
			for (it4 = (*it3).getBuckets()->begin(); it4 != (*it3).getBuckets()->end();) {
				if (it4 == (*it3).getBuckets()->end())
					break;
				if ((*it4).getSize() == 0)
					++it4;
				if (it4 == (*it3).getBuckets()->end())
					break;

				if ((*it2).getSize() <= (*it4).getSize()) {
					(*it2).setTimestamp((*it4).getTimestamp());
					(*it4).decrementSize((*it2).getSize());
				}
				else {
					int diff = (*it2).getSize() - (*it4).getSize();
					(*it2).setTimestamp((*it4).getTimestamp());
					++it4;
					if (it4 != (*it3).getBuckets()->end())
						(*it4).decrementSize(diff);
				}

				 ++it2;
				 if (it2 == (*it).getBuckets()->end()) {
				 	++it;
				 	if (it == histogram->end())
				 		goto wrapup;
				 	it2 = (*it).getBuckets()->begin();
				 }
			}
			//if (it4 != (*it3).getBuckets()->end())
			//	break;
		}

		//delete oldHistogram;
wrapup:
		index = 0;
		total += buffSum;
		buff[index] = v;
		index++;
		clearTime = time - 1;
	}
}

void ME_Alg::insert_time(int t) {
	time= t;

	while (!histogram->empty() && time - histogram->back().getBuckets()->back().getTimestamp() >= N) {
		if (histogram->back().getSize() > 1) {
			histogram->back().getBuckets()->pop_back();
			histogram->back().decrementSize();
		}
		else
			histogram->pop_back();
	}

	list<Info>::iterator it = histogram->begin();
	int _size = 1;
	if (it->getSize() == 0) {
		Info info = Info();
		info.getBuckets()->push_front(Bucket(_size, time));
		info.incrementSize();
		histogram->push_back(info);

	} else {
		it->getBuckets()->push_front(Bucket(_size, time));
		it->incrementSize();
	}

	//printHistogram();
	it = histogram->begin();


	while (it->getSize() >= k2 && it != histogram->end()){
		//merge
		int tt = it->getBuckets()->back().getTimestamp();
		int olderStamp;
		if (it->getSize() >= 2){

			Bucket first = it->getBuckets()->back();
			Bucket last = it->getBuckets()->back();

			olderStamp = min(first.getTimestamp(), last.getTimestamp());
			it->getBuckets()->pop_back();
			it->getBuckets()->pop_back();
			it->decrement2Size();
		} else {
			Bucket first = it->getBuckets()->back();
			olderStamp = (first.getTimestamp());
			it->getBuckets()->pop_back();
			it->decrementSize();
		}

		++it;

		_size <<= 1;
		if (it != histogram->end()) {
			it->getBuckets()->push_front(Bucket(_size, olderStamp));
			it->incrementSize();
		} else {
			Info mergedinfo = Info();
			mergedinfo.getBuckets()->push_front(Bucket(_size, olderStamp));
			mergedinfo.incrementSize();
			histogram->push_back(mergedinfo);
		}
	}

	//printHistogram();

}

void ME_Alg::printHistogram()
{
	list<Info>::iterator it = histogram->begin();
	while (it != histogram->end()) {
		printf("info size: %d, \n",it->getSize());
		list<Bucket>::iterator it2 = (*it).getBuckets()->begin();
		while (it2 != (*it).getBuckets()->end()){
			printf("bucket timestamp: %d, \n",it2->getTimestamp());
			++it2;
		}
		++it;
	}
}

unsigned int ME_Alg::query() {
	int buffSum = 0;
	for (int i = 0; i < B; i++)
		buffSum += buff[i];
	int size = histogram->size();
	int last_size = pow(2, size - 1);
	if (histogram->back().getSize() == 0)
		last_size /= 2;
	return total + buffSum - last_size/2;
}

unsigned int ME_Alg::slow_query(int w) { //1's in the last w bits
	int buffSum = 0;
	int tmp = 0;
	list<Info>::iterator it = histogram->begin();
	//printHistogram();
	while (it != histogram->end()){
		list<Bucket>::iterator it2 = (*it).getBuckets()->begin();
//		printf("it2 backuts iterator timestamp: %d, time: %d, queried_time: %d\n", it2->getTimestamp(), time, w);
		while (it2->getTimestamp() >= time - w){
			buffSum += it2->getSize();
			++it2;
			++tmp;
		}
		//printf("number of relevant buckets is: %d, number of buckets:%d\n", tmp, it->getSize());
		//if (it2 == (*it).getBuckets()->end())//TODO
			//break;
		tmp = 0;
		++it;
	}

/*
	if (it == histogram->end())//TODO
	{
			printf("it is histogram end()\n" );
			buffSum += it->getBuckets()->back().getSize() >> 1;
	}
*/
	return buffSum;
}

unsigned int ME_Alg::query_time(int w) {
	int buffSum = 0;
	list<Info>::iterator it = histogram->begin();
	while (it->getBuckets()->back().getTimestamp() > time - w){
		buffSum += it->getSize() * it->getBuckets()->back().getSize();
		++it;
	}
	buffSum += (*it).getBuckets()->begin()->getSize() >> 1;
	list<Bucket>::iterator it2 = (*it).getBuckets()->begin();
	while (it2->getTimestamp() >= time - w){
		buffSum += it2->getSize();
		++it2;
	}
	return buffSum;
}

int ME_Alg::usedSpace() {
	return 0;
}

/*ME_Alg::~ME_Alg() {
	delete histogram;
}*/
