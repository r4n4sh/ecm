/********************************************************************
Count-Min Sketches

G. Cormode 2003,2004

Updated: 2004-06 Added a floating point sketch and support for
                 inner product point estimation
Initial version: 2003-12

This work is licensed under the Creative Commons
Attribution-NonCommercial License. To view a copy of this license,
visit http://creativecommons.org/licenses/by-nc/1.0/ or send a letter
to Creative Commons, 559 Nathan Abbott Way, Stanford, California
94305, USA.
*********************************************************************/

#include <stdlib.h>
#include <time.h>
#include <sys/timeb.h>
#include <cstring>
#include "prng.h"
#include "countmin.h"

#ifndef CLK_PER_SEC
#ifdef CLOCKS_PER_SEC
#define CLK_PER_SEC CLOCKS_PER_SEC
#endif
#endif

//#define TEST_ERROR_MEMORY
//#define TEST_UPDATE 1
//#define TEST_QUERY 1
#define TEST_QUERY_INTERVALS

#define min(x,y)	((x) < (y) ? (x) : (y))
#define max(x,y)	((x) > (y) ? (x) : (y))

double eps;	               /* 1+epsilon = approximation factor */
double delta = 0.01;               /* probability of failure */

/************************************************************************/
/* Routines to support Count-Min sketches                               */
/************************************************************************/

/*
int main() {
	printf("TEST!\n");

	struct timeb begintb, endtb;
	clock_t begint, endt;
	double time;

	//size_t stNumberOfPackets = 10000000;
	size_t stNumberOfPackets = 100000;
	size_t stRuns = 20;
	double dPhi = 0.001;
	uint32_t u32Depth = 4;
	uint32_t u32Granularity = 8;

	double dSkew = 1.0;

	uint32_t u32Width = 2.0 / dPhi;

	prng_type * prng;
	prng = prng_Init(44545, 2);
	int64_t a = (int64_t)(prng_int(prng) % MOD);
	int64_t b = (int64_t)(prng_int(prng) % MOD);
	prng_Destroy(prng);

	uint32_t u32DomainSize = 1048575;
	std::vector<uint32_t> exact(u32DomainSize + 1, 0);

	CMH_type* cmh = CMH_Init(u32Width, u32Depth, 32, u32Granularity);
	std::vector<uint32_t> data;

	Tools::Random r = Tools::Random(0xF4A54B);
	Tools::PRGZipf zipf = Tools::PRGZipf(0, u32DomainSize, dSkew, &r);
	size_t stCount = 0;

	for (int i = 0; i < stNumberOfPackets; ++i)
	{
		++stCount;
		//if (stCount % 500000 == 0)
		//	std::cerr << stCount << std::endl;

		uint32_t v = zipf.nextLong();

		uint32_t value = hash31(a, b, v) & u32DomainSize;
		data.push_back(value);
	}

	size_t stRunSize = data.size() / stRuns;
	size_t stStreamPos = 0;
	uint64_t nsecs;

	for (size_t run = 1; run <= stRuns; ++run)
	{

		begint = clock();
		ftime(&begintb);
		for (size_t i = stStreamPos; i < stStreamPos + stRunSize; ++i)
		{
			if (data[i] > 0)
				CMH_Update(cmh, data[i], 1);
		}
		endt = clock();
		ftime(&endtb);
		time = ((double)(endt - begint)) / CLK_PER_SEC;
	}
	printf("time: %d\n", time);

	cin.get();
	return 0;
}*/


int main(int argc, char * argv[]) {
	int counters = 100;
	int threshold = 1000;
	int n = 100000;
	int k = 1;
	double time;
	int i;
	int w, x, y, z;
	clock_t begint, endt;
	struct timeb begintb, endtb;
	unsigned * weights;
	unsigned * intervals;
	unsigned long * data;
	FILE * fp = NULL;
	int M = 1;
	float gamma = 4;
	double epsilon = 0.01;
	unsigned int interval_size;
	int window_size = 1600;
	int interval_1;
	int interval_2;
	double dPhi = 0.001;
	uint32_t u32Depth = 4;
	uint32_t u32Width = 2.0 / dPhi;
	int percentage = 1;

	for (int i = 1; i < argc; ++i)
	{
		if (strcmp(argv[i], "-np") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cout << "Missing number of packets." << std::endl;
				return -1;
			}
			n = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-c") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cout << "Missing epsilon" << std::endl;
				return -1;
			}
			counters = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-k") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cout << "Missing k" << std::endl;
				return -1;
			}
			k = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-t") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cout << "Missing interval size percentage (of window size)" << std::endl;
				return -1;
			}
			percentage = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-f") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cout << "Missing trace file." << std::endl;
				return -1;
			}
			fp = fopen(argv[4], "w");
		}

		else if (strcmp(argv[i], "-gamma") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cout << "Missing gamma." << std::endl;
				return -1;
			}
			gamma = atof(argv[i]);
		}
		else if (strcmp(argv[i], "-M") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cout << "Missing M." << std::endl;
				return -1;
			}
			M = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-w") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cout << "Missing window size." << std::endl;
				return -1;
			}

			window_size = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-i") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cout << "Missing interval1." << std::endl;
				return -1;
			}

			interval_1 = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-j") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cout << "Missing interval2." << std::endl;
				return -1;
			}

			interval_2 = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-d") == 0)
		{
			i++;
			if (i >= argc)
			{
				std::cout << "Missing delta" << std::endl;
				return -1;
			}

			delta = atoi(argv[i]);
		}

		else
		{
			cout << "Unknown parameter" << argv[i] << endl;
			return -1;
		}
	}

	/*if (n / counters >= threshold) {
		printf("Unacceptable parameters: eps*n >= theshold\n");
		return 0;
	}*/

	int interval_range = min(window_size, n);
	float size_precentage = percentage/100.0; // percenatge%
	interval_size = ceil(size_precentage * window_size); // 10% of window_size

	epsilon = (double)1 / (double)counters;
	data = (unsigned long *)malloc(sizeof(unsigned long) * n);
	weights = (unsigned *)malloc(sizeof(unsigned) * n);
#ifdef TEST_ERROR_MEMORY
	unsigned long* window = new unsigned long[window_size];
#endif
#if defined(TEST_QUERY) || defined(TEST_ERROR_MEMORY) || defined(TEST_QUERY_INTERVALS)
	int interval_arr_size = ceil(n / 1000);

	intervals = (unsigned *)malloc(sizeof(unsigned) * interval_arr_size);

#endif

	u32Depth = ceil(1/delta);
	u32Width = ceil(1/epsilon);

	CM_type* cm = CM_Init(u32Width, u32Depth, 32, window_size, epsilon);

	for (i = 0; i < n; i++) {
		fscanf(fp, "%d%d%d%d", &w, &x, &y, &z);
		data[i] = (unsigned long)256 * ((unsigned long)256 * ((unsigned long)256 * w + x) + y) + z;
		fscanf(fp, "%d%d%d%d", &w, &x, &y, &z);
		fscanf(fp, "%d", weights + i);
#if defined TEST_QUERY || defined (TEST_QUERY_INTERVALS)
		int interval_idx = i / 1000;
		intervals[interval_idx] = 1 + (int)rand() % (int)(0.99 * interval_range);
#endif

#if defined(TEST_ERROR_MEMORY) 
			window[i%window_size] = data[i];
#endif

	}

#ifdef TEST_UPDATE
	begint = clock();
	ftime(&begintb);

	for (i = 0; i < n; i++) {
		CM_Update(cm, data[i], i);
	}

	endt = clock();
	ftime(&endtb);

	time = ((double)(endt - begint)) / CLK_PER_SEC;
	//memory = maxmemusage();
	printf("./cm %d pairs took %lfs [%d counters %d window_size]\n", n, time, counters, window_size);

#endif

#ifdef TEST_QUERY
	/* Test Query times */

	for (i = 0; i < n; i++) {
		CM_Update(cm, data[i], i);
	}

	begint = clock();
	ftime(&begintb);
	for (i = 0; i < n; i++) {
		CM_IntervalQuery(cm, data[i], intervals[i / 1000], intervals[i / 1000] + interval_size);
	}

	endt = clock();
	ftime(&endtb);

	time = ((double)(endt - begint)) / CLK_PER_SEC;
	//memory = maxmemusage();

	printf("./cm %d pairs took %lfs [%d counters %d window_size]\n", n, time, counters, window_size);

#endif


#ifdef TEST_QUERY_INTERVALS
	for (i = 0; i < n; i++) {
		CM_Update(cm, data[i], i);
	}

	begint = clock();
	ftime(&begintb);
	for (i = 0; i < n; i++) {
		CM_IntervalQuery(cm, data[i], intervals[i / 1000], intervals[i / 1000] + interval_size);
	}

	endt = clock();
	ftime(&endtb);

	time = ((double)(endt - begint)) / CLK_PER_SEC;
	//memory = maxmemusage();

	printf("./cm %d pairs took %lfs [%d counters %d window_size %d interval_size]\n", n, time, counters, window_size, percentage);
#endif

#ifdef TEST_ERROR_MEMORY
	/* Test Query times */
	double estimated, curr_error = 0;
	double exact = 0;
	double emp_error = 0;

	for (i = 0; i < n; i++) {
		CM_Update(cm, data[i], i);
	}


    for (i = 0; i < n; i++)  {
		double exact = 0;
		int first = 0;
		int last = interval_size;

		for (int k = first; k< last; ++k) {
			if (window[k] == data[i])
				exact += 1;
		}

		estimated = CM_IntervalQuery(cm, data[i], first, last);

		//cout << "estimated: " << estimated << " exact: " << exact << endl;
		curr_error = exact - estimated;
		curr_error = pow(curr_error, 2);
		emp_error += curr_error;
    }

	emp_error = sqrt((emp_error/n));

	printf( "./cm %d pairs emp error: %lf [%d counters %d window_size]\n", n, emp_error, counters, window_size);

#endif

#ifdef TEST_ERROR_MEMORY
	delete[] window;
#endif

	CM_Destroy(cm);
	return 0;
}

CM_type * CM_Init(int width, int depth, int seed, int w, double epsilon)
{     // Initialize the sketch based on user-supplied size
  CM_type * cm;
  int j;
  prng_type * prng;

  cm=(CM_type *) malloc(sizeof(CM_type));
  prng=prng_Init(-abs(seed),2);
  // initialize the generator to pick the hash functions

  if (cm && prng)
  {
      cm->depth=depth;
      cm->width=width;
      cm->count=0;
			cm->histogram_size = depth * width;
	  //cm->histograms = (ME_Alg **)calloc(sizeof(ME_Alg *), cm->depth);
	  //cm->histograms[0] = (ME_Alg *)calloc(sizeof(ME_Alg), cm->depth*cm->width);
			cm->hasha=(unsigned int *)calloc(sizeof(unsigned int),cm->depth);
			cm->hashb=(unsigned int *)calloc(sizeof(unsigned int),cm->depth);
			cm->histograms = new ME_Alg*[cm->depth];
			if (!cm->histograms || !cm->hasha || !cm->hashb) {
				return NULL;
			}

			for (int i = 0; i < cm->depth; i++) {
				cm->hasha[i]=prng_int(prng) & MOD;
				cm->hashb[i]=prng_int(prng) & MOD;
				cm->histograms[i] = new ME_Alg[cm->width];
				if (!cm->histograms[i])
					return NULL;
				for (int j = 0; j < cm->width; j++) {
					cm->histograms[i][j].init(w, epsilon);
				}
			}
		}
		return cm;
}

CM_type * CM_Copy(CM_type * cmold)
{     // create a new sketch with the same parameters as an existing one
  CM_type * cm;
  int j;

  if (!cmold) return(NULL);
  cm=(CM_type *) malloc(sizeof(CM_type));
  if (cm)
    {
      cm->depth=cmold->depth;
      cm->width=cmold->width;
      cm->count=0;
	  cm->histograms = (ME_Alg **)calloc(sizeof(ME_Alg *), cm->depth);
	  cm->histograms[0] = (ME_Alg *)calloc(sizeof(ME_Alg), cm->depth*cm->width);
      cm->hasha=(unsigned int *)calloc(sizeof(unsigned int),cm->depth);
      cm->hashb=(unsigned int *)calloc(sizeof(unsigned int),cm->depth);
      if (cm->histograms && cm->histograms[0] && cm->hasha && cm->hashb)
	{
	  for (j=0;j<cm->depth;j++)
	    {
	      cm->hasha[j]=cmold->hasha[j];
	      cm->hashb[j]=cmold->hashb[j];
		  cm->histograms[j] = (ME_Alg *)cm->histograms[0] + (j*cm->width);
	    }
	}
      else cm=NULL;
    }
  return cm;
}

void CM_Destroy(CM_type * cm)
{     // get rid of a sketch and free up the space
  if (!cm) return;

  if (cm->histograms)
  {
	  if (cm->histograms[0]) free(cm->histograms[0]);
	  free(cm->histograms);
	  cm->histograms = NULL;
  }
  if (cm->hasha) free(cm->hasha); cm->hasha=NULL;
  if (cm->hashb) free(cm->hashb); cm->hashb=NULL;
  free(cm);  cm=NULL;
}

int CM_Size(CM_type * cm)
{ // return the size of the sketch in bytes
  int counts, hashes, admin, histograms;
  if (!cm) return 0;
  admin=sizeof(CM_type);
  counts=cm->width*cm->depth*sizeof(int);
  histograms = cm->width*cm->depth * sizeof(ME_Alg);
  hashes=cm->depth*2*sizeof(unsigned int);
  return(admin + hashes + counts + histograms);
}

void CM_Update(CM_type * cm, unsigned int item, int diff)
{
  int j;
  int hist_val;
  if (!cm) return;
  ++cm->count;
  for (j=0;j<cm->depth;j++) {
		cm->histograms[j][hash31(cm->hasha[j], cm->hashb[j], item) % cm->width].insert_time(diff);
	}

}

int CM_PointEst(CM_type * cm, unsigned int query)
{
  // return an estimate of the count of an item by taking the minimum
  int j, histans;

  if (!cm) return 0;
  histans = cm->histograms[0][hash31(cm->hasha[0], cm->hashb[0], query) % cm->width].query();
  for (j = 1; j<cm->depth; j++)
	histans = min(histans, cm->histograms[j][hash31(cm->hasha[j], cm->hashb[j], query) % cm->width].query());
  // this can be done more efficiently if the width is a power of two
  //return (ans);//** Rana **/
  return (histans);
}

int CM_PointEstN(CM_type * cm, unsigned int query, unsigned int n)
{
	int j, histans;
	if (!cm) return 0;
	histans = cm->histograms[0][hash31(cm->hasha[0], cm->hashb[0], query) % cm->width].slow_query(n);
	for (j = 1; j<cm->depth; j++)
		histans = min(histans, cm->histograms[j][hash31(cm->hasha[j], cm->hashb[j], query) % cm->width].slow_query(n));
	return (histans);
}

/* assume n2>n1*/
int CM_IntervalQuery(CM_type* cm, unsigned int query, unsigned int n1, unsigned int n2)
{
	return CM_PointEstN(cm, query, n2) - CM_PointEstN(cm, query, n1);
}

int CM_PointMed(CM_type * cm, unsigned int query)
{
  // return an estimate of the count by taking the median estimate
  // useful when counts can become negative
  // depth needs to be larger for this to work well
  int j, * ans, result=0, *histans;

  if (!cm) return 0;
  ans=(int *) calloc(1+cm->depth,sizeof(int));
  histans = (int *)calloc(1 + cm->depth, sizeof(int));
  for (j=0;j<cm->depth;j++)
	histans[j + 1] = cm->histograms[j][hash31(cm->hasha[j], cm->hashb[j], query) % cm->width].query();


  if (cm->depth==1)
    //result=ans[1]; //** Rana **//
	result=histans[1];
  else
    if (cm->depth==2)
      {
	//result=(ans[1]+ans[2])/2;
	//if (abs(ans[1]) < abs(ans[2])) //** Rana **//
	if (abs(histans[1]) < abs(histans[2]))

	  //result=ans[1]; else result=ans[2]; //** Rana **//
	  result = histans[1]; else result = histans[2];

	// special tweak for small depth sketches
      }
    else
      //result=(MedSelect(1+cm->depth/2,cm->depth,ans));  //** Rana **//
	  result = (MedSelect(1 + cm->depth / 2, cm->depth, histans));

  return result;
  // need to adjust for routine starting at 1
}

int CM_Compatible(CM_type * cm1, CM_type * cm2)
{ // test whether two sketches are comparable (have same parameters)
  int i;
  if (!cm1 || !cm2) return 0;
  if (cm1->width!=cm2->width) return 0;
  if (cm1->depth!=cm2->depth) return 0;
  for (i=0;i<cm1->depth;i++)
    {
      if (cm1->hasha[i]!=cm2->hasha[i]) return 0;
      if (cm1->hashb[i]!=cm2->hashb[i]) return 0;
    }
  return 1;
}

int64_t CM_InnerProd(CM_type * cm1, CM_type * cm2)
{ // Estimate the inner product of two vectors by comparing their sketches
  int i,j;
  int64_t result, tmp;

  result=0;
  if (CM_Compatible(cm1,cm2))
    {
      for (i=0;i<cm1->width;i++)
	//result+=cm1->counts[0][i]*cm2->counts[0][i]; //** Rana **//
	result += cm1->histograms[0][i].query() * cm2->histograms[0][i].query();

      for (j=1;j<cm1->depth;j++)
	{
	  tmp=0;
	  for (i=0;i<cm1->width;i++)
	    //tmp+=cm1->counts[j][i]*cm2->counts[j][i]; //** Rana **//
	    tmp+=cm1->histograms[j][i].query()*cm2->histograms[j][i].query();

	  result=min(tmp,result);
	}
    }
  return result;
}

int64_t CM_F2Est(CM_type * cm)
{ // Estimate the second frequency moment of the stream
  int i,j;
  int64_t result, tmp, *ans;

  if (!cm) return 0;
  ans=(int64_t *) calloc(1+cm->depth,sizeof(int64_t));

  for (j=0;j<cm->depth;j++)
    {
      result=0;
      for (i=0;i<cm->width;i+=2)
	{
	  //tmp=(cm->counts[j][i]-cm->counts[j][i+1]); //** Rana **//
	  tmp = (cm->histograms[j][i].query() - cm->histograms[j][i + 1].query());
	  result+=tmp*tmp;
	}
      ans[j+1]=result;
    }
  result=LLMedSelect((cm->depth+1)/2,cm->depth,ans);
  return result;
}

int CM_Residue(CM_type * cm, unsigned int * Q)
{
// CM_Residue computes the sum of everything left after the points
// from Q have been removed
// Q is a list of points, where Q[0] gives the length of the list

  char * bitmap;
  int i,j;
  int estimate=0, nextest;

  if (!cm) return 0;
  bitmap=(char *) calloc(cm->width,sizeof(char));
  for (j=0;j<cm->depth;j++)
    {
      nextest=0;
      for (i=0;i<cm->width;i++)
	bitmap[i]=0;
      for (i=1;i<Q[0];i++)
	bitmap[hash31(cm->hasha[j],cm->hashb[j],Q[i]) % cm->width]=1;
      for (i=0;i<cm->width;i++)
	if (bitmap[i] == 0) nextest += cm->histograms[j][i].query();

      estimate=max(estimate,nextest);
    }
  return(estimate);
}
