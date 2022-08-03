#ifndef aligningFromBam_H_INCLUDED 
#define aligningFromBam_H_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

#include "contig.h"

using namespace std;
using namespace BamTools;


typedef struct AligningResult{
	int leftContigIndex;//左边contig的位置
	int rightContigIndex;//右边contig的位置
	long int leftPosition;//左边read的的位置
	long int rightPosition;//右边read的的位置
	bool leftOrientation;//左边比对的方向
	bool rightOrientation;//右边比对的方向
}AligningResult;
//定义了一个结构体

typedef struct AligningResultHead{
	AligningResult * aligningResult;
	int aligningCount;
}AligningResultHead;

AligningResultHead * GetAligningResultHead(char * leftBamFile, char * rightBamFile, char * alignmentResultFile,char * sortalignmentResultFile, ContigSetHead * contigSetHead, bool aligning);

void sortAligningResult(AligningResult * a, long int left, long int right);

long int MIN(long int a, long int b);

long int MAX(long int a, long int b);

void OutputAligningResultHead(AligningResultHead * aligningResultHead, ContigSetHead * contigSetHead);

void OptimizeAligningResultHead(AligningResultHead * aligningResultHead);

#endif
