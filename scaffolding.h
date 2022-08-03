#ifndef SCAFFOLDING_H_INCLUDED 
#define SCAFFOLDING_H_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "contig.h"
#include "aligningFromBam.h"
#include "scaffoldgraph.h"
#include "lp/lp_lib.h"

using namespace std;

typedef struct ContigSequence{
    int index;
    int orientation;
    int gapDistance;
    ContigSequence * next;
}ContigSequence;

typedef struct ScaffoldSet{
    int length;
    ContigSequence * contigSequence;
	int * contigIndex;
	int contigNum;
	int sequenceCount;
    ScaffoldSet * next;
}ScaffoldSet;

typedef struct ScaffoldSetHead{
    ScaffoldSet * scaffoldSet;
}ScaffoldSetHead;

long int RemoveOrientationContradictionInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, long int contigCount, bool * contigOrientation);

long int DetermineContigOrderInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, bool * contigOrientation, long int * contigPosition);

void OptimizeScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead);

ScaffoldSetHead * GetScaffoldingResult(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead);

void GetSingleScaffold(ScaffoldSetHead * scaffoldSetHead, ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, long int startNodeIndex, bool * visited);

ContigSequence *  DeterminNextNode(ScaffoldSetHead * scaffoldSetHead, ScaffoldGraphHead * scaffoldGraphHead, ContigSequence * tempContigSequence, bool * visited, bool & tail);

long int GetGapDistance(double weight);

long int SelectNodeFromContigSet(ContigSetHead * contigSetHead, bool * visited);

void OutputOneContig(FILE * fp, char * contig, char * tempLine, long int length);

void OutPutScaffoldTag(ScaffoldSet * scaffoldSet, char * result);

void OutPutScaffoldSet(ScaffoldSet * scaffoldSet, ContigSetHead * contigSetHead, char * result);


#endif
