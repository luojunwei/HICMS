#ifndef SCAFFOLDGRAPH_H_INCLUDED 
#define SCAFFOLDGRAPH_H_INCLUDED 

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

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamReader.h"
#include "api/BamAlignment.h"

using namespace BamTools;

using namespace std;


typedef struct ScaffoldGraphEdge{
	int nodeIndex;
	long int * gapDistance;
	long int * pairedReadCount;
	long int allPairedReadCount;
	double weight;
	int orientationType;
	long int * * matrix;
	long int leftBinCount;
	long int rightBinCount;
}ScaffoldGraphEdge;


typedef struct ScaffoldGraph{
	ScaffoldGraphEdge * edge;
	long int edgeCount;
}ScaffoldGraph;


typedef struct ScaffoldGraphHead{
	ScaffoldGraph * scaffoldGraph;
	long int scaffoldGraphNodeCount;
	long int binSize;
}ScaffoldGraphHead;

ScaffoldGraphHead * GetScaffoldGraphHeadFromAligningResultHead(ContigSetHead * contigSetHead, AligningResultHead * aligningResultHead, long int contigLengthThreshold, long int minPairedReadCount, double ratio, char * scaffoldGraphFile ,char * RscaffoldGraphFile);

bool DetermineAddEdgeInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, AligningResultHead * aligningResultHead, long int startIndex, long int endIndex, double minWeight, long int minPairedReadCount, long int contigLengthThreshold, double ratio);

bool RecoverTwoEdgeInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead);

bool InsertEdgeInScaffold(ScaffoldGraphHead * scaffoldGraphHead, long int leftIndex, long int edgeIndex, long int rightIndex);

bool RemoveEdgeInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, long int leftIndex, long int rightIndex);

int GetLinkTypeBetweenTwoContigs(long int leftPosition, long int rightPosition, long int leftContigLength, long int rightContigLength);

void OutputScaffoldGraphHead(ScaffoldGraphHead * scaffoldGraphHead);

void OutputScaffoldGraphHead(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, long int * contigPosition);

void OutputscaffoldGraphHead(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, char* scaffoldGraphFile);

//void OutputScaffoldGraphHead(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, char* RscaffoldGraphFile);

int GetMaxPairedReadCountIndex(long int * pairedReadCount, int count, double minRatio);

void RemoveEdgeWithLowPairedReadCount(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, int minEdgeCount);

void RemoveEdgeWithLowRatio(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, double minRatio);

void dfs(ScaffoldGraphHead * scaffoldGraphHead, long int index, bool * visited, long int * contigPosition);

void DFS(ScaffoldGraphHead * scaffoldGraphHead, long int * contigPosition);

void SortEdgeWeight(ScaffoldGraphHead * scaffoldGraphHead);

void KeepEdgeWithMaxWeight(ScaffoldGraphHead * scaffoldGraphHead, long int count);

void RemoveEdgeWithLowWeight(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, double minWeight);

double Computeweight(long int ** matrix,long int rowCount, long int leftbincount, long int rightbincount, long int leftTailLength,long int rightTailLength , long int index, long int index1);

//double ComputeDI(long int ** matrix,long int rowCount, long int leftbincount, long int rightbincount, long int leftTailLength,long int rightTailLength);
#endif
