#ifndef BUILDSCAFFOLDGRAPH_CPP_INCLUDED 
#define BUILDSCAFFOLDGRAPH_CPP_INCLUDED 

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <exception>

#include "contig.h"
#include "aligningFromBam.h"
#include "scaffoldgraph.h"



using namespace std;

ScaffoldGraphHead * GetScaffoldGraphHeadFromAligningResultHead(ContigSetHead * contigSetHead, AligningResultHead * aligningResultHead, long int contigLengthThreshold, long int minPairedReadCount, double ratio , char * scaffoldGraphFile ,char * RscaffoldGraphFile){
	
	
	long int leftContigIndex = -1;
	long int rightContigIndex = -1;
	long int count = 0;

	long int typeIndex = -1;
	long int min = -1;
	int tempLinkType = -1;
	
	long int * edgeCount = (long int *)malloc(sizeof(long int)*contigSetHead->contigCount);


	for(long int i = 0; i < contigSetHead->contigCount; i++){
		
		edgeCount[i] = 0;

	}

	for(long int i = 0; i < aligningResultHead->aligningCount; i++){
		
		if(!((aligningResultHead->aligningResult[i].leftContigIndex == leftContigIndex && aligningResultHead->aligningResult[i].rightContigIndex == rightContigIndex) || (aligningResultHead->aligningResult[i].leftContigIndex == rightContigIndex && aligningResultHead->aligningResult[i].rightContigIndex == leftContigIndex))){
			leftContigIndex = aligningResultHead->aligningResult[i].leftContigIndex;
			rightContigIndex = aligningResultHead->aligningResult[i].rightContigIndex;
			edgeCount[leftContigIndex]++;
			edgeCount[rightContigIndex]++;
		}
	}
 
	ScaffoldGraphHead * scaffoldGraphHead = (ScaffoldGraphHead *)malloc(sizeof(ScaffoldGraphHead));
 
	scaffoldGraphHead->scaffoldGraphNodeCount = contigSetHead->contigCount;
 
	scaffoldGraphHead->scaffoldGraph = (ScaffoldGraph *)malloc(sizeof(ScaffoldGraph) * scaffoldGraphHead->scaffoldGraphNodeCount);
 
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		scaffoldGraphHead->scaffoldGraph[i].edgeCount = edgeCount[i];
		scaffoldGraphHead->scaffoldGraph[i].edge = (ScaffoldGraphEdge *)malloc(sizeof(ScaffoldGraphEdge)*scaffoldGraphHead->scaffoldGraph[i].edgeCount);
		for(long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++){
			scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex = -1;
			scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance = (long int *)malloc(sizeof(long int)*4);
			scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount = (long int *)malloc(sizeof(long int)*4);
			for(long int p = 0; p < 4; p++){
				scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance[p] = 0;
				scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[p] = 0;
			}
			scaffoldGraphHead->scaffoldGraph[i].edge[j].allPairedReadCount = 0;
			scaffoldGraphHead->scaffoldGraph[i].edge[j].weight = 0;
			scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType = -1;
		}
		scaffoldGraphHead->scaffoldGraph[i].edgeCount = 0;
	}
 
	leftContigIndex = -1;
	rightContigIndex = -1;
	long int startIndex = -1;
	long int endIndex = -1;
	
  double minWeight = 10;

	for(long int i = 0; i < aligningResultHead->aligningCount; i++){
		if(aligningResultHead->aligningResult[i].leftContigIndex == leftContigIndex && aligningResultHead->aligningResult[i].rightContigIndex == rightContigIndex){
			endIndex++;
		}else{
			DetermineAddEdgeInScaffoldGraph(scaffoldGraphHead, contigSetHead, aligningResultHead, startIndex, endIndex, minWeight, minPairedReadCount, contigLengthThreshold, ratio);
			startIndex = i;
			endIndex = i;
			leftContigIndex = aligningResultHead->aligningResult[i].leftContigIndex;
			rightContigIndex = aligningResultHead->aligningResult[i].rightContigIndex;
		}
	}
	
	//RemoveEdgeWithLowPairedReadCount(scaffoldGraphHead, contigSetHead, 5);
 
  //RemoveEdgeWithLowWeight(scaffoldGraphHead, contigSetHead, 10);
  
	//RemoveEdgeWithLowRatio(scaffoldGraphHead, contigSetHead, 3)
 
 //OutputScaffoldGraphHead(scaffoldGraphHead);
 
 // OutputscaffoldGraphHead(scaffoldGraphHead, contigSetHead, scaffoldGraphFile);

  cout<<"--Out put Scaffold Graph Head--"<<endl;
  
  //OutputScaffoldGraphHead(scaffoldGraphHead);
  
 // OutputscaffoldGraphHead(scaffoldGraphHead, contigSetHead, RscaffoldGraphFile);

	KeepEdgeWithMaxWeight(scaffoldGraphHead, 6);
 
	//DFS(scaffoldGraphHead, contigSetHead);
 
	cout<<"------------------"<<endl;
 
  OutputScaffoldGraphHead(scaffoldGraphHead);
	//RecoverTwoEdgeInScaffoldGraph(scaffoldGraphHead);
	
	return scaffoldGraphHead;
	
}


bool DetermineAddEdgeInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, AligningResultHead * aligningResultHead, long int startIndex, long int endIndex, double minWeight, long int minPairedReadCount, long int contigLengthThreshold, double ratio){
	
	if(endIndex - startIndex + 1 < minPairedReadCount){
		return false;
	}

 
	long int leftContigIndex = aligningResultHead->aligningResult[startIndex].leftContigIndex;
	long int rightContigIndex = aligningResultHead->aligningResult[startIndex].rightContigIndex;
 	long int leftContigLength = contigSetHead->contigSet[leftContigIndex].contigLength;
	long int rightContigLength = contigSetHead->contigSet[rightContigIndex].contigLength;
  long int LeftPosition = aligningResultHead->aligningResult[startIndex].leftPosition;
  long int RightPosition = aligningResultHead->aligningResult[startIndex].rightPosition;
  bool LeftOrientation = aligningResultHead->aligningResult[startIndex].leftOrientation;
  bool RightOrientation = aligningResultHead->aligningResult[startIndex].rightOrientation;
  
 	if(leftContigLength < contigLengthThreshold || rightContigLength < contigLengthThreshold){
		return false;
	}
 
	long int pairedReadCount = endIndex - startIndex + 1;
 
	long int leftTailLength = MIN(300000, leftContigLength/2);
	long int rightTailLength = MIN(300000, rightContigLength/2);
	long int leftPosition = 0;
	long int rightPosition = 0;
 
  double linkdensity = (double) pairedReadCount/(double)(leftContigLength+rightContigLength);
 //four Orientation
	for(long int i = startIndex; i <= endIndex; i++){
		leftPosition = aligningResultHead->aligningResult[i].leftPosition;
		rightPosition = aligningResultHead->aligningResult[i].rightPosition;
		//首首
		if(leftPosition < leftTailLength && rightPosition < rightTailLength){
			scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].pairedReadCount[0]++;
		}
   //尾尾
		if(leftPosition > leftContigLength - leftTailLength && rightPosition > rightContigLength - rightTailLength){
			scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].pairedReadCount[1]++;
		}
   //首尾
		if(leftPosition < leftTailLength && rightPosition > rightContigLength - rightTailLength){
			scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].pairedReadCount[2]++;
		}
   //尾首
		if(leftPosition > leftContigLength - leftTailLength && rightPosition < rightTailLength){
			scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].pairedReadCount[3]++;
		}
   
	}

	long int maxOrientationPairedReadCount = -1;
	long int maxOrientationIndex = -1;
	long int secondMaxOrientationPairedReadCount = -1;
	long int secondMaxOrientationIndex = -1;
	
	for(long int i = 0; i <= 3; i++){
		if(scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].pairedReadCount[i] > maxOrientationPairedReadCount)
    {
		maxOrientationPairedReadCount = scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].pairedReadCount[i];
			maxOrientationIndex = i;
		}
	}


	for(long int i = 0; i <= 3; i++){
		if(scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].pairedReadCount[i] > secondMaxOrientationPairedReadCount && i != maxOrientationIndex){
			secondMaxOrientationPairedReadCount = scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].pairedReadCount[i];
			secondMaxOrientationIndex = i;
		}
	}



	if(maxOrientationPairedReadCount < minPairedReadCount || (secondMaxOrientationPairedReadCount != 0 && double(maxOrientationPairedReadCount/secondMaxOrientationPairedReadCount) < ratio)){
		scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].pairedReadCount[0] = 0;
		//scaffoldGraphHead->scaffoldGraph[rightContigIndex].edge[scaffoldGraphHead->scaffoldGraph[rightContigIndex].edgeCount].pairedReadCount[0] = 0;
		scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].pairedReadCount[1] = 0;
		//scaffoldGraphHead->scaffoldGraph[rightContigIndex].edge[scaffoldGraphHead->scaffoldGraph[rightContigIndex].edgeCount].pairedReadCount[1] = 0;
		scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].pairedReadCount[2] = 0;
		//scaffoldGraphHead->scaffoldGraph[rightContigIndex].edge[scaffoldGraphHead->scaffoldGraph[rightContigIndex].edgeCount].pairedReadCount[2] = 0;
		scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].pairedReadCount[3] = 0;
		//scaffoldGraphHead->scaffoldGraph[rightContigIndex].edge[scaffoldGraphHead->scaffoldGraph[rightContigIndex].edgeCount].pairedReadCount[3] = 0;
		return false;
	}

 
 	scaffoldGraphHead->binSize =30000;
  
  int leftbincount;
  int rightbincount;

	if((leftTailLength%scaffoldGraphHead->binSize) == 0){
     scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].leftBinCount = leftTailLength/scaffoldGraphHead->binSize;	
	   }else{
		 scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].leftBinCount = leftTailLength/scaffoldGraphHead->binSize + 1; 
	   }
	   
 if((rightTailLength%scaffoldGraphHead->binSize) == 0){
      scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].rightBinCount = rightTailLength/scaffoldGraphHead->binSize;
	   }else{
		 scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].rightBinCount = rightTailLength/scaffoldGraphHead->binSize + 1;
	   }  
  
	 leftbincount=scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].leftBinCount;
	 rightbincount=scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].rightBinCount;

	long int rowCount =  leftbincount+rightbincount;

	long int ** matrix = (long int **)malloc(sizeof(long int *)*rowCount);
	for(long int i = 0; i < rowCount; i++){
		matrix[i] = (long int *)malloc(sizeof(long int)*rowCount);
		for(long int j = 0; j < rowCount; j++){
			matrix[i][j] = 0;
		}
	}

	//Matrix element fill
	for(long int i = startIndex; i <= endIndex; i++){
		leftPosition = aligningResultHead->aligningResult[i].leftPosition;
		rightPosition = aligningResultHead->aligningResult[i].rightPosition;
      
		if(leftPosition < leftTailLength && rightPosition < rightTailLength){
			if(maxOrientationIndex != 0){
				continue;
			}
			long int leftBinIndex = leftbincount - 1 - leftPosition/scaffoldGraphHead->binSize;
			long int rightBinIndex = rightPosition/scaffoldGraphHead->binSize;
			rightBinIndex = rightBinIndex + leftbincount;
			matrix[leftBinIndex][rightBinIndex]++;
			matrix[rightBinIndex][leftBinIndex]++;		
		}		
		if(leftPosition > leftContigLength - leftTailLength && rightPosition > rightContigLength - rightTailLength){
			if(maxOrientationIndex != 1){
				continue;
			}
			long int leftBinIndex = (leftTailLength - leftContigLength + leftPosition)/scaffoldGraphHead->binSize;
			long int rightBinIndex = rightbincount - 1 - (rightTailLength - rightContigLength + rightPosition)/scaffoldGraphHead->binSize;
			rightBinIndex = rightBinIndex + leftbincount;
			matrix[leftBinIndex][rightBinIndex]++;
			matrix[rightBinIndex][leftBinIndex]++;	
      }		
		if(leftPosition < leftTailLength && rightPosition > rightContigLength - rightTailLength){
			if(maxOrientationIndex != 2){
				continue;
			}
			long int leftBinIndex = leftbincount - 1 - leftPosition/scaffoldGraphHead->binSize;
			long int rightBinIndex = rightbincount - 1 - (rightTailLength - rightContigLength + rightPosition)/scaffoldGraphHead->binSize;
			rightBinIndex = rightBinIndex + leftbincount;
			matrix[leftBinIndex][rightBinIndex]++;
			matrix[rightBinIndex][leftBinIndex]++;			
		}
		if(leftPosition > leftContigLength - leftTailLength && rightPosition < rightTailLength){
			if(maxOrientationIndex != 3){
				continue;
			}
			long int leftBinIndex = (leftTailLength - leftContigLength + leftPosition)/scaffoldGraphHead->binSize;
			long int rightBinIndex = rightPosition/scaffoldGraphHead->binSize;
			rightBinIndex = rightBinIndex + leftbincount;
			matrix[leftBinIndex][rightBinIndex]++;
			matrix[rightBinIndex][leftBinIndex]++;	
		}
	}
 
	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].matrix = matrix; 

	long int leftbinIndex = leftbincount - 1;
	long int rightbinIndex = leftbincount;
		
  cout<<"infor: "<<leftContigIndex<<", "<<rightContigIndex<<", "<<leftContigLength<<", "<<rightContigLength<<", "<<LeftPosition<<", "<<RightPosition<<", "<<LeftOrientation<<", "<<RightOrientation<<endl;
  cout<<"tailLength: "<<leftTailLength<<","<<rightTailLength<<", "<<"pairedReadCount: "<<pairedReadCount<<endl;  
  cout<<"leftbincount: "<<leftbincount<<", "<<"rightbincount: "<<rightbincount<<", "<<"rowCount: "<<rowCount<<", "<<"linkdensity: "<<linkdensity<<endl;

//输出矩阵
 	for(long int i = 0; i < rowCount; i++){
		for(long int j = 0; j < rowCount; j++){
			cout<<matrix[i][j]<<",";
		}
		cout<<endl;
	}
 
  double weight = Computeweight(matrix, rowCount, leftbincount, rightbincount,  leftTailLength, rightTailLength , leftbinIndex, rightbinIndex);  
  //double weight = ComputeDI(matrix,rowCount, leftbincount, rightbincount, leftTailLength, rightTailLength);
  
  cout<<"Computeweight:"<<leftContigIndex<<","<<rightContigIndex<<",  "<<"weight: "<<weight<<endl; 
  cout<<" "<<endl;

	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].orientationType = maxOrientationIndex; 
	//scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].nodeIndex = leftContigIndex;
  scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].nodeIndex = rightContigIndex;
	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].weight = weight;
	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].allPairedReadCount = pairedReadCount;
	//scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].gapDistance = gapDistance;
 	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].matrix= matrix;
 	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].leftBinCount = leftbincount;
	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edge[scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount].rightBinCount = rightbincount;
 
	scaffoldGraphHead->scaffoldGraph[leftContigIndex].edgeCount++;

	return true;
	
}
/*
double Computeweight(long int ** matrix,long int rowCount, long int leftbincount, long int rightbincount, long int leftTailLength,long int rightTailLength,long int index, long int index1){
	double a = 0;
	double b = 0;
	 
	for(long int i = index + 1; i < rowCount; i++){
		b = b + matrix[index][i];
	}
	
	for(long int i = 0; i < index1; i++){
		a = a + matrix[i][index1];
	}
	
	double e = (a + b)/2;
 
	cout<<"a="<<a<<","<<"b="<<b<<","<<"e="<<e<<endl;

  if(a==0 && b==0){
  double weight =0;
  return weight;
  }else{
  double weight =(pow(a - e, 2) + pow(b -e , 2))/(leftbincount*rightbincount);
  return weight;
  }
}
*/
/*
double Computeweight(long int ** matrix,long int rowCount, long int leftbincount, long int rightbincount, long int leftTailLength,long int rightTailLength,long int index, long int index1){
	double a = 0;
	double b = 0;
  long int e = 0;
   
	for(long int i = index + 1; i < rowCount; i++){
		b = b + matrix[index][i];
	}
	
	for(long int i = 0; i < index1; i++){
		a = a + matrix[i][index1];
	}
	
 	for(long int i = 0; i < rowCount; i++){
		for(long int j = 0; j < rowCount; j++){
			e = e + matrix[i][j];
		}
	}
   
  long int sum= e/3;
 
 double weight = e;
 return weight;
}
*/

double Computeweight(long int ** matrix,long int rowCount, long int leftbincount, long int rightbincount, long int leftTailLength,long int rightTailLength,long int index, long int index1){
	
	double a = 0;
	double b = 0;
  double z = 0;
  int minWeight = 1;
  
 	for(long int i = 0; i < index1; i++){
		a = a + matrix[i][index1];
	}
  
	for(long int i = index + 1; i < rowCount; i++){
		b = b + matrix[index][i];
	}
 
 double x = 0;
  for(long int i = 0 ; i < rowCount ; i++){
		for(long int j = 0 ; j < rowCount; j++){
		   x = x + matrix[i][j] ;
		}
	}
  double sum = x/2;
	double e = sum/(leftbincount*rightbincount);
 
 
 for(long int i = 0 ; i < rowCount ; i++){
		for(long int j = 0 ; j < rowCount; j++){
		  if(matrix[i][j]!=0){
      z=z+1;
      }
		}
	}
  double w = z/(2*(leftbincount*rightbincount));
 	cout<<"a="<<a<<","<<"b="<<b<<","<<"w="<<w<<","<<"sum="<<sum<<","<<"e="<<e<<","<<"w="<<w<<endl;
  double weight = e*w ;
  return weight;
}

bool RecoverTwoEdgeInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead){
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		for(long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++){
			if(i < scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex){
				InsertEdgeInScaffold(scaffoldGraphHead, i, j, scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex);
			}
		}
	}
	return true;
}

bool InsertEdgeInScaffold(ScaffoldGraphHead * scaffoldGraphHead, long int leftIndex, long int edgeIndex, long int rightIndex){
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraph[leftIndex].edgeCount; i++){
		if(scaffoldGraphHead->scaffoldGraph[rightIndex].edge[i].nodeIndex == leftIndex){
			return false;
		}
	}
	long int edgeCount = scaffoldGraphHead->scaffoldGraph[rightIndex].edgeCount;
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].nodeIndex = leftIndex;
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].weight= scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].weight;
	if(scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].orientationType == 0 || scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].orientationType == 1){
		scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].orientationType = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].orientationType;
	}else if(scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].orientationType == 2){
		scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].orientationType = 3;
	}else{
		scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].orientationType = 2;
	}
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].pairedReadCount[0]= scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].pairedReadCount[0];
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].pairedReadCount[1]= scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].pairedReadCount[1];
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].pairedReadCount[2]= scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].pairedReadCount[3];
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].pairedReadCount[3]= scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].pairedReadCount[2];
	scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].allPairedReadCount = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].allPairedReadCount;
  scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].gapDistance = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].gapDistance;
  scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].matrix = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].matrix;
  scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].leftBinCount = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].leftBinCount;
  scaffoldGraphHead->scaffoldGraph[rightIndex].edge[edgeCount].rightBinCount = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[edgeIndex].rightBinCount;
	scaffoldGraphHead->scaffoldGraph[rightIndex].edgeCount++;
	return true;
	
}


bool RemoveEdgeInScaffoldGraph(ScaffoldGraphHead * scaffoldGraphHead, long int leftIndex, long int rightIndex){
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraph[leftIndex].edgeCount; i++){
		if(scaffoldGraphHead->scaffoldGraph[leftIndex].edge[i].nodeIndex == rightIndex){
			for(long int j = i; j < scaffoldGraphHead->scaffoldGraph[leftIndex].edgeCount - 1; j++){
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].nodeIndex= scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].nodeIndex;
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].weight= scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].weight;
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].orientationType= scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].orientationType;
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].pairedReadCount[0]= scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].pairedReadCount[0];
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].pairedReadCount[1]= scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].pairedReadCount[1];
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].pairedReadCount[2]= scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].pairedReadCount[2];
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].pairedReadCount[3]= scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].pairedReadCount[3];
				scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].allPairedReadCount = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].allPairedReadCount;
        scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].gapDistance = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].gapDistance;
        scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].matrix = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].matrix;
        scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].leftBinCount = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].leftBinCount;
        scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j].rightBinCount = scaffoldGraphHead->scaffoldGraph[leftIndex].edge[j + 1].rightBinCount;
				
			}
			scaffoldGraphHead->scaffoldGraph[leftIndex].edgeCount--;
			return true;
		}
	}
	return false;
}


void RemoveEdgeWithLowPairedReadCount(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, int minEdgeCount){
	long int count = 0;
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		count = 0;
		for(long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++){
			if(scaffoldGraphHead->scaffoldGraph[i].edge[j].allPairedReadCount < minEdgeCount){
				count++;
				continue;
			}else{
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].nodeIndex= scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex;
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount[0]= scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[0];
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount[1]= scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[1];
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount[2]= scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[2];
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount[3]= scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[3];
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].allPairedReadCount = scaffoldGraphHead->scaffoldGraph[i].edge[j].allPairedReadCount;
			}
			
		}
		scaffoldGraphHead->scaffoldGraph[i].edgeCount = scaffoldGraphHead->scaffoldGraph[i].edgeCount - count;
	}
	
}

void RemoveEdgeWithLowWeight(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, double minWeight){
	long int count = 0;
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		count = 0;
		for(long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++){
			if(scaffoldGraphHead->scaffoldGraph[i].edge[j].weight < minWeight){
				count++;
				continue;
			}else{
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].nodeIndex= scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex;
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].gapDistance = scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance;
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount[0]= scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[0];
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount[1]= scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[1];
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount[2]= scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[2];
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount[3]= scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[3];
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount = scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount;				
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].allPairedReadCount = scaffoldGraphHead->scaffoldGraph[i].edge[j].allPairedReadCount;
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].weight = scaffoldGraphHead->scaffoldGraph[i].edge[j].weight;
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].orientationType = scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType;
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].matrix = scaffoldGraphHead->scaffoldGraph[i].edge[j].matrix;
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].leftBinCount = scaffoldGraphHead->scaffoldGraph[i].edge[j].leftBinCount;
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].rightBinCount = scaffoldGraphHead->scaffoldGraph[i].edge[j].rightBinCount;
				
			}
			
		}
		scaffoldGraphHead->scaffoldGraph[i].edgeCount = scaffoldGraphHead->scaffoldGraph[i].edgeCount - count;
	}
	
}



void RemoveEdgeWithLowRatio(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, double minRatio){
	long int count = 0;
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		count = 0;
		for(long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++){
			if(GetMaxPairedReadCountIndex(scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount, 4, minRatio) < 0){
				count++;
				continue;
			}else{
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].nodeIndex= scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex;
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount[0]= scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[0];
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount[1]= scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[1];
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount[2]= scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[2];
				scaffoldGraphHead->scaffoldGraph[i].edge[j - count].pairedReadCount[3]= scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[3];
			}
			
		}
		scaffoldGraphHead->scaffoldGraph[i].edgeCount = scaffoldGraphHead->scaffoldGraph[i].edgeCount - count;
	}
	
}

void SortEdgeWeight(ScaffoldGraphHead * scaffoldGraphHead){
	long int temp = 0;
	double weight = 0;
  long int * gap = 0;
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		for(long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount - 1; j++){
			for(long int p = j + 1; p < scaffoldGraphHead->scaffoldGraph[i].edgeCount; p++){
				if(scaffoldGraphHead->scaffoldGraph[i].edge[j].weight < scaffoldGraphHead->scaffoldGraph[i].edge[p].weight){
					temp = scaffoldGraphHead->scaffoldGraph[i].edge[j].allPairedReadCount;
					scaffoldGraphHead->scaffoldGraph[i].edge[j].allPairedReadCount = scaffoldGraphHead->scaffoldGraph[i].edge[p].allPairedReadCount;
					scaffoldGraphHead->scaffoldGraph[i].edge[p].allPairedReadCount = temp;
					
					temp = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex;
					scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex = scaffoldGraphHead->scaffoldGraph[i].edge[p].nodeIndex;;
					scaffoldGraphHead->scaffoldGraph[i].edge[p].nodeIndex = temp;
					
					temp = scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType;
					scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType = scaffoldGraphHead->scaffoldGraph[i].edge[p].orientationType;;
					scaffoldGraphHead->scaffoldGraph[i].edge[p].orientationType = temp;
					
					temp = scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[0];
					scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[0] = scaffoldGraphHead->scaffoldGraph[i].edge[p].pairedReadCount[0];
					scaffoldGraphHead->scaffoldGraph[i].edge[p].pairedReadCount[0] = temp;
					
					temp = scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[1];
					scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[1] = scaffoldGraphHead->scaffoldGraph[i].edge[p].pairedReadCount[1];
					scaffoldGraphHead->scaffoldGraph[i].edge[p].pairedReadCount[1] = temp;
					
					temp = scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[2];
					scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[2] = scaffoldGraphHead->scaffoldGraph[i].edge[p].pairedReadCount[2];
					scaffoldGraphHead->scaffoldGraph[i].edge[p].pairedReadCount[2] = temp;
					
					temp = scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[3];
					scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[3] = scaffoldGraphHead->scaffoldGraph[i].edge[p].pairedReadCount[3];
					scaffoldGraphHead->scaffoldGraph[i].edge[p].pairedReadCount[3] = temp;

					weight = scaffoldGraphHead->scaffoldGraph[i].edge[j].weight;
					scaffoldGraphHead->scaffoldGraph[i].edge[j].weight = scaffoldGraphHead->scaffoldGraph[i].edge[p].weight;
					scaffoldGraphHead->scaffoldGraph[i].edge[p].weight = weight;

          gap = scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance;
					scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance = scaffoldGraphHead->scaffoldGraph[i].edge[p].gapDistance;
					scaffoldGraphHead->scaffoldGraph[i].edge[p].gapDistance = gap;			
                                                              
				}
			}
		}
	}
}

void KeepEdgeWithMaxWeight(ScaffoldGraphHead * scaffoldGraphHead, long int count){
	
	SortEdgeWeight(scaffoldGraphHead);
	
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		if(scaffoldGraphHead->scaffoldGraph[i].edgeCount > count){
			scaffoldGraphHead->scaffoldGraph[i].edgeCount = count;
		}
	}
	/*
	long int temp = 0;
	long int index = 0;
	bool real = false;
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		temp = 0;
		for(long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++){
			index = scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex;
			real = false;
			for(long int p = 0; p < scaffoldGraphHead->scaffoldGraph[index].edgeCount; p++){
				if(scaffoldGraphHead->scaffoldGraph[index].edge[p].nodeIndex == i){
					real = true;
					break;
				}
			}
			if(real != true){
				scaffoldGraphHead->scaffoldGraph[index].edge[scaffoldGraphHead->scaffoldGraph[index].edgeCount].nodeIndex = i;
				scaffoldGraphHead->scaffoldGraph[index].edge[scaffoldGraphHead->scaffoldGraph[index].edgeCount].weight = scaffoldGraphHead->scaffoldGraph[i].edge[j].weight;
				scaffoldGraphHead->scaffoldGraph[index].edge[scaffoldGraphHead->scaffoldGraph[index].edgeCount].allPairedReadCount = scaffoldGraphHead->scaffoldGraph[i].edge[j].allPairedReadCount;
				scaffoldGraphHead->scaffoldGraph[index].edge[scaffoldGraphHead->scaffoldGraph[index].edgeCount].pairedReadCount[0] = scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[0];
				scaffoldGraphHead->scaffoldGraph[index].edge[scaffoldGraphHead->scaffoldGraph[index].edgeCount].pairedReadCount[1] = scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[1];
				scaffoldGraphHead->scaffoldGraph[index].edge[scaffoldGraphHead->scaffoldGraph[index].edgeCount].pairedReadCount[2] = scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[2];
				scaffoldGraphHead->scaffoldGraph[index].edge[scaffoldGraphHead->scaffoldGraph[index].edgeCount].pairedReadCount[3] = scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[3];
				scaffoldGraphHead->scaffoldGraph[index].edgeCount++;
			}
		}
	}
	*/
}

int GetLinkTypeBetweenTwoContigs(long int leftPosition, long int rightPosition, long int leftContigLength, long int rightContigLength){
	int a = 0;
	int b = 0;
	long int minLength = 100000;
	
	long int c = 1;
	long int leftMin = MIN(leftContigLength / c, minLength);
	long int rightMin = MIN(rightContigLength / c, minLength);
	
	if(leftPosition < leftMin){
		a = 0;
	}else if(leftPosition > leftContigLength - leftMin){
		a = 1;
	}else{
		return -1;
	}
	
	if(rightPosition < rightMin){
		b = 0;
	}else if(rightPosition > rightContigLength - rightMin){
		b = 1;
	}else{
		return -1;
	}
	
	return a*2 + b;
	
}

int GetMaxPairedReadCountIndex(long int * pairedReadCount, int count, double minRatio){
	
	int index = -1;
	int index1 = -1;
	long int max = 0;
	for(long int i = 0; i < count; i++){
		if(pairedReadCount[i] > max){
			max = pairedReadCount[i];
			index = i;
		}
	}
	max = 0;
	for(long int i = 0; i < count; i++){
		if(pairedReadCount[i] > max && i != index){
			max = pairedReadCount[i];
			index1 = i;
		}
	}
	if(max == 0){
		return index;
	}
	double a = double(pairedReadCount[index])/double(pairedReadCount[index1]);
	
	
	if(a > minRatio){
		return index;
	}else{
		return -1;
	}
	
}


void OutputScaffoldGraphHead(ScaffoldGraphHead * scaffoldGraphHead){
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		if(scaffoldGraphHead->scaffoldGraph[i].edgeCount > 4){
			cout<<"------------------------------------"<<endl;
		}
		for(long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++){
			cout<<i<<"--"<<scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex<<"--"<<scaffoldGraphHead->scaffoldGraph[i].edge[j].weight<<"--"<<scaffoldGraphHead->scaffoldGraph[i].edge[j].allPairedReadCount<<"--type:"<<scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType<<endl;
			for(long int p = 0; p < 4; p++){
				cout<<"("<<scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance[p]<<",";
				cout<<scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[p]<<"),";
			}
			cout<<endl;
		}
	}
	
}



void OutputScaffoldGraphHead(ScaffoldGraphHead * scaffoldGraphHead, ContigSetHead * contigSetHead, long int * contigPosition){
  
  for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		for(long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++){
			cout<<i<<"--"<<scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex<<"--"<<scaffoldGraphHead->scaffoldGraph[i].edge[j].weight<<"--"<<scaffoldGraphHead->scaffoldGraph[i].edge[j].allPairedReadCount<<"--type:"<<scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType<<endl;
			cout<<"position:"<<contigPosition[i]<<"--"<<contigPosition[scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex]<<endl;
			cout<<"length:"<<contigSetHead->contigSet[i].contigLength<<"--"<<contigSetHead->contigSet[scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex].contigLength<<endl;
			for(long int p = 0; p < 4; p++){
				cout<<"("<<scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance[p]<<",";
				cout<<scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount[p]<<"),";
			}
			cout<<endl;
		}
	}	
}


void OutputscaffoldGraphHead(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, char* scaffoldGraphFile) {
	FILE* fp;
	if ((fp = fopen(scaffoldGraphFile, "w")) == NULL) {
		printf("%s, does not exist!", scaffoldGraphFile);
		exit(0);
	}
	for (long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		for (long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++) {
			fprintf(fp, "%ld--%d--%d--%lf--%ld--%ld\n", i, scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex, scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType, scaffoldGraphHead->scaffoldGraph[i].edge[j].weight,*scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount, *scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance);
		}
	}
	fclose(fp);
 }
/*
void OutputScaffoldGraphHead(ScaffoldGraphHead* scaffoldGraphHead, ContigSetHead* contigSetHead, char* RscaffoldGraphFile) {
	FILE* fp2;
	if ((fp2 = fopen(RscaffoldGraphFile, "w")) == NULL) {
		printf("%s, does not exist!", RscaffoldGraphFile);
		exit(0);
	}
	for (long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++) {
		for (long int j = 0; j < scaffoldGraphHead->scaffoldGraph[i].edgeCount; j++) {
			fprintf(fp2, "%ld--%d--%d--%lf--%ld--%ld\n", i, scaffoldGraphHead->scaffoldGraph[i].edge[j].nodeIndex, scaffoldGraphHead->scaffoldGraph[i].edge[j].orientationType, scaffoldGraphHead->scaffoldGraph[i].edge[j].weight,*scaffoldGraphHead->scaffoldGraph[i].edge[j].pairedReadCount, *scaffoldGraphHead->scaffoldGraph[i].edge[j].gapDistance);
		}
	}
	fclose(fp2);
 }

*/
void dfs(ScaffoldGraphHead * scaffoldGraphHead, long int index, bool * visited, long int * contigPosition){
	visited[index] = true;
	cout<<index<<"--"<<contigPosition[index]<<endl;
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraph[index].edgeCount; i++){
		if(visited[scaffoldGraphHead->scaffoldGraph[index].edge[i].nodeIndex] != true){
			dfs(scaffoldGraphHead, scaffoldGraphHead->scaffoldGraph[index].edge[i].nodeIndex, visited, contigPosition);
		}
	}
}

void DFS(ScaffoldGraphHead * scaffoldGraphHead, long int * contigPosition){
	bool * visited = (bool *)malloc(sizeof(bool)*scaffoldGraphHead->scaffoldGraphNodeCount);
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		visited[i] = false;
	}
	for(long int i = 0; i < scaffoldGraphHead->scaffoldGraphNodeCount; i++){
		if(visited[i] == false){
			cout<<"start:"<<endl;
			dfs(scaffoldGraphHead, i, visited, contigPosition);
			cout<<endl;
		}
	}
	
}







#endif
