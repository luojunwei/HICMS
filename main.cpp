#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <malloc.h>
#include <getopt.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>

#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <string.h>
#include <math.h>

#include "contig.h"
#include "aligningFromBam.h"
#include "scaffoldgraph.h"
#include "scaffolding.h"

using namespace std;

int main(int argc,char** argv) {   //argc -- 命令行参数 总个数，包括 可执行程序名。argv[i] --第 i 个参数
	int maxSize = 2000;
	char * line = (char *)malloc(sizeof(char)*maxSize);   //动态申请一个char类型大小的内存
	
	char * resultOutPutDirectory = (char *)malloc(sizeof(char)*50);
	strcpy(resultOutPutDirectory, "./Result-OutPut-Directory/");
	int contigLengthThreshold =80000;
	int minPairedReadCount =10;
	double ratio = 1;
	char * alignmentResultFile = NULL;
	char * sortalignmentResultFile = NULL;
  char * scaffoldGraphFile = NULL;
  char * RscaffoldGraphFile = NULL;
  
	//字符型的空指针，为了给指针一个分配，防止指向了一些重要的空间，进行了错误操作对重要文件的影响
	char aligning = false;
	
	char * contigSetFile = NULL;//字符型的空指针
	char * leftReadBamFile = NULL;//字符型的空指针
	char * rightReadBamFile = NULL;//字符型的空指针
	int ch = 0;
	while ((ch = getopt(argc, argv, "c:r:l:o:t:m:p:a:")) != -1) {
		switch (ch) {//swith判断字符语句
			case 'c': contigSetFile = (char *)(optarg); break;
			case 'r': rightReadBamFile = (char *)(optarg); break;
			case 'l': leftReadBamFile = (char *)optarg; break;
			case 'o': resultOutPutDirectory = (char *)optarg; break;
			case 't': contigLengthThreshold = atoi(optarg); break;
			case 'm': minPairedReadCount = atoi(optarg); break;
			case 'p': ratio = atof(optarg); break;
			case 'a': alignmentResultFile = (char *)(optarg); break;
			default: break; 
		}
	}
	
	if(opendir(resultOutPutDirectory) == NULL){  
		mkdir(resultOutPutDirectory, 0777);       
    }
	
	int fileNameLen = strlen(resultOutPutDirectory);

	ContigSetHead * contigSetHead = GetContigSet(contigSetFile, contigLengthThreshold);
	cout<<"contig count:"<<contigSetHead->contigCount<<endl;
	

	if(alignmentResultFile == NULL){
		alignmentResultFile = (char *)malloc(sizeof(char)*fileNameLen + 50);
		strcpy(alignmentResultFile, resultOutPutDirectory);
		strcat(alignmentResultFile, "/alignmentResultFile.fa");
	}else{
		aligning = true;
	}

 	if(sortalignmentResultFile == NULL){
		sortalignmentResultFile = (char *)malloc(sizeof(char)*fileNameLen + 50);
		strcpy(sortalignmentResultFile, resultOutPutDirectory);
		strcat(sortalignmentResultFile, "/sortalignmentResultFile.fa");
	}
 
 
	if (scaffoldGraphFile == NULL) {
		scaffoldGraphFile = (char*)malloc(sizeof(char) * fileNameLen + 50);
		strcpy(scaffoldGraphFile, resultOutPutDirectory);
		strcat(scaffoldGraphFile, "/scaffoldGraphFile.fa");
	}

	if (RscaffoldGraphFile == NULL) {
		RscaffoldGraphFile = (char*)malloc(sizeof(char) * fileNameLen + 50);
		strcpy(RscaffoldGraphFile, resultOutPutDirectory);
		strcat(RscaffoldGraphFile, "/RscaffoldGraphFile.fa");
	}


 
	cout<<"Start extracting information for comparison"<<endl;
	AligningResultHead * aligningResultHead = GetAligningResultHead(leftReadBamFile, rightReadBamFile, alignmentResultFile, sortalignmentResultFile, contigSetHead, aligning);
 
	cout<<"Generate comparison information file：alignmentResultFile.fa"<<endl;
	ScaffoldGraphHead * scaffoldGraphHead = GetScaffoldGraphHeadFromAligningResultHead(contigSetHead, aligningResultHead, contigLengthThreshold, minPairedReadCount, ratio,scaffoldGraphFile,RscaffoldGraphFile);
 
 
 // cout<<"Out put Scaffold Graph"<<endl;
 // OutputScaffoldGraphHead(scaffoldGraphHead);
 //OutputscaffoldGraphHead(scaffoldGraphHead, contigSetHead, scaffoldGraphFile);
 // OutputScaffoldGraphHead( scaffoldGraphHead,  contigSetHead, contigPosition);
  
 	cout<<"Optimize Scaffold Graph"<<endl;
	OptimizeScaffoldGraph(scaffoldGraphHead, contigSetHead);

	ScaffoldSetHead * scaffoldSetHead = GetScaffoldingResult(scaffoldGraphHead, contigSetHead);
	char * scaffoldTagFileName = (char *)malloc(sizeof(char)*fileNameLen + 50);
	strcpy(scaffoldTagFileName, resultOutPutDirectory);
	strcat(scaffoldTagFileName, "/scaffold");
 
	OutPutScaffoldSet(scaffoldSetHead->scaffoldSet, contigSetHead, scaffoldTagFileName);
	cout<<"Out put final result: scaffold_Set"<<endl;
	
}
