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

int main(int argc,char** argv) {   //argc -- �����в��� �ܸ��������� ��ִ�г�������argv[i] --�� i ������
	int maxSize = 2000;
	char * line = (char *)malloc(sizeof(char)*maxSize);   //��̬����һ��char���ʹ�С���ڴ�
	
	char * resultOutPutDirectory = (char *)malloc(sizeof(char)*50);
	strcpy(resultOutPutDirectory, "./Result-OutPut-Directory/");
	int contigLengthThreshold =80000;
	int minPairedReadCount =10;
	double ratio = 1;
	char * alignmentResultFile = NULL;
	char * sortalignmentResultFile = NULL;
  char * scaffoldGraphFile = NULL;
  char * RscaffoldGraphFile = NULL;
  
	//�ַ��͵Ŀ�ָ�룬Ϊ�˸�ָ��һ�����䣬��ָֹ����һЩ��Ҫ�Ŀռ䣬�����˴����������Ҫ�ļ���Ӱ��
	char aligning = false;
	
	char * contigSetFile = NULL;//�ַ��͵Ŀ�ָ��
	char * leftReadBamFile = NULL;//�ַ��͵Ŀ�ָ��
	char * rightReadBamFile = NULL;//�ַ��͵Ŀ�ָ��
	int ch = 0;
	while ((ch = getopt(argc, argv, "c:r:l:o:t:m:p:a:")) != -1) {
		switch (ch) {//swith�ж��ַ����
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
 
	cout<<"Generate comparison information file��alignmentResultFile.fa"<<endl;
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
