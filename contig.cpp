#ifndef CONTIG_CPP_INCLUDED 
#define CONTIG_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "contig.h"

using namespace std;

ContigSetHead * GetContigSet(char * contigSetFile, int contigLengthThreshold){
   
    ContigSetHead * contigSetHead = (ContigSetHead *)malloc(sizeof(ContigSetHead));
	
    contigSetHead->contigSet = NULL;

    contigSetHead->contigCount = 0;
  
    long int maxSize = 90000;
    char * contig = NULL;
	
    if(NULL == (contig = (char*)malloc(sizeof(char)*maxSize))){
    		perror("malloc�?error!");                          
        exit(1);
    }                                                         

    
    FILE * fp; 
    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }
	
    while((fgets(contig, maxSize, fp)) != NULL){                                         
       if(contig[0] == '>'){  
           contigSetHead->contigCount++; 
       }  
    }  
    fclose(fp);
    
    contigSetHead->contigSet = (Contig *)malloc(sizeof(Contig)*contigSetHead->contigCount);

    for(long int i = 0; i < contigSetHead->contigCount; i++){
        contigSetHead->contigSet[i].contig = NULL;
        contigSetHead->contigSet[i].contigLength = 0;
    }
    //初始化contigset
    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }//当contig文件是空的，输出不存�?
    
    long int allocateLength = 0;
    long int contigIndex = -1;
    while((fgets(contig, maxSize, fp)) != NULL){
       if(contig[0] == '>'){

		   if(strlen(contig) == maxSize-1){              
               while((fgets(contig, maxSize, fp)) != NULL){
                   if(strlen(contig) != maxSize-1){
                       break;
                   }
               }        
           }
//读取文件的每一�?
		   contigIndex++;//建立索引
		   long int len = strlen(contig);
           continue;
       }
       
       
       long int extendLength = strlen(contig);
	   //获取contig的长�?
	  
       if(contig[extendLength-1] == '\n' || contig[extendLength-1] == '\t' || contig[extendLength-1] == '\r'){
           extendLength--;
       }//字符串结束符'\0',转义字符以反斜线"\"开头，后跟一个或几个字符。转义字符具有特定的含义，不同于字符原有的意义，故称“转义”字�?
	   //
	   
	   /*
	    * if(contig[extendLength-1] == '\n' || contig[extendLength-1] == '\t' || contig[extendLength-1] == '\r'){
           extendLength--;   
       }
       */
	   
       long int contigLength = 0;
       char * tempContig = NULL;
       if(contigSetHead->contigSet[contigIndex].contig != NULL){
           if(contigSetHead->contigSet[contigIndex].contigLength + extendLength >= allocateLength){//当空间长�?contig长度>=临时长度�?
		       contigLength = contigSetHead->contigSet[contigIndex].contigLength;    
               contigSetHead->contigSet[contigIndex].contig = (char *)realloc(contigSetHead->contigSet[contigIndex].contig, allocateLength + maxSize + 1);

               allocateLength = allocateLength + maxSize + 1;
               
               strncpy(contigSetHead->contigSet[contigIndex].contig + contigLength, contig, extendLength);
               contigSetHead->contigSet[contigIndex].contig[contigLength + extendLength] = '\0';    
               contigSetHead->contigSet[contigIndex].contigLength = contigLength + extendLength;
                       
           }else{
               strncpy(contigSetHead->contigSet[contigIndex].contig + contigSetHead->contigSet[contigIndex].contigLength, contig, extendLength);
               contigSetHead->contigSet[contigIndex].contig[contigSetHead->contigSet[contigIndex].contigLength + extendLength] = '\0';
               contigSetHead->contigSet[contigIndex].contigLength = contigSetHead->contigSet[contigIndex].contigLength + extendLength;
           }   
           
       }else{
           contigSetHead->contigSet[contigIndex].contig = (char *)malloc(sizeof(char)*(maxSize+1));
           strncpy(contigSetHead->contigSet[contigIndex].contig, contig, extendLength);
           contigSetHead->contigSet[contigIndex].contig[extendLength] = '\0';
           contigSetHead->contigSet[contigIndex].contigLength = extendLength;
           allocateLength = maxSize + 1;
       }  
    }  
	
    fflush(fp);
    fclose(fp);
	if(contigSetHead->contigCount <= 0){
		cout<<"Contig count is zero, please check the contig file!"<<endl;
		exit(0);
	}
	
	long int minContigCount = 0;
	for(long int i = 0; i < contigSetHead->contigCount; i++){
		if(contigSetHead->contigSet[i].contigLength < 100000){
			minContigCount++;
		}
	}
	cout<<"minContigCount(contigLength < 100000):"<<minContigCount<<endl;

  return contigSetHead;
}


void DeleteTailOfContigSet(ContigSetHead * contigSetHead, int tailLength){
	for(int i = 0; i < contigSetHead->contigCount; i++){
		char * tempContig = (char *)malloc(sizeof(char)*(contigSetHead->contigSet[i].contigLength - 2*tailLength + 1));
		strncpy(tempContig, contigSetHead->contigSet[i].contig + tailLength, contigSetHead->contigSet[i].contigLength - 2*tailLength);
		tempContig[contigSetHead->contigSet[i].contigLength - 2*tailLength] = '\0';
		free(contigSetHead->contigSet[i].contig);
		contigSetHead->contigSet[i].contig = tempContig;
		tempContig = NULL;
		contigSetHead->contigSet[i].contigLength = contigSetHead->contigSet[i].contigLength - 2*tailLength;

	}
}

void OutputContigSet(ContigSetHead * contigSetHead){  
	char * file = (char *)malloc(sizeof(char)*50);
	strcpy(file, "./contigSet.fa");

	FILE * fp; 
    if((fp = fopen(file, "w")) == NULL){
        printf("%s, does not exist!", file);
        exit(0);
    }

	for(long int i = 0; i < contigSetHead->contigCount; i++){
		fprintf(fp, ">%ld--%d/n", i, contigSetHead->contigSet[i].contigLength);
		fprintf(fp, "%s\n",contigSetHead->contigSet[i].contig);
	}
	
}

char * ReverseComplement(char * temp){
    long int len = strlen(temp);
	char * rcTemp = (char *)malloc(sizeof(char)*(len + 1));

    for(long int i = 0; i < len; i++){
        if(temp[i] == 'A' || temp[i] == 'a'){
            rcTemp[len - 1 - i] = 'T';
        }else if(temp[i] == 'T' || temp[i] == 't'){
            rcTemp[len - 1 - i] = 'A';
        }else if(temp[i] == 'G' || temp[i] == 'g'){
            rcTemp[len - 1 - i] = 'C';
        }else if(temp[i] == 'C' || temp[i] == 'c'){
            rcTemp[len - 1 - i] = 'G';
        }else if(temp[i] == 'N' || temp[i] == 'n'){
            rcTemp[len - 1 - i] = 'N';
        } 
    }
    rcTemp[len]='\0';
    return rcTemp;
}











































#endif