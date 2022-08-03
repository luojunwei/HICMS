#ifndef CONTIG_H_INCLUDED 
#define CONTIG_H_INCLUDED 


#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

using namespace std;


typedef struct Contig{
	char * contig;//结构体指针就是指向结构体变量的指针，表示的是这个结构体变量占内存中的起始位置，
    int contigLength;   //如果把一个结构体变量的起始地址存放在一个指针变量中，那么，这个指针变量就指向该结构体变量	
}Contig;//Contig就是结构体

typedef struct ContigSetHead{
	Contig * contigSet;
	long int contigCount;	
}ContigSetHead;//ContigSetHead就是结构体




ContigSetHead * GetContigSet(char * contigSetFile, int contigLengthThreshold);
//定义指向ContigSetHead类型数据的指针变量GetContigSet函数，GetContigSet有两个参数，首先它是一个函数，只不过这个函数的返回值是一个地址值。
//函数返回值必须用同类型的指针变量来接受，也就是说，指针函数一定有函数返回值，而且，在主调函数中，函数返回值必须赋给同类型的指针变量

void DeleteTailOfContigSet(ContigSetHead * contigSetHead, int tailLength);
//定义了一个函数，包含两个参数

char * ReverseComplement(char * temp);
//定义了一个函数,只不过这个函数的返回值是一个地址值,有一个参数。

void OutputContigSet(ContigSetHead * contigSetHead);
//定义了一个函数，包含一个参数。




































#endif