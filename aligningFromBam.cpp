#ifndef aligningFromBam_CPP_INCLUDED 
#define aligningFromBam_CPP_INCLUDED 
#include <cstdlib>
#include <iostream>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include<vector>


#include "aligningFromBam.h"

using namespace std;

AligningResultHead * GetAligningResultHead(char * leftBamFile, char * rightBamFile, char * alignmentResultFile, char *  sortalignmentResultFile, ContigSetHead * contigSetHead, bool aligning){
	
	FILE * fp;
	long int aligningCount = 0;	
	if(aligning == false){
		if((fp = fopen(alignmentResultFile, "w")) == NULL){
			printf("%s, does not exist!", alignmentResultFile);
			exit(0);
		}
		BamReader bamReaderLeft;
		BamReader bamReaderRight;

		bamReaderLeft.Open(leftBamFile);
		bamReaderRight.Open(rightBamFile);

		BamAlignment alignmentLeft;
		BamAlignment alignmentRight;
		
		long int minScore = 10;
		
		std::vector< int > clipSizes;
		std::vector< int > readPositions;
		std::vector< int > genomePositions;

		long int leftRefID = -1;
		long int rightRefID = -1;
		long int leftPosition = -1;
		long int rightPosition = -1;
		long int leftOrientation = -1;
		long int rightOrientation = -1;
		
		bool token = false;
		
		while(bamReaderLeft.GetNextAlignment(alignmentLeft) && bamReaderRight.GetNextAlignment(alignmentRight)){

			while((alignmentLeft.AlignmentFlag & 0x900) != 0){
				token = false;
				bamReaderLeft.GetNextAlignment(alignmentLeft);
				continue;
			} 
			while((alignmentRight.AlignmentFlag & 0x900) != 0){
				token = false;
				bamReaderRight.GetNextAlignment(alignmentRight);
				continue;
			}
	
			
			//if(alignmentLeft.IsMapped() && alignmentRight.IsMapped() && alignmentLeft.RefID != alignmentRight.RefID && alignmentLeft.MapQuality > minScore && alignmentRight.MapQuality > minScore && alignmentLeft.GetSoftClips(clipSizes, readPositions, genomePositions) != true && alignmentRight.GetSoftClips(clipSizes, readPositions, genomePositions) != true){
			if(alignmentLeft.IsMapped() && alignmentRight.IsMapped() && alignmentLeft.RefID != alignmentRight.RefID && alignmentLeft.MapQuality > minScore && alignmentRight.MapQuality > minScore){

				fprintf(fp, "%d,%d,%d,%d,%d,%d\n", alignmentLeft.RefID, alignmentRight.RefID, alignmentLeft.Position, alignmentRight.Position, alignmentLeft.IsReverseStrand(), alignmentRight.IsReverseStrand());
				aligningCount++;
			}
		}
		
		cout<<"aligningCount:"<<aligningCount<<endl;
		fflush(fp);
		fclose(fp);

	}
    aligningCount = 0;
	//exit(0);
	
	
	if((fp = fopen(alignmentResultFile, "r")) == NULL){
        printf("%s, does not exist!", alignmentResultFile);
        exit(0);
    }
	char * p;
	const char * split = ","; 
	long int maxSize = 10000;
	char * line = (char *)malloc(sizeof(char)*maxSize);
  long int i = 0;
	while((fgets(line, maxSize, fp)) != NULL){ 
		aligningCount++;
	}
	fclose(fp);

	AligningResultHead * aligningResultHead = (AligningResultHead *)malloc(sizeof(AligningResultHead));
  	aligningResultHead->aligningCount = aligningCount;
	aligningResultHead->aligningResult = (AligningResult *)malloc(sizeof(AligningResult)*aligningCount);
	
    if((fp = fopen(alignmentResultFile, "r")) == NULL){
        printf("%s, does not exist!", alignmentResultFile);
        exit(0);
    }
	//cout<<"ff"<<endl;
	
	long int temp = 0;
	while((fgets(line, maxSize, fp)) != NULL){ //读取line中的每一行，不为空时
		
		p = strtok(line,split);
		aligningResultHead->aligningResult[i].leftContigIndex = atoi(p);//atoi把字符串转换为整形
		
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[i].rightContigIndex = atoi(p);
		
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[i].leftPosition = atoi(p);
		
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[i].rightPosition = atoi(p);
		
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[i].leftOrientation = atoi(p);
		
		p = strtok(NULL,split);
		aligningResultHead->aligningResult[i].rightOrientation = atoi(p);
		
		if(aligningResultHead->aligningResult[i].leftContigIndex > aligningResultHead->aligningResult[i].rightContigIndex){
			temp = aligningResultHead->aligningResult[i].leftContigIndex;//左边的id大于右边时，把左边id赋值给temp
			aligningResultHead->aligningResult[i].leftContigIndex = aligningResultHead->aligningResult[i].rightContigIndex;
			aligningResultHead->aligningResult[i].rightContigIndex = temp;//左边的id等于右边时，把右边id赋值给temp
			
			temp = aligningResultHead->aligningResult[i].leftPosition;
			aligningResultHead->aligningResult[i].leftPosition = aligningResultHead->aligningResult[i].rightPosition;
			aligningResultHead->aligningResult[i].rightPosition = temp;
			
			temp = aligningResultHead->aligningResult[i].leftOrientation;
			aligningResultHead->aligningResult[i].leftOrientation = aligningResultHead->aligningResult[i].rightOrientation;
			aligningResultHead->aligningResult[i].rightOrientation = temp;
			//对每一行的ID进行排序，使左边的ID小于右边的ID
		}
		
		i++;
	}
	
	cout<<"Sort the Aligning results"<<endl;
	sortAligningResult(aligningResultHead->aligningResult, 0, aligningCount - 1);

	//OptimizeAligningResultHead(aligningResultHead);
	//OutputAligningResultHead(aligningResultHead, contigSetHead);
	//exit(0);
 
 //排序文件
/*
AligningResult* a =  aligningResultHead->aligningResult;
FILE * fp1;
if ((fp1 = fopen(sortalignmentResultFile, "w")) == NULL) {
		printf("%s, does not exist!", sortalignmentResultFile);
		exit(0);
	}
 for (int i = 0; i < aligningCount; i++)
{
	fprintf(fp1, "%d,%d,%ld,%ld,%d,%d\n", aligningResultHead->aligningResult[i].leftContigIndex, aligningResultHead->aligningResult[i].rightContigIndex, aligningResultHead->aligningResult[i].leftPosition, aligningResultHead->aligningResult[i].rightPosition, aligningResultHead->aligningResult[i].leftOrientation, aligningResultHead->aligningResult[i].rightOrientation);
}
	fflush(fp1);
	fclose(fp1);
	cout<<"Generate sortalignmentResultFile.fa !"<<endl;    
 */ 
	return aligningResultHead;
  
}

long int MIN(long int a, long int b){
	if(a > b){
		return b;
	}else{
		return a;
	}
}

long int MAX(long int a, long int b){
	if(a < b){
		return b;
	}else{
		return a;
	}
}


void sortAligningResult(AligningResult * a, long int left, long int right){
   /* FILE * fp; //*fp鏄寚鍚戞枃浠剁粨鏋勪綋鐨勬寚閽堝彉閲忥紝閫氳繃fp鍙壘鍒板瓨鏀炬煇涓枃浠朵俊鎭殑缁撴瀯鍙橀噺锛屾牴鎹繖涓粨鏋勫彉閲忕殑淇℃伅鎵惧埌璇ユ枃浠讹紝瀹炴柦瀵规枃浠剁殑鎿嶄綔銆俧p閫氬父琚垚涓轰竴涓寚鍚戞枃浠剁殑鎸囬拡銆?    if((fp = fopen(contigSetFile, "r")) == NULL){
        printf("%s, does not exist!", contigSetFile);
        exit(0);
    }*/
//FILE *fp = fopen("a.txt", "a+")
	if(left >= right){
        return ;
    }

    long int i = left;
    long int j = right;
    long int leftKey = a[left].leftContigIndex;
   	long int rightKey = a[left].rightContigIndex;
  	long int leftPosition = a[left].leftPosition;
	  long int rightPosition = a[left].rightPosition;
    bool leftOrientation = a[left].leftOrientation;
  	bool rightOrientation = a[left].rightOrientation;
	
    while(i < j){
        while(i < j && (MIN(leftKey, rightKey) < MIN(a[j].leftContigIndex, a[j].rightContigIndex) || (MIN(leftKey, rightKey) == MIN(a[j].leftContigIndex, a[j].rightContigIndex) && MAX(leftKey, rightKey) < MAX(a[j].leftContigIndex, a[j].rightContigIndex)))){
            j--;
        }
		
		if(i < j){
			a[i].leftContigIndex = a[j].leftContigIndex;
			a[i].rightContigIndex = a[j].rightContigIndex;
			a[i].leftPosition = a[j].leftPosition;
			a[i].rightPosition = a[j].rightPosition;
			a[i].leftOrientation = a[j].leftOrientation;
			a[i].rightOrientation = a[j].rightOrientation;
			i++;
		}
      
        while(i < j && (MIN(leftKey, rightKey) > MIN(a[i].leftContigIndex, a[i].rightContigIndex) || (MIN(leftKey, rightKey) == MIN(a[i].leftContigIndex, a[i].rightContigIndex) && MAX(leftKey, rightKey) > MAX(a[i].leftContigIndex, a[i].rightContigIndex)))){
            i++;
        }
		
		if(i < j){
			a[j].leftContigIndex = a[i].leftContigIndex;
			a[j].rightContigIndex = a[i].rightContigIndex;
			a[j].leftPosition = a[i].leftPosition;
			a[j].rightPosition = a[i].rightPosition;
			a[j].leftOrientation = a[i].leftOrientation;
			a[j].rightOrientation = a[i].rightOrientation;
			j--;
		}
    }

    a[i].leftContigIndex = leftKey;
	a[i].rightContigIndex = rightKey;
	a[i].leftPosition = leftPosition;
	a[i].rightPosition = rightPosition;
	a[i].leftOrientation = leftOrientation;
	a[i].rightOrientation = rightOrientation;
	
    sortAligningResult(a, left, i - 1);
    sortAligningResult(a, i + 1, right); 
		
	//fprintf(fp, "%d,%d,%d,%d,%d,%d\n", alignmentLeft.RefID, alignmentRight.RefID, alignmentLeft.Position, alignmentRight.Position, alignmentLeft.IsReverseStrand(), alignmentRight.IsReverseStrand());
	//fclose(fp)
}

void OutputAligningResultHead(AligningResultHead * aligningResultHead, ContigSetHead * contigSetHead){
	
	long int start = 0;
	long int end = -1;
	long int leftContigIndex = -1;
	long int rightContigIndex = -1;
	
	long int m = 0;
	long int n = 0;
	
	long int temp = 0;
	
	for(long int i = 0; i < aligningResultHead->aligningCount; i++){

		if(aligningResultHead->aligningResult[i].leftContigIndex == 2725 && aligningResultHead->aligningResult[i].rightContigIndex == 2726){
		
			cout<<aligningResultHead->aligningResult[i].rightContigIndex<<","<<aligningResultHead->aligningResult[i].leftContigIndex<<","<<aligningResultHead->aligningResult[i].rightPosition<<","<<aligningResultHead->aligningResult[i].leftPosition<<","<<aligningResultHead->aligningResult[i].rightOrientation<<","<<aligningResultHead->aligningResult[i].leftOrientation<<endl;
		}
		
	}
	
}

void OptimizeAligningResultHead(AligningResultHead * aligningResultHead){
	
	int * type = (int *)malloc(sizeof(int)*4);
	long int leftContigIndex = -1;
	long int rightContigIndex = -1;
	long int count = 0;
	long int startIndex = -1;
	long int typeIndex = -1;
	long int max = 0;
	
	for(long int i = 0; i < aligningResultHead->aligningCount; i++){
		if(aligningResultHead->aligningResult[i].leftContigIndex == leftContigIndex && aligningResultHead->aligningResult[i].rightContigIndex == rightContigIndex){
			if(aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 0){
				type[0]++;
			}
			if(aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 1){
				type[1]++;
			}
			if(aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 0){
				type[2]++;
			}
			if(aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 1){
				type[3]++;
			}
			count++;
		}else if(aligningResultHead->aligningResult[i].leftContigIndex == rightContigIndex && aligningResultHead->aligningResult[i].rightContigIndex == leftContigIndex){
			if(aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 0){
				type[0]++;
			}
			if(aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 1){
				type[2]++;
			}
			if(aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 0){
				type[1]++;
			}
			if(aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 1){
				type[3]++;
			}
			count++;
		}else{
			if(count > 1){
				max = 0;
				for(long int j = 0; j < 4; j++){
					if(type[j] > max){
						max = type[j];
						typeIndex = j;
					}
				}
				
				for(long int j = startIndex; j < startIndex + count; j++){
					if(aligningResultHead->aligningResult[j].leftContigIndex == leftContigIndex && aligningResultHead->aligningResult[j].rightContigIndex == rightContigIndex){
						if(aligningResultHead->aligningResult[j].leftOrientation == 0 && aligningResultHead->aligningResult[j].rightOrientation == 0){
							if(type[0] != max){
								aligningResultHead->aligningResult[j].leftContigIndex = -1;
							}
						}
						if(aligningResultHead->aligningResult[j].leftOrientation == 0 && aligningResultHead->aligningResult[j].rightOrientation == 1){
							if(type[1] != max){
								aligningResultHead->aligningResult[j].leftContigIndex = -1;
							}
						}
						if(aligningResultHead->aligningResult[j].leftOrientation == 1 && aligningResultHead->aligningResult[j].rightOrientation == 0){
							if(type[2] != max){
								aligningResultHead->aligningResult[j].leftContigIndex = -1;
							}
						}
						if(aligningResultHead->aligningResult[j].leftOrientation == 1 && aligningResultHead->aligningResult[j].rightOrientation == 1){
							if(type[3] != max){
								aligningResultHead->aligningResult[j].leftContigIndex = -1;
							}
						}
					}else if(aligningResultHead->aligningResult[j].leftContigIndex == rightContigIndex && aligningResultHead->aligningResult[j].rightContigIndex == leftContigIndex){
						if(aligningResultHead->aligningResult[j].leftOrientation == 0 && aligningResultHead->aligningResult[j].rightOrientation == 0){
							if(type[0] != max){
								aligningResultHead->aligningResult[j].leftContigIndex = -1;
							}
						}
						if(aligningResultHead->aligningResult[j].leftOrientation == 0 && aligningResultHead->aligningResult[j].rightOrientation == 1){
							if(type[2] != max){
								aligningResultHead->aligningResult[j].leftContigIndex = -1;
							}
						}
						if(aligningResultHead->aligningResult[j].leftOrientation == 1 && aligningResultHead->aligningResult[j].rightOrientation == 0){
							if(type[1] != max){
								aligningResultHead->aligningResult[j].leftContigIndex = -1;
							}
						}
						if(aligningResultHead->aligningResult[j].leftOrientation == 1 && aligningResultHead->aligningResult[j].rightOrientation == 1){
							if(type[3] != max){
								aligningResultHead->aligningResult[j].leftContigIndex = -1;
							}
						}
					}
				}
				
			}
			
			type[0] = 0;
			type[1] = 0;
			type[2] = 0;
			type[3] = 0;
			leftContigIndex = aligningResultHead->aligningResult[i].leftContigIndex;
			rightContigIndex = aligningResultHead->aligningResult[i].rightContigIndex;
			startIndex = i;
			count = 1;
			if(aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 0){
				type[0]++;
			}
			if(aligningResultHead->aligningResult[i].leftOrientation == 0 && aligningResultHead->aligningResult[i].rightOrientation == 1){
				type[1]++;
			}
			if(aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 0){
				type[2]++;
			}
			if(aligningResultHead->aligningResult[i].leftOrientation == 1 && aligningResultHead->aligningResult[i].rightOrientation == 1){
				type[3]++;
			}
			
		}
		
	}
	
	
	count = 0;
	for(long int i = 0; i < aligningResultHead->aligningCount; i++){
		if(aligningResultHead->aligningResult[i].leftContigIndex == -1){
			count++;
			continue;
		}else{
			aligningResultHead->aligningResult[i - count].leftContigIndex = aligningResultHead->aligningResult[i].leftContigIndex;
			aligningResultHead->aligningResult[i - count].rightContigIndex = aligningResultHead->aligningResult[i].rightContigIndex;
			aligningResultHead->aligningResult[i - count].leftPosition = aligningResultHead->aligningResult[i].leftPosition;
			aligningResultHead->aligningResult[i - count].rightPosition = aligningResultHead->aligningResult[i].rightPosition;
			aligningResultHead->aligningResult[i - count].leftOrientation = aligningResultHead->aligningResult[i].leftOrientation;
			aligningResultHead->aligningResult[i - count].rightOrientation = aligningResultHead->aligningResult[i].rightOrientation;
		}
	}
	aligningResultHead->aligningCount = aligningResultHead->aligningCount - count;
	
	
}






#endif
