# HICMS
HICMS：a scaffolding method based on Hi-C contact Matrix
=========


Scaffolder: HICMS
=================

1) Introduction
```
    HICMS is an scaffolder which according to hic comparison informationaims to determine the orientations and orders of contigs. 
    The contigs can be produced by any assembler.
    The input data of Two bam files after bwa alignment(bam format),and the contigs (fasta format). 
```
2) Before installing and running
```
    Please install BWA from https://github.com/lh3/bwa.
	Please install Samtools from https://sourceforge.net/projects/samtools/files/samtools/.
	Please build and install Bamtools from https://github.com/pezmaster31/bamtools.
```

3) Installing.
```
    HICMS should run on Linux operating sysetm with gcc. We test HICMS using gcc9.4.0 on Ubuntu.
    Create a main directory HICMS . Copy all source code to this directory.
	cd HICMS
	export BAMTOOLS_HOME_INCLUDE=/path_bamtools_include_api_shared/
	export BAMTOOLS_HOME_LIB=/path_bamtools_lib_libbamtools.a/
	make all
```

4) Running.
```
	step 1: bwa index contigs.fasta
	step 2: bwa mem -t 10 contigs.fasta Hicreads1.fastq > Hicalign1.sam
	step 3: samtools view -Sb Hicalign1.sam > Hicalign1.bam
	step 4: bwa mem -t 10 contigs.fasta Hicreads2.fasta > Hicalign2.sam
	step 5: samtools view -Sb Hicalign2.sam > Hicalign2.bam
	step 6: HICMS -c contig -l Hicalign1.bam -r Hicalign2.bam -o resultOutPutDirectory
	
	-c <contigs.fa>: 
	    The file includes contigs produced by one assembler.
	-l <align1.bam>:
	    The aligning result between the contigs and the Hi-C left reads(Hicalign1.bam).
	-r <align2.bam>:
	    The aligning result between the contigs and the Hi-C right reads(Hicalign2.bam).
	-o <outdir>: 
	   The output path of the file.
		

5) Output.
```
    The output file "scaffold_set.fa" is the scaffolding result. 
```


