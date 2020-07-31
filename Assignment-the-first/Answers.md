# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read1 |
| 1294_S1_L008_R2_001.fastq.gz | Index1 |
| 1294_S1_L008_R3_001.fastq.gz | Index2 |
| 1294_S1_L008_R4_001.fastq.gz | Read2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.  
    Read1
    ![Read 1](https://github.com/2020-bgmp/demultiplexing-derrik-gratz/blob/master/Assignment-the-first/read1.png)      
    Read2
    ![Read 2](https://github.com/2020-bgmp/demultiplexing-derrik-gratz/blob/master/Assignment-the-first/read2.png)      
    Index1   
    ![Index 1](https://github.com/2020-bgmp/demultiplexing-derrik-gratz/blob/master/Assignment-the-first/index1.png)      
    Index2 
    ![Index 2](https://github.com/2020-bgmp/demultiplexing-derrik-gratz/blob/master/Assignment-the-first/index2.png)  
    2. ```Your answer here```
    3. ```zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | awk 'NR % 4 == 2' | grep 'N' -c```
    
    7304664
    
## Part 2
1. Define the problem
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
