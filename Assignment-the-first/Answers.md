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
    
    2. Based on a 2016 paper by Wright and Vetsigian, a quality score cutoff of 26 for index reads was the best tradeoff for maintaining correct reads and removing invalid reads. This value may vary per individual dataset, and they suggest a way of emperically finding the best value for your own dataset based on the number of reads and the number of mis-assignments, but I probably won't perform that on this data set as that would require all the data to be read multiple times. 
    >Wright, E.S., Vetsigian, K.H. Quality filtering of Illumina index reads mitigates sample cross-talk. BMC Genomics 17, 876 (2016). https://doi.org/10.1186/s12864-016-3217-x
   
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
