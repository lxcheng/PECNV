# PECNV
code for the paper "Detection of copy number variants and loss of heterozygosis from impure tumor samples using whole exome sequencing data"
workflow(namelist,out_path,GC_name,map_name);
1,namelist
here namelist stands for a file containg absolute path of read count and BAF file,for example:
C:\Users\normal.count  C:\Users\normal.depth
C:\Users\s1.count  C:\Users\s1.depth
C:\Users\s2.count  C:\Users\s2.depth
The first line stands for the path of normal samples,separated by tab.
While the following lines contain the matched tumor samples.

2,out_path
out_path specify the output path;
 
3,GC_name,map_name
GC_name,map_name stands for the filename of GC content and mappbility ,which can be downloaded from UCSC.
