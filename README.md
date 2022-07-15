# CRCLM


Runtime environment: R 3.61

System requirements: At least 32Gb RAM and 32 core CPU

R Package requirements:

1. Seurat 3.2.0 Single cell analysis cite: Butler et al., Nature Biotechnology 2018.
2. monocle 2.16.0 Pseudotime analysis cite: Trapnell C, Cacchiarelli D, Grimsby J, Pokharel P, Li S, Morse M, Lennon NJ, Livak KJ, Mikkelsen TS, Rinn JL (2014). “The dynamics and regulators of cell fate decisions are revealed by pseudo-temporal ordering of single cells.” Nature Biotechnology.
3. dplyr 1.0.2 data processing
4. tidyr 1.1.2 data processing
5. reshape2 1.4.4 Data processing

Install and Run method:

1. Place each file in the CRCLM/ folder under the same folder named "work_folder/"
2. Create two folders "outs/filtered_feature_bc_matrix/" and "outputs/" in work_folder/. Put the read file of 10X single cell sequencing into
outs/filtered_feature_bc_matrix/
3. Install all dependent packages in R, and run R code
4. The results of the operation will be automatically output to outputs/
5. If want to use the data that has been analyzed to perform other analyses in our paper, please run corresponding R files, please also find the dependent files in CRCLM/.

License:

This software is licensed for any non-commercial use after the official paper is published.
