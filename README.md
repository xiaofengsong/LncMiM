# LncMiM

# Introduction
Increasing lncRNA-associated competing triplets were found to play important roles in cancers. With the accumulation of high-throughput sequencing data in public databases, the big size of available tumor samples introduces new challenges to identify competing triplets. Here, we present a novel method, called LncMiM, to identify the competing triplets for large-scale expression data of lncRNA, miRNA, and mRNA. In LncMiM, non-linear correlation analysis is used to cover the problem of weak correlations between miRNA-target pairs, which is due to the difference in the magnitude of the expression level. In addition, besides the miRNA, the impact of lncRNA and mRNA on the interactions in triplets are also considered to improve the identification sensitivity of LncMiM without reducing its accuracy. Original codes and detailed manuals are included in the downloading URL: https://github.com/xiaofengsong/LncMiM/archive/v1.0.0.tar.gz.

Authors: Jian Zhao (zhaojian@nuaa.edu.cn), Xiaofeng Song (xfsong@nuaa.edu.cn) <br>
Maintainer: Jian Zhao (zhaojian@nuaa.edu.cn)


# Usage
>* ### Software Prerequisites:
>> LncMiM is implemented in R (version 3.3.1) under Linux system.

>* ### Command: 
>>      Rscript LncMiM.r miRNA_expr lncRNA_expr mRNA_expr candidate_triplet outputfile mode widow_size cor_method  
	
>* ### Arguments:
>>      1)  miRNA_expr
>>>		miRNA expression file
>>      2)  lncRNA_expr
>>>		lncRNA expression file
>>      3)  mRNA_expr
>>>		mRNA expression file
>>      4)  candidate_triplet
>>>		Candidate triplets (miRNA lncRNA  mRNA)
>>      5)  outputfile
>>>		Results of the identified competing triplets
>>      6)  mode
>>>		Three modes: miRNA (default), lncRNA, mRNA, corresponding to the candidate tripelt construction method
>>      7)  window_size
>>>		The sliding window size (0~1), default value is 0.25 (one quarter of the total samples).
>>      8)  cor_method
>>>		Three method: pearson, spearman (default), and kendall

# Note
> * 'miRNA-centered' triplets: miRNA is negatively correlated to mRNA and lncRNA. <br>
> * 'lncRNA-centered' triplets: lncRNA is negatively correlated to miRNA and is positively correlated to mRNA. <br>
> * 'mRNA-centered' triplets: mRNA is negatively correlated to miRNA and is positively correlated to lncRNA. <br>
> * Mode: miRNA, lncRNA, and mRNA are respectely corresponding to 'miRNA-centered', 'lncRNA-centered', and 'mRNA-centered' condidate triplets. <br> 

# Example
	Rscript LncMiM.r example/miRNA_expr.txt example/lncRNA_expr.txt example/pc_expr.txt example/candidate_miRNA_centered_triplets.txt result_miRNA_centered_triplets.txt miRNA

	Rscript LncMiM.r example/miRNA_expr.txt example/lncRNA_expr.txt example/pc_expr.txt example/candidate_lncRNA_centered_triplets.txt result_lncRNA_centered_triplets.txt lncRNA

	Rscript LncMiM.r example/miRNA_expr.txt example/lncRNA_expr.txt example/pc_expr.txt example/candidate_mRNA_centered_triplets.txt result_mRNA_centered_triplets.txt mRNA

# Reference
[Jian Zhao, Xiaofeng Song*, Tianyi Xu, Bin Jiang, Jing Wu*. Identification of potential prognostic competing triplets in high-grade serous ovarian cancer. Submitted](https://github.com/xiaofengsong/LncMiM)
# Contact
	Please contact Jian Zhao (zhaojian@nuaa.edu.cn) for questions and comments.

Copyright (C) 2020 Xiaofeng Song.
