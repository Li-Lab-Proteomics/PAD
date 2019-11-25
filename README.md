## PAD：An information resource for drug repurposing by propagation in a heterogeneous multi-level network

PAD (Predictive Database for Drug Repurposing) is a public database for efficient and reliable drug repurposing (http://lilab.life.sjtu.edu.cn:8080/pad/index.html). The predicted ranking results for drug repurposing should be useful for discovering new therapeutic effects for existing drugs. 

This text you see here includes supplementary codes and sample data. We will show you :

  **1. Description of sample data**
  
  **2. Explanation on parameters**
  
  **3. Implementation**

**1. Description of sample data**


| TXT | Description |
| ------ | ------ |
|f_disease_list.txt|encoding of disease |
|f_protein_list.txt|encoding of protein|
|f_drug_list.txt|encoding of drug|
|f_disease_disease_similarity.txt | disease similarity where the first two columns are diseases and the third is corresponding similarity|
|f_drug_disease.txt  | relation between drug and disease where columns are drug and disease|
|f_drug_drug_similarity.txt|drug similarity where the first two columns are drugs and the third is corresponding similarity|
|f_drug_target.txt|interaction between drug and protein where columns are drug and protein|
|f_gene_disease.txt|relation between gene and disease where columns are corresponding gene-encoded protein and disease |
|f_protein_protein.txt|interaction beween proteins where the first two columns are proteins and the third is corresponding interaction|

**2. Explanation on parameters**

There are four parameters in our model: restart probability *r*, weighting parameters of disease, target protein and drug *a*, *b* and *c* (the sum of *a*, *b* and *c* shall be 1), which control the impact of disease, protein and drug nodes in both initial and transmission process.

**3. Implementation**

  - **Protocal**
  
   You can implement the code as follows:
   
   ```sh
    python calculate_predict_all.py [a] [b] [c] [r]
   ```
    
   It should be noticed that `a+b+c=1, r<1`.
   For instance:
    
   ```sh
    python calculate_predict_all.py 0.5 0.4 0.1 0.7
   ```
    
   where `a=0.5 b=0.4 c=0.1 r=0.7`
  - **Output**   
    The output includes following documents:  
    *---Final_matrix_a_b_c_r.txt* which contains the main result
    
       The *Final_matrix_a_b_c_r.txt* contains the resulting probability matrix, in which rows and columns follow the sorted order of disease, protein and durg which are sorted in input files. Based on this matrix, the probability that the drug will appraoch corresponding disease could be obtained.
    
    *---process_a_b_c_r.txt*   which records the process
    


It should be mentioned that if you want to change the absolute path of input data,it should be correctly specified in original file.

**If you have problems using PAD, please contact jing.li@sjtu.edu.cn**
