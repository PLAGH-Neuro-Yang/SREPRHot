# SREPRHot: a SMOTE and Random grouping strategies-based Ensemble learning model for Protein-RNA binding Hotspot prediction

The performance process includes three steps: feature extraction, feature integration and prediction.

A protein-RNA complex, 4JVH_A, is used as an example to show the process. Here, 'A' is the chain identifer and the experimental complex structure has been solved with PDB code '4JVH'.

SREPRHot uses the following dependencies:
* python 3.8
* numpy
* pandas
* joblib

## Step 1: feature extraction

### 1. Physicochemical characteristics of amino acids

Run "Physicochemical_characteristics.py" with "4JVH_A.pdb" as input file (keep "Physicochemical_characteristics.py" and "4JVH_A.pdb" in the same folder).

```{bash}
python ./Physicochemical_characteristics.py 4JVH_A.pdb
```

Then you can obtain a file called "4JVH_A_PC.csv".

### 2. PSSM
Go to the website https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome to get a pssm file.

Here, select "PSI-BLAST" for "Program Selection" item, and use default parameters for other items.

### 3. CX/DPX
Go to the website https://sourceforge.net/projects/psaia/ to download and install the program "psaia.exe".

a. Run "psaia.exe";

b. Step by step, selecet "Structure Analyser" tab control, then find "Analysis Types" and check both "Analyse as Bound" and "Analyse by Chain". All parameters are set to default;

c. Input the pdb file "4JVH_A.pdb" (hydrogen atoms removed) to the program, and click "run" to get the result, i.e., file "4JVH_A_unbound.tbl" and "4JVH_A_bound.tbl".

### 4. ASA
Go to the website http://www.bioinf.manchester.ac.uk/naccess/ to download and install the program "Naccess V2.1.1" on Linux.

a. Delete the nucleic coordinates to get a unbound pdb file "4JVH_A_unbound.pdb";

b. Run "Naccess V2.1.1" on Linux and input corresponding pdb file (e.g. 4JVH_A.pdb) to get a result file named "4JVH_A_bound.rsa";

c. Repeat step b and get another file "4JVH_A_unbound.rsa".

### 5. Solvent exposure
Get Solvent exposure featues from the website https://sunflower.kuicr.kyoto-u.ac.jp/~sjn/hse/webserver.html.Input query sequence, and select "PSI-BLAST+PSIPRED+AA+W+L" for "SVR Models" item. Copy the results from the e-mail and name it as 4JVH_A_HSE.txt.

### 6. IP
Run "IP.py" with "4JVH_A.pdb" as input file (keep "IP.py" and "4JVH_A.pdb" in the same folder).

```{bash}
python ./IP.py ./4JVH_A.pdb
```

Then you can get a file called "4JVH_A_IP.csv".

### 7. ANN
Go to the website http://sysbio.suda.edu.cn/NACEN/index.html to download "dssp-3.0.0.exe" and R package "NACEN", and see the instruction on the website to perform the following operations.

a. Run R, install the package "NACEN" and library(NACEN);

b. Use function "NACENConstructor" to construct NACEN with the parameter "WeightType = 'SAS'" (node-weighted AAN based on residue solvent accessibility);

c. Calculate by using function "NACEAnanlyzer" and obtain a result file named "4JVH_A_ANN_SAS.txt";

d. Repeat the above procedures with "WeightType = 'SAS'" replaced by "WeightType = 'Hydrophobicity'" (node-weighted AAN based on residue hydrophobicity) to obtain a result file named "4JVH_A_ANN_Hydrophobicity.txt".

## Step 2: feature integration

After all features are extracted, all the feature files obtained above include "4JVH_A_PC.csv", "4JVH_A.asn_matrix.txt", "4JVH_A_unbound.tbl", "4JVH_A_bound.tbl", "4JVH_A_bound.rsa", "4JVH_A_unbound.rsa", "4JVH_A_HSE.txt", "4JVH_A_IP.csv", "4JVH_A_ANN_SAS.txt", "4JVH_A_ANN_Hydrophobicity.txt".

Put all the feature files mentioned above into the "Input" folder where "Feature_integration.py" is located.

Run "Feature_integration.py" with the "Input" folder and "4JVH_A.pdb" in the same folder where "Feature_integration.py" is located.

```{bash}
python ./Feature_integration.py ./4JVH_A.pdb
```

## Step 3: prediction

Put "4JVH_A_feature.csv" in the "Input_data" folder, and keep "predict.py", "Model_file" folder and "Input_data" folder in the same path.

```{bash}
python ./predict.py ./Input_data/4JVH_A_feature.csv
```

The output is "./Result/4JVH_A_hotspot_result.csv" where the predicted hotspots are stored.

It should be pointed out that although SREPRHot can give the prediction results for all the residues, only the results corresponding to the interfacial residues are significant. 

## Help

For any questions, please contact us by chunhuali@bjut.edu.cn.
