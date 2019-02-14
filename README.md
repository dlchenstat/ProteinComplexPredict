# ProteinComplexPredict
R code developed by Donglai Chen for manuscript "A Label Free Mass spectrometry Method to Predict Endogenous Protein Complex Composition".

This code returns cluster IDs for split fitted profile clustering for concatenated profiles of two IEX alone, two SEC alone and IEX+SEC.

Execute file "ProteinComplexPredictFunctions.r" first and then use "ProteinComplexPredictMain.r" for clustering results.

Input files are two IEX datasets, two SEC datasets, list of SEC cytosolic proteins, list of protein contaminents in IEX, (if proteins are already filtered, those two lists are not needed) Gaussian fitting results of four datasets, list of known protein complexes.

# List of input files
The IEX input file
- IEX_bio1_common_cytosol.csv
- IEX_bio2_common_cytosol.csv
- CytoContaminents.csv (remove low quality IEX proteins)

The IEX peak detection file (generated from Gaussian fitting Matlab code)
- peakloc-iex-bio1-uma-2015aug.csv
- peakloc-iex-bio2-uma-2015aug.csv

The SEC input file
- SEC_Bio1_nov.csv
- SEC_Bio2_nov.csv
- SEC_Bio1_Bio2_cytosol_list.csv (select cytosolic proteins in SEC)

The SEC peak detection file  (generated from Gaussian fitting Matlab code)
- peakloc-sec1-uma-2015nov.csv
- peakloc-sec2-uma-2015nov.csv
