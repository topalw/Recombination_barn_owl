# The recombination landscape of the barn owl, from families to populations 

This is a directory that contains key scripts for the publication of "https://doi.org/10.1093/genetics/iyae190".

The initial division is in data in `1.data`, metadata in `2.metadata` and scripts in `x.scripts`. 

The `1.data` is missing from the repository since it contains some files of large size. 

It will be uploaded as a dryad repository at the time of publication. 

With the combination of scripts and data anyone can reproduce the figures and analyses we performed in the paper. 

Some earlier scripts like running LepMap3 and pyrho might require tinkering to run in separate cores or in a cluster due to
the large amount of data they use. 

However, the commands in `x.scripts` and the details in the methods of the paper should suffice in recreating the results. 

If something is amiss let me know. 

## 1.data

Folder found in [this link](https://zenodo.org/records/13982583) as a zenodo repository. More details on the structure can be found there.

The data in the folder contain raw results of LepMap3, pyrho and smc++ along with windowed stats of the assembly etc. 

If thats not raw enought, fastq data can be found in NCBI. Here is the data availability from the manuscript: 
"All sample information is found in Supplementary File 2. Sequence data used in the study from previous publications are available on NCBI under BioProject codes PRJNA700797, PRJNA727915, PRJNA727977, PRJNA774943, and PRJNA925445. Sequence data generated for this study are available on NCBI under BioProject code PRJNA1172395." 

## 2.metadata 

Files here detail sample information and data used for the a couple of supplementary tables (S1 & S2). 

## synteny 

Files to create the synteny plot of Chicken-Barnowl-Chicken. 

## x.scripts 

The division is clear (for me - lol). 

- All files that start with `1.` work with LepMAP3 data. 
- All files that start with `2.` work with pyrho data. 
- All other files are miscellaneous scripts used in testing or as parts of other scripts.
- Files that start with `5` are needed to recreate the figures of the manuscript. 
- There is a folder called plots with ... well plots. 


