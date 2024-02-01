# x.scripts 

This folder contains scripts. 

## LepMAP3 - 1.#

The first file `1.1_run_LepMap.sh` is the pipeline of commands used to obtain the maps and the orders from a filtered VCF.

The rest plot or prune the data.

`1.2_plot_maps.r` plots the different maps for different LOD scores to help you choose the best LOD score.
It involves some manual curation of the results tailored to our study!

`1.3_prune_orders.r` prunes the end of the orders and makes the orders ascending.
It also removes points that are 1/2 SD away from the mean of the absolute residuals of a regression of the physical order on the genetic order.

`1.4_choose_bestL.r` compares the likelihood of each order generated from the 3 runs and finds the best one to include in the final file.

`1.5_compare_functions.r` compares the 3 mapping functions and is not mandatory.

`1.6.gams.r`,`1.6.1.gams_male.r` and `1.6.2.gams_female.r` make ascending gam models for the sex-averaged, male and female maps respectively. 
These models are saved as rds files in this folder. One for each scaffold. 

`1.7.gam_gather.r` gathers the rds files and merges them in one file. 

## pyrho

The first file `2.0.run_SMCpp.sh` has the commands we used to run SMC++. 
The second file `2.1.run_pyrho.sh` has the commands we used to run pyrho.
he parts of the these files were ran separately on an HPC cluster using multithreading.
Tread cautiously

## making windows 

Results of pyrho are on between SNP intervals and thus we create non-overlapping windows along the genome. 
Their lengths can differ from 1kb to 5Mb. We created the windows using bedtools but here we overlap them with the pyrho file and calculate 
the recombination rate in each window using `3.0.make_windows.sh` and `3.0.1.make_windows.r`. 

Because the output of pyrho can be affected by selection and changes in $N_e$ we scale the output to have the total length of the genetic map 
in `3.1.scale.r`. 

Because we used a mappability mask and removed SNPs from regions of the genome we need to take that in account in calculations of pi so we scale
pi estimates in `3.2.scale_pi.r`.

Finally we merge all windowed statistics (pi, recombination rate, TSS, etc) in big masterfiles using `3.3.merge.r`. 

## pre plotting analyses 

Some last steps before plotting stuff are under the number `4` in the script folder. 

Most of them were exploratory analyses so only intermittent stuff remain. 

File `4.5.COs_sexes.r` calculates where the crossover happen in each sex using the LepMap3 raw output of the ordering step.

File `4.8.rintra.r` calculates the ration of intrachromosomal shuffling ($\overline{r}$ or $r_{intra}$) for male and female linkage maps. 

## Plotting 

File `5.0.manuscript_figures.r` makes the $v1/_figureX$ pdfs in `plots/` which are then modified in illustrator to produce the $v2/_$ final versions of figures. 

Filte `5.1.supplementary_tables_figures.r` makes the supplementary figure in `plots/SUPP*` and a couple of tables, in `../2.metadata`. 

## AUX

`functions.r` provides some functions used in other scripts. 
