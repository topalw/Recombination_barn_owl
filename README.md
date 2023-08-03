# Chapter_A

A directory that contains key data and scripts for my first chapter in case i die, get fired or decide to publish. 

The initial division is in data in `1.data` and scripts in `x.scripts`. 

## x.scripts

These files are probably the only ones kept in github since data are a bit too big to include here. 

The division is clear. 

- All files that start with `1.` work with LepMAP3 data. 
- All files that start with `2.` work with pyrho data. 
- All other files are miscellaneous scripts used in testing or as parts of other scripts.

### LepMAP3

The first file is the pipeline of commands used to obtain the maps and the orders from a filtered VCF. 
The rest plot or prune the data. 

`1.2_plot_maps.r` plots the different maps for different LOD scores to help you choose the best LOD score. 
It involves some manual curation of the results ! 

`1.3_prune_orders.r` prunes the end of the orders and makes them ascending. 

It also removes points that are 1/2 SD away from the mean of the absolute residuals of a regression of the physical order on the genetic order.
This complicated sentence is it removes SNPs that are placed with uncertainty. 

`1.4_choose_bestL.r` compares the likelihood of each order from the 3 runs and finds the best one to include in the final file. 

`1.5_compare_functions.r` compares the 3 mapping functions and is not mandatory. 

### pyrho

blahblah

## 1.data

The directory contains data. 

It has a folder for LepMAP `1.LepMAP`. 

This in turn contains the `maps` and `orders`. 

### maps 

Contains the LG maps created with LepMAP3. 

The most useful thing here is the `Call2_tf1_maf05_miss05_1kb_pruned_f.call_lod15_theta0.03_map_manually_pruned.txt` that contains the map used to create the orders. 

For a description of how these maps came to be check my mid-thesis. 

### orders

3 different ordering runs (`1st.ordering`, `2nd.ordering`, `3rd.ordering`) and 2 with different mapping functions (Haldane is default) (`other_functions`). 
Wish i had a consistent naming system. 

This also contains the `likelihoods` folder that contains likelihoods of each order for each LG. 


