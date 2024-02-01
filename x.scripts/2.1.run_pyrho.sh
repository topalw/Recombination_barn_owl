# 1. make table 
# for other pops remove --moran_pop_size argument 
# and provide number of haplotypes (26) in sample size 

pyrho make_table --samplesize 152 --approx --moran_pop_size 200 \
 --numthreads 12 --mu 4.6e-9 --outfile ../1.data/2.pyrho/1.make_table/ch_m200 \
 --smcpp_file ../1.data/2.pyrho/smcpp/3.plots/CH.csv  

# 2. hyperparam 
# run for every pop and set sample size accordingly 

pyrho hyperparam --samplesize 152 --tablefile ../1.data/2.pyrho/1.make_table/ch_m200 --mu 4.6e-9 --ploidy 2 --smcpp_file CH.csv --outfile ../1.data/2.pyrho/2.hyperparam/ch_m200.hyperparam --numthreads 12 --windowsize 30,40,50,60,70,80,90 --blockpenalty 15,20,25,30,35,40,45,50

# 3. optimize
# run this for every scaffold 
# run this for every pop

ss="Super-Scaffold_1"

vcf="../0.data/CH_*${ss}.vcf.gz"
out="../3.optimize/CH/CH_m200_${ss}"
win=70
bp=15

pyrho optimize --vcffile ${vcf} --windowsize ${win} --blockpenalty ${bp} --tablefile ../1.data/2.pyrho/1.make_table/ch_m200 --ploidy 2 --numthreads 12 --outfile ${out}
