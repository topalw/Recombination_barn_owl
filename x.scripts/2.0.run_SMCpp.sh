# 1. vcf2smc 
# run for each population 
# and for each scaffold 
# by providing the masked regions 
# and the distinguished individuals per population

pop=CH
ss=Super-Scaffold_1
vcf=`ls ../1.data/0.vcfs/${pop}*${ss}.vcf.gz`
chr_mask="../1.data/2.pyrho/smcpp/chrmasks/${ss}_mask.bed.gz"

while read -r line; # run for each distinguished individual 
do
d=${line} # distinguished individual 
smc_out="../1.data/2.pyrho/smcpp/1.smcfiles/d_${d}_${pop}_${ss}"
smpc++ vcf2smc ${vcf} ${smc_out} -d ${d} ${d} -m ${chr_mask} ${ss} ${pop}:`cat ../1.data/2.pyrho/smcpp/1.smcfiles/${pop}_samples.csv`
done < ../1.data/2.pyrho/smcpp/${pop}.distinguished

# 2. estimate the model 
# we ran it with 12 cores 

outdir="../1.data/2.pyrho/smcpp/2.estimates/${pop}"
files=`ls ../1.data/SMCpp/1.smcfiles/*${pop}_${ss}* `

smc++ estimate --spline piecewise -o ${outdir} --cores 1 4.6e-9 ${files}

# 3. plot 

smc++ plot -c ../1.data/2.pyrho/smcpp/3.plots/${pop}.pdf -c ${outdir}/model.final.json


