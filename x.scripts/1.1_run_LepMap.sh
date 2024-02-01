# Global parameters 

ped='../1.data/1.LepMAP/data/ped.txt'
vcf='../1.data/0.vcfs/LM3.vcf'
out1='../1.data/1.LepMAP/data/ParentCall.call'
out2='../1.data/1.LepMAP/data/Filtering.call'
out3="../1.data/1.LepMAP/data/SepChr_lod${lod}.map.txt"
out4="../1.data/1.LepMAP/data/Ordered_lg${lg}.txt"
# 1. ParentCall 
# to go from vcf to LM3 format 

java -cp ~/Lep_MAP3/bin ParentCall2 halfSibs=1 data=${ped} vcfFile=${vcf} removeNonInformative=1 > ${out1}

# 2. Filtering 
# Now filter markers based on segregation distortion 

java -cp ~/Lep_MAP3/bin Filtering2 data=${out1} dataTolerance=0.01 removeNonInformative=1 missingLimit=28 familyInformativeLimit=28 > ${out2}

# make a file with snps present in filtered dataset 
awk '(NR>=7)' ${out2} | cut -f 1,2 > ${out2}.snps.txt
# make a file with the distribution of snps from each scaffold that made it 
tail -n +2 ${out2}.snps.txt | cut -f 1 | sort | uniq -c | sed 's/^ *//g' |sed 's/ /,/g' > ${out2}_dist.csv

# 3. Separate Chromosomes to identify Linkage groups 
# might want to run for different lods, we tested 11-21 
# this might take a while ! (we ran with multiple threads) 
lod=15

java -cp ~/Lep_MAP3/bin SeparateChromosomes2 data=${out2} lodLimit=${lod} numThreads=1 sizeLimit=1 theta=0.03 > ${out3}

### makes counts per linkage group 
fixedfile="${out3}.counts.fixed"
if [ -e  $fixedfile ]
then
        echo "Fixed positions found"
else
        sort $out3 | uniq -c | sort > "${out3}.counts"
        tail -n +2 "${out3}.counts" | sed -e 's/^ *//g' > ${fixedfile}
        rm "${out3}.counts"
fi

# 4. Order markers per linkage group 
lg=1

java -cp ~/Lep_MAP3/bin OrderMarkers2 data=${out2} map=${out3} numThreads=1 chromosome=${lod} usePhysical=1 > ${out4}

# make orders look nice with physical positions 
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' ${out2}.snps.txt ${out4} | tail -n +4 | cut -f 1,2,3,4 > ${out4}.perf
