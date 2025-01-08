compare_pair=lbg1A_lbg1B
species_prefix_1=LB_BJ_A
species_prefix_2=LB_BJ_B
work_dir=/home/fsq2/Projects/49_mummer_lb/new_folder
lbg1_hic_fasta=/home/fsq2/Projects/37.1_lbg123_annotation_new/lbg1_dip/lbg1_dip/2.RepeatMasker/rearrange/lbg1.AB.fasta.masked
lbg2_hic_fasta=/home/fsq2/Projects/37.1_lbg123_annotation_new/lbg2_dip/lbg2_dip/2.RepeatMasker/rearrange/lbg2.AB.fasta.masked
lbg3_hic_fasta=/home/fsq2/Projects/37.1_lbg123_annotation_new/lbg3_dip/lbg3_dip/2.RepeatMasker/lbg3.AB.fasta.masked
l=20
c=65
#i=50
threads=60

<<prepare
################data prepare################
#copy
cd $work_dir
mkdir data
cd data
cat $lbg1_hic_fasta $lbg2_hic_fasta $lbg3_hic_fasta > lbg123.AB.fasta

#
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\n":$0 }' lbg123.AB.fasta  > tmp
mv tmp lbg123.AB.fasta

#name list
grep '>' lbg123.AB.fasta | sed 's/>//g' > names

#extract
for id in $(cat names); do
grep -A1 $id lbg123.AB.fasta > $id.fasta
done
################data prepare################
prepare

<<eof

################nucmer alignment################
#激活环境 fsq2下的mummer4
#conda activate mummer4

###Basic mummer
for a in $(seq 1 8); do
cd $work_dir/new_l20c65
mkdir $compare_pair
cd $compare_pair
#mummer
nucmer -l $l -c $c -p chr$a -t $threads $work_dir/data/"$species_prefix_1"_$a.fasta $work_dir/data/"$species_prefix_2"_$a.fasta
#delta-filter
delta-filter -l $l -1 chr$a.delta > chr$a.delta.filter
#plot
/pub/miniconda3/bin/mummerplot/mummerplot -p chr$a.delta.filter chr$a.delta.filter -t png
#show-coords
show-coords -r -c -l chr$a.delta.filter > chr$a.delta.filter.coord
LINKVIEW.py -t2 -o chr$a chr$a.delta.filter.coord
#dnadiff
dnadiff -d chr$a.delta.filter -p chr$a.delta.filter.diff
# $a.diff.snps contains the SNVs and INDELs
mv chr$a.delta.filter.diff.snps chr$a.delta.filter.diff.var
# split variations into two files (snp and indel).
python3 /home/fsq2/software/marmoset/marmoset-master/Mummer_alignment/split_snp_indel.py chr$a.delta.filter.diff.var chr$a.delta.filter.diff
# change indel format for easy reading
python3 /home/fsq2/software/marmoset/marmoset-master/Mummer_alignment/indel_statistic_RefPos.py chr$a.delta.filter.diff.indel chr$a.delta.filter.diff
# split small and longer indels
awk '$5<=50' chr$a.delta.filter.diff.indel.report.txt >chr$a.delta.filter.diff.sm.indel
awk '$5>50' chr$a.delta.filter.diff.indel.report.txt >chr$a.delta.filter.diff.lg.indel
sed "s/.delta.filter.diff//g;s/chr/"$species_prefix_1"_/g" chr$a.delta.filter.diff.sm.indel | grep "$species_prefix_1" | awk '{print $1"\t"$2"\t"$2}' > chr$a.delta.filter.diff.sm.indel.bed
sed "s/.delta.filter.diff//g;s/chr/"$species_prefix_1"_/g" chr$a.delta.filter.diff.lg.indel | grep "$species_prefix_1" | awk '{print $1"\t"$2"\t"$2}' > chr$a.delta.filter.diff.lg.indel.bed
done

###Intersect 50k
cd $work_dir/new_l20c65/$compare_pair
#生成窗口文件， 窗口大小50Kb
seqtk comp $lbg1_hic_fasta |grep "^$species_prefix_1" |awk '{print $1"\t"$2}' > genome.len
bedtools  makewindows -w 50000 -g  genome.len > genome.window.bed
cd $work_dir/new_l20c65/$compare_pair

#calculation
for a in $(seq 1 8); do
###SNP
awk '{print $11"\t"$1"\t"$1}' chr$a.delta.filter.diff.snp >> chr$a.delta.filter.diff.snp.bed
###INDEL
#All indel bed
sed "s/.delta.filter.diff//g;s/chr/"$species_prefix_1"_/g" chr$a.delta.filter.diff.indel.report.txt | grep "$species_prefix_1" > chr$a.delta.filter.diff.indel.modified.report.txt
awk '{print $1"\t"$2"\t"$2}' chr$a.delta.filter.diff.indel.modified.report.txt > chr$a.delta.filter.diff.indel.modified.bed
done

cd $work_dir/new_l20c65/$compare_pair
#SNP
cat chr*.delta.filter.diff.snp.bed > All.delta.filter.diff.snp.bed
bedtools intersect -a genome.window.bed -b All.delta.filter.diff.snp.bed -c > All.l$l.c$c.snp.txt
#indel
#all
cat chr*.delta.filter.diff.indel.modified.bed > All.delta.filter.diff.indel.modified.bed
bedtools intersect -a genome.window.bed -b All.delta.filter.diff.indel.modified.bed -c > All.l$l.c$c.indel.txt
#sm,lg
cat chr*.delta.filter.diff.sm.indel.bed > All.delta.filter.diff.sm.indel.bed
cat chr*.delta.filter.diff.lg.indel.bed > All.delta.filter.diff.lg.indel.bed
bedtools intersect -a genome.window.bed -b All.delta.filter.diff.sm.indel.bed -c > All.l$l.c$c.sm.indel.txt
bedtools intersect -a genome.window.bed -b All.delta.filter.diff.lg.indel.bed -c > All.l$l.c$c.lg.indel.txt


##plots
cd $work_dir/new_l20c65/$compare_pair
mkdir plots
cd plots
cp ../All.l$l.c$c.snp.txt ./
cp ../All.l$l.c$c.indel.txt ./
cp ../All.l$l.c$c.sm.indel.txt ./
cp ../All.l$l.c$c.lg.indel.txt ./

for a in $(seq 1 8); do
grep "_$a" All.l$l.c$c.snp.txt > chr$a.snp.txt
grep "_$a" All.l$l.c$c.indel.txt > chr$a.indel.txt
grep "_$a" All.l$l.c$c.sm.indel.txt > chr$a.sm.indel.txt
grep "_$a" All.l$l.c$c.lg.indel.txt > chr$a.lg.indel.txt
done


cd $work_dir/new_l20c65/$compare_pair/plots
for a in $(seq 1 8); do
printf "type\tlength\tnumber\n" > name.tmp
awk '{print $1"\t"$2"\t"$4}' chr$a.snp.txt | sed "s/$species_prefix_1../snp/g" > snp.tmp
awk '{print $1"\t"$2"\t"$4}' chr$a.indel.txt | sed "s/$species_prefix_1../indel/g" > indel.tmp
awk '{print $1"\t"$2"\t"$4}' chr$a.sm.indel.txt | sed "s/$species_prefix_1../small_indel/g" > sm.indel.tmp
awk '{print $1"\t"$2"\t"$4}' chr$a.lg.indel.txt | sed "s/$species_prefix_1../large_indel/g" > lg.indel.tmp

cat name.tmp snp.tmp indel.tmp sm.indel.tmp lg.indel.tmp > chr$a.4kinds.txt
cat name.tmp snp.tmp indel.tmp > chr$a.snp_indel.txt
done



cd $work_dir/new_l20c65/$compare_pair/plots
rename "s/chr/$compare_pair.chr/" *

for a in $(seq 1 8); do
awk '{print $1}' $compare_pair.chr$a.snp_indel.txt | sed "s/snp/$compare_pair/g;s/indel/$compare_pair/g;s/type/compare/g" > 1list.tmp
paste 1list.tmp $compare_pair.chr$a.snp_indel.txt > 4line.$compare_pair.chr$a.snp_indel.txt
done

cd $work_dir/new_l20c65/$compare_pair/plots
for a in $(seq 1 8); do
awk '{print $4}' $compare_pair.chr$a.snp.txt  > $compare_pair.chr$a.snp.number
awk '{print $4}' $compare_pair.chr$a.indel.txt  > $compare_pair.chr$a.indel.number
paste $compare_pair.chr$a.snp.number $compare_pair.chr$a.indel.number > $compare_pair.chr$a.snp.indel.number
done

printf "snp\tindel\n" > name.tmp
cat name.tmp *.snp.indel.number > $compare_pair.All.snp.indel.number.txt


cd $work_dir/new_l20c65/$compare_pair/plots

for a in $(seq 1 8); do
printf "type\tlength\tnumber\n" > name.tmp
awk '{print $1"\t"$4}' $compare_pair.chr$a.snp.txt | sed "s/$species_prefix_1../snp/g" >> $compare_pair.snp.tmp
awk '{print $1"\t"$4}' $compare_pair.chr$a.indel.txt | sed "s/$species_prefix_1../indel/g" >> $compare_pair.indel.tmp
done


cat $compare_pair.snp.tmp $compare_pair.indel.tmp > $compare_pair.snp.indel.tmp
printf "type\tnumber\n" > name.tmp
cat name.tmp $compare_pair.snp.indel.tmp > $compare_pair.snp.indel.txt

eof
###vcf
cd $work_dir/new_l20c65/$compare_pair

for a in $(seq 1 8); do
        all2vcf mummer --snps chr$a.delta.filter.diff.indel --reference /home/fsq2/Projects/49_mummer_lb/data/LB_BJ_A_$a.fasta > LB_BJ_AB_$a.indel.vcf
        all2vcf mummer --snps chr$a.delta.filter.diff.snp --reference /home/fsq2/Projects/49_mummer_lb/data/LB_BJ_A_$a.fasta > LB_BJ_AB_$a.snp.vcf
done    

cat LB_BJ_AB_*.indel.vcf > LB_BJ_AB_all.indel.vcf
cat LB_BJ_AB_*.snp.vcf > LB_BJ_AB_all.snp.vcf

cd $work_dir/new_l20c65/$compare_pair
mkdir snp_indel_snpeff
cd snp_indel_snpeff

cp /home/fsq2/Projects/49_mummer_lb/old_2/l20c65_no_i/lbg1A_lbg1B/snp_indel_snpeff/snpEff.config ./
cp -rf /home/fsq2/Projects/49_mummer_lb/old_2/l20c65_no_i/lbg1A_lbg1B/snp_indel_snpeff/data ./

java -Xmx4g -jar /pub/miniconda3/pkgs/snpeff-5.0-hdfd78af_1/share/snpeff-5.0-1/snpEff.jar -c snpEff.config -ud 1500 -csvStats all.indel.csv -htmlStats all.indel.html -o vcf lbg1_dip ../LB_BJ_AB_all.indel.vcf   > $compare_pair.final.indel.vcf

java -Xmx4g -jar /pub/miniconda3/pkgs/snpeff-5.0-hdfd78af_1/share/snpeff-5.0-1/snpEff.jar -c snpEff.config -ud 1500 -csvStats all.snp.csv -htmlStats all.snp.html -o vcf lbg1_dip ../LB_BJ_AB_all.snp.vcf   > $compare_pair.final.snp.vcf

<<eof

#calibrated with the sequences length

cd $work_dir/new_l20c65/$compare_pair
for a in $(seq 1 8); do
grep '| LB' chr$a.delta.filter.coord  | awk '{print $18"\t"$1"\t"$2}' > chr$a.coord.bed
bedtools intersect -a genome.window.bed -b chr$a.coord.bed -wo > chr$a.sequence.numbers.bed
done

printf "length\tuseful_length\n" > tmp
awk '$1=="LB_BJ_A_1"{ seen[$2] += $7 } END { for (i in seen) print i, seen[i] }' chr1.sequence.numbers.bed > chr1.numbers.tmp
awk '$1=="LB_BJ_A_2"{ seen[$2] += $7 } END { for (i in seen) print i, seen[i] }' chr2.sequence.numbers.bed > chr2.numbers.tmp
awk '$1=="LB_BJ_A_3"{ seen[$2] += $7 } END { for (i in seen) print i, seen[i] }' chr3.sequence.numbers.bed > chr3.numbers.tmp
awk '$1=="LB_BJ_A_4"{ seen[$2] += $7 } END { for (i in seen) print i, seen[i] }' chr4.sequence.numbers.bed > chr4.numbers.tmp
awk '$1=="LB_BJ_A_5"{ seen[$2] += $7 } END { for (i in seen) print i, seen[i] }' chr5.sequence.numbers.bed > chr5.numbers.tmp
awk '$1=="LB_BJ_A_6"{ seen[$2] += $7 } END { for (i in seen) print i, seen[i] }' chr6.sequence.numbers.bed > chr6.numbers.tmp
awk '$1=="LB_BJ_A_7"{ seen[$2] += $7 } END { for (i in seen) print i, seen[i] }' chr7.sequence.numbers.bed > chr7.numbers.tmp
awk '$1=="LB_BJ_A_8"{ seen[$2] += $7 } END { for (i in seen) print i, seen[i] }' chr8.sequence.numbers.bed > chr8.numbers.tmp

cat tmp chr1.numbers.tmp >"$compare_pair".chr1.number
cat tmp chr2.numbers.tmp >"$compare_pair".chr2.number
cat tmp chr3.numbers.tmp >"$compare_pair".chr3.number
cat tmp chr4.numbers.tmp >"$compare_pair".chr4.number
cat tmp chr5.numbers.tmp >"$compare_pair".chr5.number
cat tmp chr6.numbers.tmp >"$compare_pair".chr6.number
cat tmp chr7.numbers.tmp >"$compare_pair".chr7.number
cat tmp chr8.numbers.tmp >"$compare_pair".chr8.number
