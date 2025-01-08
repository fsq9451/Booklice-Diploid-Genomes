#1. extract lbg1_A ID
cd /home/fsq2/Projects/53_lbgdip_haplotype_orthofinder/data
for a in $(seq 1 8); do
grep '>' LB_BJ_A_$a.fasta | sed 's/>//' > LB_BJ_A_$a.names
done

#2. extract 1 vs 1 (lbg1A lbg1B lbg2A lbg2B lbg3A lbg3B LbruA)
cd /home/fsq2/Projects/54_lbdip_orthologs
cp /home/fsq2/Projects/52_lbdip_synteny/*/lbg1_A*.anchors.new ./
cp /home/fsq2/Projects/53_lbgdip_haplotype_orthofinder/data/*names ./

for a in *new; do
cat LB_BJ_A_1.names | xargs -i echo "grep  {} $a >> $a.chr1.orthologs " | sh
cat LB_BJ_A_2.names | xargs -i echo "grep  {} $a >> $a.chr2.orthologs " | sh
cat LB_BJ_A_3.names | xargs -i echo "grep  {} $a >> $a.chr3.orthologs " | sh
cat LB_BJ_A_4.names | xargs -i echo "grep  {} $a >> $a.chr4.orthologs " | sh
cat LB_BJ_A_5.names | xargs -i echo "grep  {} $a >> $a.chr5.orthologs " | sh
cat LB_BJ_A_6.names | xargs -i echo "grep  {} $a >> $a.chr6.orthologs " | sh
cat LB_BJ_A_7.names | xargs -i echo "grep  {} $a >> $a.chr7.orthologs " | sh
cat LB_BJ_A_8.names | xargs -i echo "grep  {} $a >> $a.chr8.orthologs " | sh
done

for a in *new; do
printf "LB_BJ_A\t$a\tscore\n" |sed 's/lbg1_A.//;s/.anchors.new//' > name.tmp
cat name.tmp $a.chr1.orthologs > $a.chr1.title.orthologs
cat name.tmp $a.chr2.orthologs > $a.chr2.title.orthologs
cat name.tmp $a.chr3.orthologs > $a.chr3.title.orthologs
cat name.tmp $a.chr4.orthologs > $a.chr4.title.orthologs
cat name.tmp $a.chr5.orthologs > $a.chr5.title.orthologs
cat name.tmp $a.chr6.orthologs > $a.chr6.title.orthologs
cat name.tmp $a.chr7.orthologs > $a.chr7.title.orthologs
cat name.tmp $a.chr8.orthologs > $a.chr8.title.orthologs
done


for a in $(seq 1 8); do
join -j 1 -o 1.1,1.2,2.2 <(sort -k1 lbg1_A.lbg1_B.anchors.new.chr$a.orthologs) <(sort -k1 lbg1_A.lbg2_A.anchors.new.chr$a.orthologs) > tmp1
join -j 1 -o 1.1,1.2,1.3,2.2 tmp1 <(sort -k1 lbg1_A.lbg2_B.anchors.new.chr$a.orthologs) > tmp2
join -j 1 -o 1.1,1.2,1.3,1.4,2.2 tmp2 <(sort -k1 lbg1_A.lbg3_A.anchors.new.chr$a.orthologs) > tmp3
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,2.2 tmp3 <(sort -k1 lbg1_A.lbg3_B.anchors.new.chr$a.orthologs) > tmp4
join -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.2 tmp4 <(sort -k1 lbg1_A.lbgLbru_A.anchors.new.chr$a.orthologs) > chr$a.all.orthologs
done


#3.extract data for each chrom
cat /home/fsq2/Projects/52_lbdip_synteny/data/*pep.fa /home/fsq2/Projects/52_lbdip_synteny/1ALbruA/lbgLbru_A.pep > all.pep
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\n":$0 }' all.pep | awk '{print $1}'  > all.pep.fasta
rm  all.pep

cat /home/fsq2/Projects/52_lbdip_synteny/data/*cds.fa /home/fsq2/Projects/37_hic_lb_annotation/Lbru/Lbru/3.maker/4.final-annotation/Lbru.rna > all.rna
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\n":$0 }' all.rna | awk '{print $1}'  > all.rna.fasta
rm  all.rna

for i in *.all.orthologs; do
mkdir /home/fsq2/Projects/54_lbdip_orthologs/$i.pepcds
cd /home/fsq2/Projects/54_lbdip_orthologs/$i.pepcds
cat ../$i | awk '{print $1}' > lbg1.name
cat lbg1.name | xargs -i echo "grep {} ../$i > {}.list" | sh
sed -i 's/ /\n/g' *list
ls -l *list | awk '$5!=105 {print}' | awk '{print $9}' | xargs -i echo 'rm {}' | sh
ls -l *list | wc -l > $i.number
ls *list | sed 's/.list//g' > gene.list
ls *.list | xargs -i echo "bioawk -cfastx 'BEGIN{while((getline k <\"{}\")>0)i[k]=1}{if(i[\$name])print \">\"\$name\"\n\"\$seq}' ../all.pep.fasta > {}.pep" |sh &
ls *.list | xargs -i echo "bioawk -cfastx 'BEGIN{while((getline k <\"{}\")>0)i[k]=1}{if(i[\$name])print \">\"\$name\"\n\"\$seq}' ../all.rna.fasta > {}.rna" |sh &
done

cd /home/fsq2/Projects/54_lbdip_orthologs
cat *orthologs.pep/*number > chrall.orthologs.number

#4.Rename all genes
cd /home/fsq2/Projects/52_lbdip_synteny/data
grep '>' lbg1.A.pep.fa | sed 's/>//g' > 1A.name
grep '>' lbg1.B.pep.fa | sed 's/>//g' > 1B.name
grep '>' lbg2.A.pep.fa | sed 's/>//g' > 2A.name
grep '>' lbg2.B.pep.fa | sed 's/>//g' > 2B.name
grep '>' lbg3.A.pep.fa | sed 's/>//g' > 3A.name
grep '>' lbg3.B.pep.fa | sed 's/>//g' > 3B.name
grep '>' Lbru.rna | sed 's/>//g' > Lbru.name
cp *name /home/fsq2/Projects/54_lbdip_orthologs

cd /home/fsq2/Projects/54_lbdip_orthologs/chr1.all.orthologs.pepcds 
for a in $(cat gene.list); do
awk '/^>/{print ">name" ++i; next}{print}' $a.list.rna > $a.list.rna.rename
awk '/^>/{print ">name" ++i; next}{print}' $a.list.pep > $a.list.pep.rename
done
cd /home/fsq2/Projects/54_lbdip_orthologs/chr2.all.orthologs.pepcds 
for a in $(cat gene.list); do
awk '/^>/{print ">name" ++i; next}{print}' $a.list.rna > $a.list.rna.rename
awk '/^>/{print ">name" ++i; next}{print}' $a.list.pep > $a.list.pep.rename
done
cd /home/fsq2/Projects/54_lbdip_orthologs/chr3.all.orthologs.pepcds 
for a in $(cat gene.list); do
awk '/^>/{print ">name" ++i; next}{print}' $a.list.rna > $a.list.rna.rename
awk '/^>/{print ">name" ++i; next}{print}' $a.list.pep > $a.list.pep.rename
done
cd /home/fsq2/Projects/54_lbdip_orthologs/chr4.all.orthologs.pepcds 
for a in $(cat gene.list); do
awk '/^>/{print ">name" ++i; next}{print}' $a.list.rna > $a.list.rna.rename
awk '/^>/{print ">name" ++i; next}{print}' $a.list.pep > $a.list.pep.rename
done
cd /home/fsq2/Projects/54_lbdip_orthologs/chr5.all.orthologs.pepcds 
for a in $(cat gene.list); do
awk '/^>/{print ">name" ++i; next}{print}' $a.list.rna > $a.list.rna.rename
awk '/^>/{print ">name" ++i; next}{print}' $a.list.pep > $a.list.pep.rename
done
cd /home/fsq2/Projects/54_lbdip_orthologs/chr6.all.orthologs.pepcds 
for a in $(cat gene.list); do
awk '/^>/{print ">name" ++i; next}{print}' $a.list.rna > $a.list.rna.rename
awk '/^>/{print ">name" ++i; next}{print}' $a.list.pep > $a.list.pep.rename
done
cd /home/fsq2/Projects/54_lbdip_orthologs/chr7.all.orthologs.pepcds 
for a in $(cat gene.list); do
awk '/^>/{print ">name" ++i; next}{print}' $a.list.rna > $a.list.rna.rename
awk '/^>/{print ">name" ++i; next}{print}' $a.list.pep > $a.list.pep.rename
done
cd /home/fsq2/Projects/54_lbdip_orthologs/chr8.all.orthologs.pepcds 
for a in $(cat gene.list); do
awk '/^>/{print ">name" ++i; next}{print}' $a.list.rna > $a.list.rna.rename
awk '/^>/{print ">name" ++i; next}{print}' $a.list.pep > $a.list.pep.rename
done

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/chr$a.all.orthologs.pepcds
sed -i 's/name1/LB1A/' *
sed -i 's/name2/LB1B/' *
sed -i 's/name3/LB2A/' *
sed -i 's/name4/LB2B/' *
sed -i 's/name5/LB3A/' *
sed -i 's/name6/LB3B/' *
sed -i 's/name7/Lbru/' *
done

#5.msa
cd /home/fsq2/Projects/54_lbdip_orthologs/
mkdir msa
cd msa

for a in $(seq 1 8); do
mkdir chr$a
cat /home/fsq2/Projects/54_lbdip_orthologs/chr$a.all.orthologs.pepcds/gene.list | xargs -i echo "perl /home/fsq2/bin/translatorx_vLocal.pl -i /home/fsq2/Projects/54_lbdip_orthologs/chr$a.all.orthologs.pepcds/{}.list.rna.rename -o chr$a/{} -c 1 -p M" > chr$a.sh
sh chr$a.sh &
done

#6. chromosome_trees
cd /home/fsq2/Projects/54_lbdip_orthologs/
mkdir trees
cd trees

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/trees
seqkit concat /home/fsq2/Projects/54_lbdip_orthologs/msa/chr$a/*.aa_ali.fasta > chr$a.aa.fasta
treebest nj chr$a.aa.fasta -W -t jtt -b 1000 > chr$a.treebest.out && sed 's/:\([0-9.]\+\)\[&&NHX:B=\([0-9]\+\)\]/\2:\1/' chr$a.treebest.out   | awk '{printf $0}' > chr$a.treebest.nwk &
done

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/trees
iqtree -s chr$a.aa.fasta -m MFP -bb 1000 -nt 30 -pre chr$a.iqtree &
done

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/trees
seqkit concat /home/fsq2/Projects/54_lbdip_orthologs/msa/chr$a/*.nt_ali.fasta > chr$a.nt.fasta
iqtree -s chr$a.nt.fasta -m MFP -bb 1000 -nt 64 -pre chr$a.nt.iqtree
done

#7. all_trees
#aa trees
cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees
mkdir chr$a
cd chr$a
cat /home/fsq2/Projects/54_lbdip_orthologs/chr$a.all.orthologs.pepcds/gene.list | xargs -i echo "iqtree -s /home/fsq2/Projects/54_lbdip_orthologs/msa/chr$a/{}.aa_ali.fasta -m VT -bb 1000 -nt 10 -pre {}.iqtree" > run$a.sh
done

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees/chr$a
nohup sh run*.sh &
done

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees
cd chr$a
cat /home/fsq2/Projects/54_lbdip_orthologs/chr$a.all.orthologs.pepcds/gene.list | xargs -i echo "pxrr -t {}.iqtree.contree -g Lbru -r > {}.iqtree.contree.reroot" | sh
cat /home/fsq2/Projects/54_lbdip_orthologs/chr$a.all.orthologs.pepcds/gene.list | xargs -i echo "java -jar /pub/software/paretree/PareTree1.0.2.jar -t O -nbs -topo -f {}.iqtree.contree.reroot &" | sh
done

#produce    *_nbs.iqtree.contree.reroot

#nt trees
cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees
mkdir nt.chr$a
cd nt.chr$a
cat /home/fsq2/Projects/54_lbdip_orthologs/chr$a.all.orthologs.pepcds/gene.list | xargs -i echo "iqtree -s /home/fsq2/Projects/54_lbdip_orthologs/msa/chr$a/{}.nt_ali.fasta -m GTR -bb 1000 -nt 10 -pre {}.iqtree" > run$a.nt.sh
done

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees/nt.chr$a
nohup sh run*.nt.sh &
done

#8. Summarize trees
for a in $(seq 1 8); do
mkdir /home/fsq2/Projects/54_lbdip_orthologs/densitrees/summary.chr$a
cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees/summary.chr$a

grep ',' ../chr$a/*_nbs.iqtree.contree.reroot | sed 's/...chr..//g;s/_nbs.iqtree.contree.reroot:/\t/g' > chr$a.all_nbs.iqtree.contree.reroot
awk '{print $2}' chr$a.all_nbs.iqtree.contree.reroot > chr$a.trees

grep ',' chr$a.all_nbs.iqtree.contree.reroot | grep 'LB1A,LB1B' >tmp1
grep ',' chr$a.all_nbs.iqtree.contree.reroot | grep 'LB1B,LB1A' >>tmp1

grep ',' tmp1 |  grep 'LB2A,LB2B' >tmp2
grep ',' tmp1 |  grep 'LB2B,LB2A' >>tmp2

grep ',' tmp2 |  grep 'LB3A,LB3B' >tmp3
grep ',' tmp2 |  grep 'LB3B,LB3A' >>tmp3

cp tmp3 contrees.final.txt
awk '{print $1}' contrees.final.txt > chr$a.same.gene.list
done

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees/summary.chr$a
printf "chr$a\n" > numbers.txt
cat contrees.final.txt | wc -l >> numbers.txt
cat chr$a.trees | wc -l >> numbers.txt
cat /home/fsq2/Projects/54_lbdip_orthologs/chr$a.all.orthologs.pepcds/gene.list | wc -l >>numbers.txt
done

cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees/
mkdir summary.all
cd summary.all
paste ../summary*/numbers.txt > summary.all.txt
cat ../summary.chr*/chr*.same.gene.list > all.same.gene.list

#9. sumterees
for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees/summary.chr$a
sumtrees.py chr$a.trees -o chr$a.log -x chr$a
done

cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees/summary.all
cat /home/fsq2/Projects/54_lbdip_orthologs/densitrees/summary.chr*/chr?.trees > all.trees
sumtrees.py all.trees -o all.trees.log -x all.trees

#10. ASTRAL tree
cd /home/fsq2/Projects/54_lbdip_orthologs/densitrees/
mkdir astral.aa
cd astral.aa

for a in $(seq 1 8); do
cat ../chr$a/*contree > ./aa.chr$a.contree
java -jar /pub/software/Astral/astral.5.7.8.jar -i ./aa.chr$a.contree -o aa.chr$a.out.tre 2>aa.chr$a.out.log
done


#11. dN/dS

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/selection
mkdir chr$a
cd chr$a
cat /home/fsq2/Projects/54_lbdip_orthologs/chr$a.all.orthologs.pepcds/gene.list | xargs -i echo "grep -A1 LB1A /home/fsq2/Projects/54_lbdip_orthologs/msa/chr$a/{}.nt_ali.fasta > {}.msa.fasta" > run_1.sh
cat /home/fsq2/Projects/54_lbdip_orthologs/chr$a.all.orthologs.pepcds/gene.list | xargs -i echo "grep -A1 LB1B /home/fsq2/Projects/54_lbdip_orthologs/msa/chr$a/{}.nt_ali.fasta >> {}.msa.fasta" > run_2.sh
sh run_1.sh
sh run_2.sh
done

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/selection
cd chr$a
cat /home/fsq2/Projects/54_lbdip_orthologs/chr$a.all.orthologs.pepcds/gene.list | xargs -i echo "sed \"s/genename/{}/g\" /home/fsq2/Projects/54_lbdip_orthologs/selection/Nsite_0.ctl > {}.ctl" > run_codeml_1.sh
cat /home/fsq2/Projects/54_lbdip_orthologs/chr$a.all.orthologs.pepcds/gene.list | xargs -i echo "codeml {}.ctl" > run_codeml_2.sh
sh run_codeml_1.sh
done

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/selection
cd chr$a
sh run_codeml_2.sh
done

for a in $(seq 1 8); do
cd /home/fsq2/Projects/54_lbdip_orthologs/selection
grep 'omega (dN/dS)' chr$a/*mlc_0 | sed 's/\/Lbbj/\t\/Lbbj/g;s/.mlc_0:omega (dN\/dS) =  /\t/g' > summary/selection_chr$a.txt
grep -v omega summary/selection_chr$a.txt > summary/selection_chr$a.final.txt
grep -v omega summary/selection_chr$a.txt | awk '$3<2{print}' > summary/selection_chr$a.final.0-2.txt
done

cd summary
cat selection_chr*.final.txt > selection_all.final.txt
cat selection_chr*.final.0-2.txt > selection_all.final.0-2.txt

awk '$3>0.6 {print $2}' selection_all.final.txt | sed 's/\///g'  > 0.6higher.txt










