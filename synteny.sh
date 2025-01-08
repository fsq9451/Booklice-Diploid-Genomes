work_folder=/home/fsq2/Projects/52_lbdip_synteny
group1=1
haplotype1=A
group2=1
haplotype2=B

#Synteny jcvi
mkdir "$group1""$haplotype1""$group2""$haplotype2"
cd "$group1""$haplotype1""$group2""$haplotype2"
cp ../data/lbg"$group1"."$haplotype1".pep.fa ./lbg"$group1"_"$haplotype1".pep
cp ../data/lbg"$group2"."$haplotype2".pep.fa ./lbg"$group2"_"$haplotype2".pep
cp ../data/lbg"$group1"_dip.round3.all.onlymaker."$haplotype1".noseq.gff ./lbg"$group1"_"$haplotype1".gff3
cp ../data/lbg"$group2"_dip.round3.all.onlymaker."$haplotype2".noseq.gff ./lbg"$group2"_"$haplotype2".gff3

python -m jcvi.formats.gff bed --type=mRNA lbg"$group1"_"$haplotype1".gff3 -o lbg"$group1"_"$haplotype1".bed
python -m jcvi.formats.gff bed --type=mRNA lbg"$group2"_"$haplotype2".gff3 -o lbg"$group2"_"$haplotype2".bed

##共线性模块鉴定
python -m jcvi.compara.catalog ortholog lbg"$group1"_"$haplotype1" lbg"$group2"_"$haplotype2" --no_strip_names --dbtype=prot

## 点图
python -m jcvi.graphics.dotplot lbg"$group1"_"$haplotype1".lbg"$group2"_"$haplotype2".anchors

## 共线性图
#得到所有染色体名称，并以,链接
awk '{print $1}' ../data/lbg"$group1"."$haplotype1".length | sort | sed ":a;N;s/\n/,/g;ta" > lbg"$group1"_"$haplotype1".names
awk '{print $1}' ../data/lbg"$group2"."$haplotype2".length | sort | sed ":a;N;s/\n/,/g;ta" > lbg"$group2"_"$haplotype2".names

cat lbg"$group1"_"$haplotype1".names lbg"$group2"_"$haplotype2".names >seqids2.lbg"$group1"_"$haplotype1".lbg"$group2"_"$haplotype2"
sed 's/XXA/lbg"$group1"_"$haplotype1"/g; s/XXB/lbg"$group2"_"$haplotype2"/g' /pub/database/jcvi/layout2 > layout2.lbg"$group1"_"$haplotype1".lbg"$group2"_"$haplotype2"

python -m jcvi.compara.synteny screen --minspan=30 --simple lbg"$group1"_"$haplotype1".lbg"$group2"_"$haplotype2".anchors lbg"$group1"_"$haplotype1".lbg"$group2"_"$haplotype2".anchors.new
python -m jcvi.graphics.karyotype seqids2.lbg"$group1"_"$haplotype1".lbg"$group2"_"$haplotype2" layout2.lbg"$group1"_"$haplotype1".lbg"$group2"_"$haplotype2"
mv karyotype.pdf karyotype.lbg"$group1"_"$haplotype1".lbg"$group2"_"$haplotype2".pdf 