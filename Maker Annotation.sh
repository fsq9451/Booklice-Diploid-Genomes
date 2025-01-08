#https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2

##MUST
work_folder=/home/fsq2/Projects/37.1_lbg123_annotation_new/lbg1_dip
species_name=lbg1_dip
genome=/home/fsq2/Projects/36.1_hic_lb_final/lbg1_dip/rearrange1/lbg1.AB.fasta
est=/home/fsq2/Projects/37_hic_lb_annotation/data/trinity-rna/trinity_G123.individual.Trinity.lbg1_Final.fasta
protein=/home/fsq2/Projects/37_hic_lb_annotation/data/protein/protein.annotation.lbg1_Final.fasta
repeat=$work_folder/$species_name/1.RepeatModeler/"$species_name"-families.fa
threads=64

##Maker parameters
#1 or 0
est2genome=1
protein2genome=1
trna=0
alt_splice=0
#gene id prefix
maker_prefix=Lbbj_
#gene id numbers
maker_number=6

<<eof
##创建工作目录
mkdir $work_folder/$species_name

#####1.RepeatModeler######
cd $work_folder/$species_name
mkdir 1.RepeatModeler
cd 1.RepeatModeler

BuildDatabase -name $species_name $genome
RepeatModeler -pa $threads -engine ncbi -database $species_name -LTRStruct


#####2.RepeatMasker######
cd $work_folder/$species_name
mkdir 2.RepeatMasker
cd 2.RepeatMasker

RepeatMasker -e ncbi -lib ../1.RepeatModeler/"$species_name"-families.fa -pa $threads -gff -dir ./ $genome


#####3.maker######
mkdir $work_folder/$species_name/3.maker


###1.round1
cd $work_folder/$species_name/3.maker
mkdir 1.round1
cd 1.round1
cp /pub/database/maker/maker* ./

echo "#-----Genome (these are always required)
genome=$genome  #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=$est #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=$protein  #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib=$repeat   #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/pub/software/maker3/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=$est2genome #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=$protein2genome #infer predictions from protein homology, 1 = yes, 0 = no
trna=$trna #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=10000 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=50 #require at least this many amino acids in predicted proteins
alt_splice=$alt_splice #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files" >maker_opts.ctl


#1.1round1 - run
/home/caupql2/mpich2/bin/mpiexec -n $threads maker -base $species_name.round1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

#Annotation merge,fasta merge (for mRNA, protein sequences), GFF w/o the sequences
cd $species_name.round1.maker.output
gff3_merge -s -d $species_name.round1_master_datastore_index.log > $species_name.round1.all.maker.gff
fasta_merge -d $species_name.round1_master_datastore_index.log
gff3_merge -s -n -d $species_name.round1_master_datastore_index.log > $species_name.round1.all.maker.noseq.gff


#1.2round1 - software training

##SNAP
cd $work_folder/$species_name/3.maker
mkdir -p snap/round1
cd snap/round1
# export 'confident' gene models from MAKER and rename to something meaningful
maker2zff -c 0.8 -e 0.8 -o 0.8 -x 0.2 -l 50 -d $work_folder/$species_name/3.maker/1.round1/$species_name.round1.maker.output/$species_name.round1_master_datastore_index.log
# gather some stats and validate
fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
fathom genome.ann genome.dna -validate > validate.log 2>&1
# collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
# create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
# assembly the HMM
hmm-assembler.pl $species_name.round1.zff.length50_aed0.2 params > $species_name.round1.zff.length50_aed0.2.hmm

##AUGUST
cd $work_folder/$species_name/3.maker
mkdir -p august/round1
cd august/round1
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' $work_folder/$species_name/3.maker/1.round1/$species_name.round1.maker.output/$species_name.round1.all.maker.noseq.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
  bedtools getfasta -fi $genome -bed - -fo $species_name.round1.all.maker.transcripts1000.fasta

#BUSCO augustus
busco -i $species_name.round1.all.maker.transcripts1000.fasta  -o $species_name.round1 -l /pub/database/busco/insecta_odb10/ \
  -m genome -c 30 --long --augustus_species fly --augustus_parameters='--progress=true' --offline

mkdir $species_name.augustus_retraining_parameters
cd $species_name.augustus_retraining_parameters
cp $work_folder/$species_name/3.maker/august/round1/$species_name.round1/run_/augustus_output/retraining_parameters/BUSCO_$species_name.round1/* ./

rename "s/BUSCO_$species_name.round1/$species_name.round1/g" *

sed -i "s/BUSCO_$species_name.round1/$species_name.round1/g" $species_name.round1_parameters.cfg
sed -i "s/BUSCO_$species_name.round1/$species_name.round1/g" $species_name.round1_parameters.cfg.orig1

# may need to sudo
mkdir $AUGUSTUS_CONFIG_PATH/species/$species_name.round1
cp $species_name.round1*  $AUGUSTUS_CONFIG_PATH/species/$species_name.round1


eof

###2.round2
cd $work_folder/$species_name/3.maker
mkdir 2.round2
cd 2.round2
cp /pub/database/maker/maker* ./

echo "#-----Genome (these are always required)
genome=$genome  #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=$est #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=$protein  #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib=$repeat   #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/pub/software/maker3/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=$work_folder/$species_name/3.maker/snap/round1/$species_name.round1.zff.length50_aed0.2.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=$species_name.round1 #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=$est2genome #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=$protein2genome #infer predictions from protein homology, 1 = yes, 0 = no
trna=$trna #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=10000 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=50 #require at least this many amino acids in predicted proteins
alt_splice=$alt_splice #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files" >maker_opts.ctl


#2.1round2 - run
/home/caupql2/mpich2/bin/mpiexec -n $threads maker -base $species_name.round2 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

#Annotation merge,fasta merge (for mRNA, protein sequences), GFF w/o the sequences
cd $species_name.round2.maker.output
gff3_merge -s -d $species_name.round2_master_datastore_index.log > $species_name.round2.all.maker.gff
fasta_merge -d $species_name.round2_master_datastore_index.log
gff3_merge -s -n -d $species_name.round2_master_datastore_index.log > $species_name.round2.all.maker.noseq.gff

#Extract only maker
awk '$2=="maker" {print $0}' $species_name.round2.all.maker.noseq.gff  > $species_name.round2.all.onlymaker.noseq.gff

#Evaluation
#A. Count the number of gene models and the gene lengths after each round.
cat $species_name.round2.all.maker.gff | awk '{ if ($3 == "gene") print $0 }' | awk '{ sum += ($5 - $4) } END { print NR, sum / NR }'  > gene.numbers.txt
#B. Visualize the AED distribution. AED ranges from 0 to 1 and quantifies the confidence in a gene model based on empirical evidence. Basically, the lower the AED, the better a gene model is likely to be. Ideally, 95% or more of the gene models will have an AED of 0.5 or better in the case of good assemblies. 
perl /pub/software/maker3/bin/AED_cdf_generator.pl -b 0.025 $species_name.round2.all.maker.gff  > AED.txt

#BUSCO evaluation
mkdir busco
cd busco
busco -i ../$species_name.round2.all.maker.transcripts.fasta  -o round2_trans -l /pub/database/busco/insecta_odb10/ -m transcriptome -c 30 --offline &
busco -i ../$species_name.round2.all.maker.proteins.fasta  -o round2_proteins -l /pub/database/busco/insecta_odb10/ -m protein -c 30 --offline
cd ..

#REname ---- final version
mkdir final_files
cd final_files
maker_map_ids --prefix $maker_prefix --justify $maker_number ../$species_name.round2.all.maker.gff > $species_name.round2.id.map

cp ../$species_name.round2.all.maker.transcripts.fasta ./
cp ../$species_name.round2.all.maker.proteins.fasta ./
cp ../$species_name.round2.all.maker.noseq.gff ./
cp ../$species_name.round2.all.onlymaker.noseq.gff ./

map_fasta_ids $species_name.round2.id.map $species_name.round2.all.maker.proteins.fasta
map_fasta_ids $species_name.round2.id.map $species_name.round2.all.maker.transcripts.fasta
map_gff_ids $species_name.round2.id.map l$species_name.round2.all.maker.noseq.gff
map_gff_ids $species_name.round2.id.map $species_name.round2.all.onlymaker.noseq.gff


#FINALLLLL
cd $work_folder/$species_name/3.maker
mkdir 4.final-annotation
cd 4.final-annotation
cp $work_folder/$species_name/3.maker/2.round2/$species_name.round2.maker.output/final_files/* ./


####DONE  !!!

