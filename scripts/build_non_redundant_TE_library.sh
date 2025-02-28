#!/bin/sh
# Overview: Build TE library from each species and then cluster the repeat libraries.
# Step 1: Build TE library from each species. Use M.bellicosus as an example.
sp=Mbel
genome=Mbel_Q1_clean.fasta
db_name=${sp}_clean
~/bin/RepeatModeler-2.0.5/BuildDatabase -name $db_name $genome
~/bin/RepeatModeler-2.0.5/RepeatModeler -database $db_name -threads 64 -LTRStruct

# Step 2: Run MC_helper (https://github.com/GonzalezLab/MCHelper).
repeats=${sp}_clean-families.fa
hmm_db=~insecta_odb10.hmm
output=${sp}_MC
python3 MCHelper.py -l $repeats -o $output -g $genome --input_type fasta -b $hmm_db -a F -t 64

# Step 3: Cluster repeat libraries.
r_repeat_lib=termite_TE.fa
nr_repeat_lib=termite_TE_clustered_90.fa

cat *fasta > $r_repeat_lib # Gather all repeat libraries from the output of MChelper.

# Run ClusterPartialMatchingSubs.pl with the interval: 100 - 500, 500 - 1000, 1000 to 2000, 2000 to 5000, 5000 to 10,000, 10,000 to 20,000, and higher.
#E.g., interval 500 to 1000
~/bin/RepeatModeler-2.0.5/util/ClusterPartialMatchingSubs.pl $r_repeat_lib -lmax 1000 -lmin 500 -n 3 -t 64 -a

cat */${r_repeat_lib}*_90 > $nr_repeat_lib # Gather the output to make a non-redundant repeat library.

# Step 4: Run RepeatMasker to annotate genomes with the non-redundant repeat library.
out_dir=${sp}_RM_90
~/bin/RepeatMasker/RepeatMasker -pa 64 -gff -lib $nr_repeat_lib -dir $out_dir $genome -e rmblast -xsmall

