# mkdir group1 group2 group3
mkdir group1 group2

cd group1
pbsim   ../../transcriptome_for_simulation_3.fa \
        --prefix Sim_G1_S \
        --data-type CLR \
        --model_qc ../../pbsim-1.0.3-Linux-amd64/data/model_qc_clr \
        --length-mean 7800 \
        --length-sd 2300 \
        --length-min 100 \
        --length-max 20000 \
        --difference-ratio 33:36:31 \
        --accuracy-mean 0.87 \
        --accuracy-min 0.8 \
        --depth 4

cd ../group2/
pbsim   ../../transcriptome_for_simulation_3.fa \
        --prefix Sim_G2_S \
        --data-type CLR \
        --model_qc ../../pbsim-1.0.3-Linux-amd64/data/model_qc_clr \
        --length-mean 7800 \
        --length-sd 2300 \
        --length-min 100 \
        --length-max 20000 \
        --difference-ratio 33:36:31 \
        --accuracy-mean 0.87 \
        --accuracy-min 0.8 \
        --depth 10

# cd ../group3/
# pbsim   ../../transcriptome_for_simulation_3.fa \
#         --prefix Sim_G3_S \
#         --data-type CLR \
#         --model_qc ../../../pbsim-1.0.3-Linux-amd64/data/model_qc_clr \
#         --length-mean 7800 \
#         --length-sd 2300 \
#         --length-min 100 \
#         --length-max 20000 \
#         --difference-ratio 33:36:31 \
#         --accuracy-mean 0.87 \
#         --accuracy-min 0.8 \
#         --depth 30

cd ..
cat group1/*.fastq > dataset_sim_dm_ONT2D_g1.fastq
cat group2/*.fastq > dataset_sim_dm_ONT2D_g2.fastq
# cat group3/*.fastq > dataset_sim_dm_ONT2D_g3.fastq

python ../tran_qname.py dataset_sim_dm_ONT2D_g1.fastq SimG1_S g1.fastq
mv g1.fastq dataset_sim_dm_ONT2D_g1.fastq
python ../tran_qname.py dataset_sim_dm_ONT2D_g2.fastq SimG2_S g2.fastq
mv g2.fastq dataset_sim_dm_ONT2D_g2.fastq
# python ../tran_qname.py dataset_sim_dm_ONT2D_g3.fastq SimG3_S g3.fastq
# mv g3.fastq dataset_sim_dm_ONT2D_g3.fastq

rm group1/*.fastq
rm group2/*.fastq
# rm group3/*.fastq

# cat dataset_sim_dm_ONT2D_g1.fastq dataset_sim_dm_ONT2D_g2.fastq dataset_sim_dm_ONT2D_g3.fastq > dataset_sim_dm_ONT2D.fastq
python ../fastq_to_fasta.py dataset_sim_dm_ONT2D_g1.fastq
python ../fastq_to_fasta.py dataset_sim_dm_ONT2D_g2.fastq