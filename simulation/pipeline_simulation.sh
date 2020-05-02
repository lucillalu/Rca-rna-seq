python Annotation_grouping.py Drosophila_melanogaster.BDGP6.94.gtf
python generate_transcriptome.py Drosophila_melanogaster.BDGP6.94_AS.gtf Drosophila_melanogaster.BDGP6.fa AS_transcriptome.fa
python generate_transcriptome.py Drosophila_melanogaster.BDGP6.94_SS.gtf Drosophila_melanogaster.BDGP6.fa SS_transcriptome.fa
python generate_transcriptome.py Drosophila_melanogaster.BDGP6.94_short.gtf Drosophila_melanogaster.BDGP6.fa short_transcriptome.fa
cat AS_transcriptome.fa SS_transcriptome.fa short_transcriptome.fa > merge_transcriptome.fa

python ../samscripts/fastqfilter.py minlen 200 merge_transcriptome.fa transcriptome_for_simulation.fa
python fa_to_3.py transcriptome_for_simulation.fa

rm Drosophila_melanogaster.BDGP6.94_AS.gtf Drosophila_melanogaster.BDGP6.94_SS.gtf Drosophila_melanogaster.BDGP6.94_short.gtf
rm AS_transcriptome.fa SS_transcriptome.fa short_transcriptome.fa merge_transcriptome.fa

mkdir ONT2D
cd ONT2D
bash ../simulate_ONT2D.sh
cd ..

mkdir ONT2D_1
cd ONT2D_1
bash ../simulate_ONT2D_1.sh
cd ..