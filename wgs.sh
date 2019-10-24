#! /bin/bash


###chmod +x TP2.sh to make it executable
###to be executed in tp_2 folder
###warning: change fastq/databases/softs paths if executed elsewhere


###################################
####analyse de l echantillon D ####
###################################


####bowtie2#####
ref_bwt2=databases/all_genome.fasta

soft/bowtie2-build $ref_bwt2 ${ref_bwt2%.*}

R1=fastq/EchD_1.fastq.gz
R2=fastq/EchD_2.fastq.gz
sam=${R1%_1.fastq.gz}.sam
soft/bowtie2 --local -x ${ref_bwt2%.*} -p 6 -1 $R1 -2 $R2 -S $sam  ###mapping


####samtools####
samtools=soft/samtools-1.6/samtools

unsorted=${R1%_1.fastq.gz}.unsorted.bam     
$samtools view -bt $ref_bwt2 -o $unsorted $sam          ###sam2bam    

bam=${R1%_1.fastq.gz}.sorted.bam                        ###bam sort
$samtools sort $unsorted -o $bam

$samtools index $bam                                    ###bam index

$samtools idxstats $bam > stats_EchD.txt

####stats#####

grep ">" databases/all_genome.fasta|cut -f 2 -d ">" >association.tsv ###get gi/subsp equivalence

cut -f1 stats_EchD.txt | while read line; do assoc=$(grep $line association.tsv); echo $assoc >> tmp.txt ; done
paste -d' ' stats_EchD.txt tmp.txt > EchD_stats.tsv ###fichier de base à retravailler sur libreoffice pour générer le fichier graph_EchD.png


####megahit####
megahit=soft/megahit

$megahit -1 $R1 -2 $R2 --k-min 21 --k-max 21 --kmin-1pass -o EchD.fasta

mv EchD.fasta/final.contigs.fa EchD.fasta/EchD.fasta

####Prodigal####
prodigal=soft/prodigal
prodigal -i EchD.fasta/EchD.fasta -d EchD_genes.fasta -s -p meta

sed "s:>:*\n>:g" EchD_genes.fasta | sed -n "/partial=00/,/*/p"|grep -v "*" > genes_full.fna


#####resfinder#####
ref_resfinder=databases/resfinder.fna
blastn=soft/blastn

makeblastdb -in $ref_resfinder -dbtype nucl

$blastn -db $ref_resfinder -query EchD_genes.fasta -evalue 0.001 -out EchD_genes_blastout.txt -qcov_hsp_perc 80 -perc_identity 80 -outfmt 6 -num_threads 6



####genes de resistances####
##tetracycline_tet(A)
##phenicol_floR
##beta-lactamase_blaCMY-2
##beta-lactamase_blaADC-25
##tetracycline_tet(M)
##beta-lactamase_blaCMY-76
##aminoglycoside_ant(6)-Ia
##beta-lactamase_blaKPC-3
##sulphonamide_sul1
##beta-lactamase_blaSHV-1
##macrolide_erm(B)
##beta-lactamase_blaSHV-11
##...









