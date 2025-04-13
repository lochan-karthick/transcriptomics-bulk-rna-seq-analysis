#! /bin/bash
#PBS -l nodes=1:ppn=4:centos7,cput=24:00:00,walltime=48:00:00
#PBS -N RNAseq_Tuxedo_assignment
#PBS -d /export/biostuds/3046323r/transcript_assignment1
#PBS -m abe
#PBS -M 3046323R@student.gla.ac.uk
#PBS -q bioinf-stud
#PBS -o /export/biostuds/3046323r/transcript_assignment1/assignment_output.log
#PBS -e /export/biostuds/3046323r/transcript_assignment1/assignment_error.log

#ensuring bin is in path:
export PATH=$HOME/bin:$PATH

#Setting Directory:
Working_dir="/export/biostuds/3046323r/transcript_assignment1"
mkdir -p ${Working_dir}
cd ${Working_dir}

#Setting Paths to Resources:
Data_dir="/export/home/gmh5n/assignment"
Genome_dir="/export/projects/polyomics/Genome/Mus_musculus/mm10/Hisat2Index/chr2"
GTF_file="/export/projects/polyomics/Genome/Mus_musculus/mm10/annotations/chr2.gtf"
Adapter_file="/export/projects/polyomics/biostuds/data/illumina_adapter.fa"

#Creating Output Directories:
hisat_dir="${Working_dir}/hisat2_results"
stringtie_dir="${Working_dir}/stringtie"
fastq_dir="${Working_dir}/fastq"
trimmed_dir="${Working_dir}/trimmed"
count_res="${Working_dir}/count_matrices"
mkdir -p ${hisat_dir} ${stringtie_dir} ${fastq_dir} ${trimmed_dir} ${count_res}

#creating list for gtf:
gtflist="gtf_list.txt"
rm -f ${gtflist}

#list of samples:
samples=("s1.c2" "s2.c2" "s3.c2" "s4.c2" "s5.c2" "s6.c2" "s7.c2" "s8.c2" "s9.c2" "s10.c2" "s11.c2" "s12.c2")

#creating soft links:
for sample in "${samples[@]}"; do
    ln -sfn ${Data_dir}/${sample}.fq ${fastq_dir}/${sample}.fq
done

#extracting 25K reads:
#for sample in "${samples[@]}"; do
    #awk 'NR<=100000' ${sample} > fastq/${sample}.fq
#done

#One loop for QC, alignment, sorting and assembly for every sample:
for sample in "${samples[@]}"; do
    echo "Processing ${sample}..."
    
    fastq="${fastq_dir}/${sample}.fq"
    trim1="${trimmed_dir}/${sample}_t1.fq"
    trim2="${trimmed_dir}/${sample}_t2.fq"
    sam="${hisat_dir}/${sample}.sam"
    bam="${hisat_dir}/${sample}.bam"
    sorted_bam="${hisat_dir}/${sample}.sorted.bam"

    #Trimming adapters using scythe:
    scythe -a ${Adapter_file} -q sanger -o ${trim1} ${fastq}
    if [ $? -ne 0 ]; then
        echo "Scythe trimming failed for ${sample}. Exiting."
        exit 1
    fi

    #Quality trimming using sickle:
    sickle se -f ${trim1} -t sanger -o ${trim2} -q 10 -l 50
    if [ $? -ne 0 ]; then
        echo "Sickle low quality base trimming failed for ${sample}. Exiting."
        exit 1
    fi
 
    #Using HISAT 2 on stranded mode to align the reads:
    hisat2 -p 4 --rna-strandness RF --phred33 -x ${Genome_dir} -U ${trim2} -S ${sam} 
    if [ $? -ne 0 ]; then
        echo "HISAT2 alignment failed for ${sample}. Exiting."
        exit 1
    fi
  
    #Using SAMtools to convert SAM to BAM and sorting:
    samtools view -b -o ${bam} ${sam}
    if [ $? -ne 0 ]; then
        echo "SAMtools view command failed for ${sample}. Exiting."
        exit 1
    fi
    samtools sort -o ${sorted_bam} ${bam}
    if [ $? -ne 0 ]; then
        echo "SAMtools sorting on .bam failed for ${sample}. Exiting."
        exit 1
    fi

    #Cleaningup files:
    rm ${sam} ${bam} ${trim1} ${trim2}

    #Tying up the individual transcriptomes using stringtie:
    str_smp_dir="${stringtie_dir}/${sample}"
    mkdir -p ${str_smp_dir}
    sample_tr_gtf="${str_smp_dir}/${sample}_transcripts.gtf"

    stringtie -p 4 -e -B -t --rf -G ${GTF_file} -o ${sample_tr_gtf} ${sorted_bam} 
    if [ $? -ne 0 ]; then
        echo "StringTie failed for ${sample}. Exiting."
        exit 1
    fi

    #Appending necessary Addresses to GTF list:
    echo -e "${sample}\t${sample_tr_gtf}" >> ${gtflist}
done

#Running prepDE.py to obtain count tables:
python ${HOME}/bin/prepDE.py -i ${gtflist}
if [ $? -ne 0 ]; then
    echo "prepDE.py failed to make count tables. Exiting."
    exit 1
fi

#Making sure output is in the right place:
mv gene_count_matrix.csv transcript_count_matrix.csv ${count_res}  
if [ $? -ne 0 ]; then
    echo "Failed to move count matrices to count_matrices directory. Exiting."
    exit 1
fi

echo "Transcriptome construction and count matrix generation complete for full set of reads!"
