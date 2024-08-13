# Bulk-ATAC-seq
```Linux
#!/bin/bash

# 1.mapping
cd './fastq'
ls *_R1.fq.gz | while read id
do 
	headname=${id%_R1.fq.gz} 
	trimmomatic PE -threads 40 ${headname}_R1.fq.gz ${headname}_R2.fq.gz ${headname}.clean_1.fastq ${headname}.clean.unpaired_1.fastq ${headname}.clean_2.fastq ${headname}.clean.unpaired_2.fastq ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36
done

if [ ! -d ../mapping ]
then
	mkdir ../mapping
fi


ls *.clean_1.fastq |while read id
do
    headname=${id%.clean_1.fastq} 
    bowtie2 --threads 40 \
	-x /data/labShare/Sequencing/GENOME/mm10/mm10 \
	-1 ${headname}.clean_1.fastq \
	-2 ${headname}.clean_2.fastq \
	-S ../mapping/${headname}.sam
	rm ${headname}.clean_1.fastq
	rm ${headname}.clean_2.fastq
	rm ${headname}.clean.unpaired_1.fastq
	rm ${headname}.clean.unpaired_2.fastq
done



# 2.SamToBam & filter low-quality mapping
cd '../mapping'
ls *.sam |while read id
do
    samtools view -@ 40 -b -F 4 -S ${id} -o ./${id%.sam}.bam
    samtools sort -@ 40 ./${id%.sam}.bam -o ./${id%.sam}.sorted.bam
	rm ${id}
done

# 3.de duplicates reads
if [ ! -d ../dedup ]
then
	mkdir ../dedup
fi


ls *.sorted.bam |while read id
do
	picard MarkDuplicates REMOVE_DUPLICATES=true I=./${id} O=../dedup/${id%.sorted.bam}.dedup.bam M=../dedup/${id%.sorted.bam}.marked_dup_metrics.txt
done

# 4.convert to BigWig
cd '../dedup'
if [ ! -d ../bg ]
then
	mkdir ../bg
fi
ls *.dedup.bam |while read id
do
    samtools index ./${id}
    bamCoverage -b ./${id} -o ../bg/${id%.dedup.bam}.bw -bs 1 --normalizeUsing BPM -p 40
done

# 5.calling peaks
if [ ! -d ../macs2 ]
then
	mkdir ../macs2
fi

ls *.dedup.bam |while read id
do
	macs2 callpeak -t ${id} -n sample --shift -75 --extsize 150 --nomodel -B --SPMR -g mm --outdir ../macs2/${id%.dedup.bam}.Macs2_out -q 0.01
	
done

# 6.Motif analysis
if [ ! -d ../peak_motif ]
then
mkdir ../peak_motif
fi

ls ../macs2/*Macs2_out/sample_peaks.narrowPeak | while read id
do
findMotifsGenome.pl ${id} mm10 ../peak_motif/ -size given -len 8,10,12 -p 40
done

# 7.TSS heatmap-deeptools
if [ ! -d ../TSS ]
then
mkdir ../TSS
fi
ls ../macs2/*Macs2_out/sample_peaks.narrowPeak | while read id
computeMatrix reference-point  --referencePoint TSS  -p 40 -b 3000 -a 3000 -R ../macs2/*.Macs2_out/sample_summits.bed -S ../bg/*.bw --skipZeros  -o *_TSS.gz --outFileSortedRegions *_TSS.bed
plotHeatmap -m *_TSS.gz -out *_TSS_Heatmap.png

echo ""
echo "WORK DONE!"
echo "Good luck,"
echo "ZhongLab Server"
```
