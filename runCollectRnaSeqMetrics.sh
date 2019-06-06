#----------------------------------------------
# Run CollectRnaSeqMetrics tool on bam files
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_analysis_CollectRnaSeqMetrics.php
#----------------------------------------------

# Convert GTF reference to GenePred
./gtfToGenePred Homo_sapiens.GRCh38.94.gtf Homo_sapiens.GRCh38.94.GenePred

# Convert GenePred to REF_FLAT (http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat)
	# converts GenePred format to REF_FLAT - the only difference is that GenePred is missing the second column of REF_FLAt ("Name of gene as it appears in Genome Browser.")
  
	#This duplicates the first column ("Name of gene") of GenePred, so it has the same number of columns as REF_FLAT: 
	awk 'BEGIN { FS="\t"; OFS="\t" } { $1=$1 "\t" $1 } 1' Homo_sapiens.GRCh38.94.GenePred > Homo_sapiens.GRCh38.94.refFlat
  
# Create IntervalList (adapted from Kamil: https://gist.github.com/slowkow/b11c28796508f03cdf4b)

	chrom_sizes=chrom.sizes.GRCh38 #downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
	genes=Homo_sapiens.GRCh38.94.gtf

	# Copy header from one of the bams (interval_list must have the same header as the bams!)
	samtools view -H bam-name.bam > gencode.GRCh38.rRNA.interval_list

	# Intervals for rRNA and rRNA_pseudogene. If something else is needed, you need to change what is grepped.
	grep 'gene_biotype "rRNA"\|gene_biotype "rRNA_pseudogene"' $genes | \
	    awk '$3 == "gene"' | \
	    cut -f1,4,5,7,9 | \
	    perl -lane '/gene_name "([^"]+)"/ or die "no gene_name on $.";
	        print join "\t", (@F[0,1,2,3], $1)' | \
	    sort -k1V -k2n -k3n \
	>> gencode.GRCh38.rRNA.interval_list
  
# Run CollectRnaSeqMetrics on all bams
	for file in *.bam 
	do	
	java -jar ./picard.jar CollectRnaSeqMetrics \
	      I=$file \
	      O=$file.RNA_Metrics_Ensembl \
	      REF_FLAT=Homo_sapiens.GRCh38.94.refFlat \
	      STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
	      RIBOSOMAL_INTERVALS=gencode.GRCh38.rRNA.interval_list
	done
 
