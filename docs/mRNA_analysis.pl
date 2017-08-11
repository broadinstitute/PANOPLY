#!/usr/bin/perl -w

# by Wenhu Guo (wenhu_guo@acgtinc.com)
#

foreach(glob("*.R1.fq.gz")){
	chomp;
	$_ =~ m/(.+)\.R1\.fq\.gz/;
	$read1 = $1.".R1.fq.gz";
	$read2 = $1.".R2.fq.gz";
	system "trim_galore -q 30 --length 50 --dont_gzip --paired $read1 $read2 2>$1.trim_galore.err";
	system "python /data/bin/MapSplice-v2.1.8/mapsplice.py -c /data/reference/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes -x /data/reference/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome --gene-gtf /data/reference/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf -1 $1.R1_val_1.fq -2 $1.R2_val_2.fq -p 8 --bam --non-canonical --fusion-non-canonical -o $1.mapsplice_out 2>$1.mapsplice.err";
	system "samtools sort -@ 8 $1.mapsplice_out/alignments.bam $1.mapsplice_out/alignments.sorted";
    system "samtools index $1.mapsplice_out/alignments.sorted.bam";
	system "cufflinks -o $1.cufflinks_out -p 8 -g /data/reference/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf --library-type fr-firststrand --no-update-check $1.mapsplice_out/alignments.sorted.bam 2>$1.cufflinks.err";
	system "samtools view -b -f 4 $1.mapsplice_out/alignments.sorted.bam -o $1.mapsplice_out/unmapped.bam";
    system "samtools mpileup -f /data/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa $1.mapsplice_out/alignments.sorted.bam | VarScan mpileup2cns --variants 1 --strand-filter 0 --output-vcf 1 --min-avg-qual 25 > $1.vcf";
	$cmd_t1 = q(sed 's/[";]//g;');
	$cmd_t2 = q(awk '{OFS="\t"; print $1, $4-1,$5,$12,0,$7,$14,$22,$10}');
    system "fgrep -w transcript $1.cufflinks_out/transcripts.gtf | $cmd_t1 | $cmd_t2 > $1.transcripts.bed";
    system "bedtools multicov -bams $1.mapsplice_out/alignments.sorted.bam -bed $1.transcripts.bed > $1.transcripts.ReadCount";
    system "./get_consensus_bam_batch_v7.pl -b $1.mapsplice_out/alignments.sorted.bam -t $1.transcripts.ReadCount";
	system "mv transcripts.fa $1.transcripts.fa";
	$cmd_t3 = q(grep -v "contig");
	$cmd_t4 = q(sed 's/scaffold/transcript/g');
	system "/data/bin/assemblathon_stats.pl $1.transcripts.fa | $cmd_t3 | $cmd_t4 > $1.transcripts.stat";
    system "/data/bin/transdecoder/TransDecoder -t $1.transcripts.fa --search_pfam /data/bin/transdecoder/pfam/Pfam-AB.hmm.bin --CPU 8";
}