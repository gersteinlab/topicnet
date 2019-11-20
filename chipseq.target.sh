#!/bin/bash


### proj directory
/ysm-gpfs/pi/gerstein/dl598/project/SKL_mm10

### list of bed files
/ysm-gpfs/pi/gerstein/from_louise/skltmp/mouse/
ls -1 /ysm-gpfs/pi/gerstein/from_louise/skltmp/mouse/*.bed > bedfile.txt

##human

more gencode.v19.tx.annotaton.gtf |grep "transcript_type \"protein_coding\"; transcript_status \"KNOWN\""|perl -lane 'BEGIN{$gname=""; $id=0;}if($_=~/gene_name "([^"]+)";/){ if($gname ne $1) {$gname=$1; $id=0;}else{$id++; } $rf=($F[6] eq"+")?"f":"r";  $tss=($F[6] eq "+")?($F[3]-1):($F[4]-1); print "$F[0]\t$tss\t".($tss+1)."\t$gname"."_utr5_$id"."_$F[0]"."_$tss"."_$rf\t0\t$F[6]";}  ' > hg19_gcv19_tss.bed

more hg19_gcv19_tss.bed|perl -lane ' print "$F[0]\t".($F[1]-2500)."\t".($F[1]+2500)."\t".join("\t",@F[3..$#F]);' > hg19_gcv19_tss2500.bed

###mouse
more mouse_mm9_biomart.txt |sed 1d| sort -k3,3 | perl -lane 'BEGIN{$gname=""; $id=0;} if($F[9] eq "protein_coding" && $F[8] eq "KNOWN") { if ($gname ne $F[2]){ $gname=$F[2]; $id=0;} else{ $id++;}  $chr="chr$F[3]"; $strand=($F[6]==1)?"+":"-"; $tss=($strand eq "+")?($F[4]-1):($F[5]-1); $flag=($strand eq "+")?"f":"r"; print "$chr\t$tss\t".($tss+1)."\t$gname"."_utr5_$id"."_$chr"."_$tss"."_$flag\t0\t$strand";} ' > mm9_biomart_tss.bed

more  mm9_biomart_tss.bed| perl -lane ' print "$F[0]\t".($F[1]-2500)."\t".($F[1]+2500)."\t".join("\t",@F[3..$#F]);' > mm9_biomart_tss2500.bed
### making annotation
awk -F'\t' 'BEGIN {OFS="\t"}{if($6=="+") print $1,$2,$2+1,$4,$5,$6; else print $1,$3-1,$3,$4,$5,$6}' ucsc_refseq_5utr.bed > ucsc_refseq_tss.bed
awk -F'\t' 'BEGIN {OFS="\t"}{if($6=="+") print $1,$2-2500,$2+2500,$4,$5,$6; else print $1,$3-2500,$3+2500,$4,$5,$6}' ucsc_refseq_5utr.bed > ucsc_refseq_tss2500.bed


##
grep chr *ENCFF*.bed | awk -F'.bed:' 'BEGIN {OFS="\t"}{print $1,$2}' | awk -F'\t' 'BEGIN {OFS="\t"}{print $2,$3,$4,$1,$6,$7,$8,$9,$10,$11}' > /ysm-gpfs/pi/gerstein/dl598/project/SKL_mm10/encode_mm10_tf_merged.bed

sort -k1,1 -k2,2n encode_mm10_tf_merged.bed > encode_mm10_tf_merged_sorted.bed

sort -k1,1 -k2,2n ucsc_refseq_tss.bed > ucsc_refseq_tss_sorted.bed

awk -F'\t' 'BEGIN {OFS="\t"}{print $4}' encode_mm10_tf_merged.bed | sort | uniq | wc -l
awk -F'\t' 'BEGIN {OFS="\t"}{print $1}' encode_mm10_tf_merged.bed | sort | uniq
awk -F'\t' 'BEGIN {OFS="\t"}{print $1}' ucsc_refseq_tss2500.bed | sort | uniq
awk -F'\t' 'BEGIN {OFS="\t"}{print $1}' ucsc_refseq_tss.bed | sort | uniq


#### memory 2fold larger than;
#encode_mm10_tf_merged.bed peakfile (shaoke provide 4column); 
bedtools window -w 2500 -a encode_mm10_tf_merged.bed -b ucsc_refseq_tss.bed > encode_mm10_tf2gene_tss2500.bed

bedtools closest -a encode_mm10_tf_merged_sorted.bed -b ucsc_refseq_tss_sorted.bed -D b > encode_mm10_tf2gene_dist2tss.bed

awk -F'\t' 'BEGIN {OFS="\t"}{print $17}' encode_mm10_tf2gene_dist2tss.bed | head

bedtools overlap -i encode_mm10_tf2gene_tss2500.bed -cols 2,3,12,13 | awk -F'\t' 'BEGIN {OFS="\t"}{print $17}' > encode_mm10_tf2gene_overlap.txt

awk -F'\t' 'BEGIN {OFS="\t"}{print $4,$14}' encode_mm10_tf2gene_tss2500.bed | awk -F'_utr5_' 'BEGIN {OFS="\t"}{print $1}' | sort | uniq > ENCODE_mm10_edgeList_TF2GENE_TSS2500.tsv







