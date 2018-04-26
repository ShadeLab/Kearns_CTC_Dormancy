# Analysis of 16S Illumina sequencing data from bean rhizosphere 

All analyses were performed in USEARCH v. 10.0.240 (64-bit) with the UPARSE OTU picking method. Subsequent analyses were performed in QIIME v. 1.8. 

## Merge Paired End Reads
```
#decompress the reads
gunzip *.gz

mkdir mergedfastq

./usearch64 -fastq_mergepairs *R1*.fastq -relabel @ -fastq_maxdiffs 10 -fastqout mergedfastq/merged.fq -fastq_merge_maxee 1.0 -fastq_minmergelen 250 -fastq_maxmergelen 300
```

## Dereplicate sequences
```
./usearch64 -fastx_uniques mergedfastq/merged.fq -fastqout mergedfastq/uniques_combined_merged.fastq -sizeout

#move raw reads to new folder
mkdir raw_reads
mv *.fastq raw_reads
```

## Remove Singeltons
```
./usearch64 -sortbysize mergedfastq/uniques_combined_merged.fastq -fastqout mergedfastq/nosigs_uniques_combined_merged.fastq -minsize 2
```

## Precluster Sequences
```
./usearch64 -cluster_fast mergedfastq/nosigs_uniques_combined_merged.fastq -centroids_fastq mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size
```

## Reference-based OTU picking against Silva v. 123
```
./usearch64 -usearch_global mergedfastq/denoised_nosigs_uniques_combined_merged.fastq -id 0.97 -db silva123.fasta  -strand plus -uc mergedfastq/ref_seqs.uc -dbmatched mergedfastq/closed_reference.fasta -notmatchedfq mergedfastq/failed_closed.fq
```

## Sort by size and then de novo OTU picking on sequences that failed to hit GreenGenes
```
./usearch64 -sortbysize mergedfastq/failed_closed.fq -fastaout mergedfastq/sorted_failed_closed.fq

#de novo pick and chimera check
./usearch64 -cluster_otus mergedfastq/sorted_failed_closed.fq -minsize 2 -otus mergedfastq/denovo_otus.fasta -relabel OTU_dn_ -uparseout denovo_out.up
```

## Combine the rep sets between de novo and reference-based OTU picking
```
cat mergedfastq/closed_reference.fasta mergedfastq/denovo_otus.fasta > mergedfastq/full_rep_set.fna
```

## Map rep_set back to pre-dereplicated sequences and make OTU tables
```
./usearch64 -usearch_global mergedfastq/merged.fq -db mergedfastq/full_rep_set.fna  -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.biom
```

# Switch to QIIME v. 1.8

## Assign taxonomy to Silva v. 123 and UCLUST
```
assign_taxonomy.py -i mergedfastq/full_rep_set.fna -o taxonomy -r silva123.fasta -t silva123_taxonomy.txt
```

## Add taxonomy to OTU table
```
biom convert -i OTU_table.txt --table-type='OTU table' -o otu_jsn.biom
biom add-metadata -i otu_jsn.biom -o otu_table_tax.biom --observation-metadata-fp=taxonomy/full_rep_set_tax_assignments.txt --sc-separated=taxonomy --observation-header=OTUID,taxonomy
```

## Filter non-bacteria/archaea and Unassigned reads form OTU table
```
filter_taxa_from_otu_table.py -i otu_table_tax.biom -o otu_table_tax_filt.biom -n o__Streptophyta, c__Chloroplast,f__mitochondria,Unassigned
```

## Align sequences to Silva with PyNast
```
align_seqs.py -i mergedfastq/full_rep_set.fna -o alignment -t silva123.fasta
```

## Filter excess gaps from alignment
```
filter_alignment.py -i alignment/full_rep_set_aligned.fasta -o alignment/filtered_alignment
```

## Make phylogeny with fasttree
```
make_phylogeny.py -i alignment/filtered_alignment/full_rep_set_aligned_pfiltered.fasta -o rep_set.tre
```

## Summarize the OTU table and rarefy OTU table to lowest sequencing depth
```
biom summarize-table -i otu_table_tax_filt.biom -o otu_table_summary.txt

single_rarefaction.py -d 13123 -o single_rare.biom -i otu_table_tax_filt.biom
```

## Calculate global alpha and beta diversity
```
beta_diversity.py -m bray_curtis,unweighted_unifrac,weighted_unifrac -i single_rare.biom -o beta_div -t rep_set.tre
principal_coordinates.py -i beta_div -o coords
```

