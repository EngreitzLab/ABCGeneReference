# download data from GENCODE
wget https://www.encodeproject.org/files/ENCFF159KBI/@@download/ENCFF159KBI.gtf.gz
mv ENCFF159KBI.gtf.gz > gencode.v29.gtf

# protein coding transcripts, columns: chr, start, end, gene, score, strand, gene type
awk '{if($3=="transcript" && $0~"level (1|2);" && $0~"basic" && $0~"protein_coding"){print $0}}' gencode.v29.gtf | awk 'OFS="\t" {print $1,$4-1,$5,$16,$22,$7,$14}' | tr -d '";' | sort -k1,1 -k2,2n > gencode.v29.transcripts.level12.basic.protein_coding.bed

# add TSS
Rscript gene_to_tss_list.R --infile gencode.v29.transcripts.level12.basic.protein_coding.bed --outfile gencode.v29.transcripts.level12.basic.protein_coding.TSS500bp.bed

# make altTSS selection file, columns: chr, TSS start, TSS end, TSS, score, strand, gene, gene start, gene end
Rscript transcripts_for_altTSS_selection.R --infile gencode.v29.transcripts.level12.basic.protein_coding.bed --outfile gencode.v29.transcripts.levels12.basic.protein_coding.TSS500bp_with_gene.forAltTSSSelection.txt

# make version with lincRNAs
awk '{if($3=="transcript" && $0~"level (1|2);" && $0~"basic" && ($0~"protein_coding" || $0~"lincRNA")){print $0}}' gencode.v29.gtf | awk 'OFS="\t" {print $1,$4-1,$5,$16,$22,$7,$14}' | tr -d '";' | sort -k1,1 -k2,2n > gencode.v29.transcripts.level12.basic.protein_coding_lincRNA.bed

# add TSS
Rscript gene_to_tss_list.R --infile gencode.v29.transcripts.level12.basic.protein_coding_lincRNA.bed --outfile gencode.v29.transcripts.level12.basic.protein_coding_lincRNA.TSS500bp.bed

# use Rscript choose_TSS.R to select single promoter/transcript per gene
Rscript choose_TSS.R --infile gencode.v29.transcripts.level12.basic.protein_coding.bed --outfile gencode.v29.transcripts.level12.basic.protein_coding.unique_TSS.bed --geneOutFile gencode.v29.transcripts.level12.basic.protein_coding.unique_TSS.gene_bodies.bed