#!/usr/bin/bash

gff_path=$1
gtf_output_path=$2

# Convert gff to gtf
gffread -F --keep-genes $gff_path > temp_file

# Remove CDS lines, nuclear_mt_pseudogene, problematic genes (gene that have two or more entries with same gene_id but different coordinates) and genes with gene_group feature
grep -vE "CDS|nuclear_mt_pseudogene|SPNCRNA.103.2|SPAC688.08.1|SPAC29E6.03c.1|SPBC16E9.16c.1|SPBC14C8.09c.1|Parent=SPBC16D10.10.1|SPNCRNA.576|SPNCRNA.584|#|gene_group" temp_file > $gtf_output_path
rm temp_file

# Change ID for gene_id in gene lines
sed -i -E 's/ID=(.*);ge/gene_id "\1";ge/g; s/;/; /g' $gtf_output_path

# Add gene_name to gene lines that don't have Name attribute
awk -i inplace '$3=="gene" && $0 !~ /Name/ {sub(/geneID=/,"gene_name ")}1' $gtf_output_path
sed -i -E 's/gene_name ([^\t]*)/gene_name "\1";/g' $gtf_output_path

# Change Name for gene_name
sed -i -E 's/Name=([^\t]*)/gene_name "\1";/g' $gtf_output_path

# Remove geneID from gene lines
sed -i -E 's/ geneID=[a-zA-Z0-9\.]*;//g' $gtf_output_path

#Invert gene lines that are problematic (are not direcly before their transcript/exon lines)
sed -zEi 's/([^\n]*gene[^\n]*SPCC18.15[^\n]*)\n([^\n]*\n[^\n]*)/\2\n\1/g' $gtf_output_path
sed -zEi 's/([^\n]*gene[^\n]*SPAC959.04c[^\n]*)\n([^\n]*\n[^\n]*)/\2\n\1/g' $gtf_output_path
sed -zEi 's/([^\n]*exon[^\n]*3391750[^\n]*SPAC959.04c[^\n]*)\n([^\n]*3393125[^\n]*SPAC959.04c[^\n]*)\n([^\n]*ID=SPAC959.04c[^\n]*)/\3\n\1\n\2/g' $gtf_output_path

# Add gene_biotype and fake transcript line to genes that don't have transcript lines nor gene_biotype
sed -zEi 's/([^\n]*PomBase.gene[^\n]*)\n([^\n]*PomBase.gene[^\n]*)/\1 gene_biotype "unknown"TEMP\n\1 gene_biotype "unknown"; transcript_biotype "unknown"\n\2/g; s/TEMP(\n[^\n]*)gene([^\n]*)gene_id ("[a-zA-Z0-9_\.\-]*)"([^\n]*transcript_biotype "unknown")/\1transcript\2transcript_id \3\.1"; gene_id \3"\4/g' $gtf_output_path

# Add gene_biotype and transcript_biotype to transcript lines and change feature for transcript
biotype=('mRNA' \
        'ncRNA' \
        'pseudogenic_transcript' \
        'rRNA' \
        'snoRNA' \
        'snRNA' \
        'tRNA')
for biot in ${biotype[@]}; do
  sed -i -E "s/($biot[^\n]*)/\1; gene_biotype \"$biot\"; transcript_biotype \"$biot\"/g" $gtf_output_path
  sed -i -E "s/\t$biot\t/\ttranscript\t/g" $gtf_output_path
done


# Add gene_biotype to gene lines
sed -zEi 's/;\n/;REPLACE\n/g; s/REPLACE\n([^\n]*(gene_biotype "[a-zA-Z0-9_\-]*"))/ \2\n\1/g' $gtf_output_path

# Add gene_name. gene_id and transcript_id to transcript lines
sed -zEi 's/(gene_name "[a-zA-Z0-9_\.\-]*"; )([^\n]*\n[^=]*)ID=([a-zA-Z0-9_\.\-]*)/\1\2transcript_id "\3"; \1/g; s/[ ]*; Parent=([[a-zA-Z0-9_\.\-]*)/ gene_id "\1"/g' $gtf_output_path

# Add gene_name to alternative transcripts within genes with more than one transcript and to all exon lines
gawk -i inplace 'BEGIN{gene_id="";gene_name=""}{if (match($0, /gene_id "([a-zA-Z0-9_\.\-]*)"/, temp) && (match($0, /gene_name "([a-zA-Z0-9_\.\-]*)"/, tempo))) {gene_id=temp[1]; gene_name=tempo[1]; print $0} else print $0"; gene_name ""\""gene_name"\""}' $gtf_output_path

# Patch gene name for the transcripts and exons of this gene (SPBC2D10.10c)
sed -Ei 's/("*SPBC2D10.10c.[^\n]*)gene_name.*/\1gene_name "fib1"/g' $gtf_output_path

# Change ID for transcript_id (for genes with multiple transcripts)
sed -Ei 's/ID=([a-zA-Z0-9\.]*) gene_id/transcript_id "\1"; gene_id/g' $gtf_output_path

# Add gene_id, gene_biotype and transcript_biotype to existing exon lines
gawk -i inplace 'BEGIN{gene_id="";gene_biotype=""}{if (match($0, /gene_id "([a-zA-Z0-9_\.\-]*)"/, temp) && (match($0, /gene_biotype "([a-zA-Z0-9_\.\-]*)"/, tempo))) {gene_id=temp[1]; gene_biotype=tempo[1]; print $0} else print $0"; gene_id ""\""gene_id"\"""; gene_biotype ""\""gene_biotype"\"""; transcript_biotype ""\""gene_biotype"\""}' $gtf_output_path

# Add exon lines to genes that miss it
sed -zEi 's/([^\n]*PomBase.transcript[^\n]*)(\n[^\n]*PomBase.gene)/\1\n\1EXONN\2/g; s/([^\n]*PomBase.)transcript([^\n]*)EXONN/\1exon\2/g' $gtf_output_path

# Change Parent for transcript_id to existing exon lines
sed -Ei 's/Parent=([a-zA-Z0-9_\.\-]*);/transcript_id "\1";/g' $gtf_output_path

# Add exon_number to exon lines
gawk -i inplace 'BEGIN{trans_id="";exon_nb=1}{if (match($0, /transcript_id "([a-zA-Z0-9_\.\-]*)"/, temp) && match($0, /exon/, tempo)) {if (temp[1]==trans_id){exon_nb+=1; print $0"; exon_number ""\""exon_nb"\""} else {trans_id=temp[1]; exon_nb=1; print $0"; exon_number ""\""exon_nb"\""}} else print $0}' $gtf_output_path

# Add exon_id to exon lines (it's the transcript id plus a .[0-9] depending on the number of exon in that transcript)
sed -Ei 's/([^\n_]*exon[^\n]*)transcript_id "([a-zA-Z0-9_\.\-]*)"([^\n]*)exon_number "([0-9]*)"/\1exon_id "\2\.\4"; transcript_id "\2"\3exon_number "\4"/g' $gtf_output_path
