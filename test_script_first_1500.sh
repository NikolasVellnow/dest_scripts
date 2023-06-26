
# create new file with only first 5000 lines for practice (24 MiB)
zcat dest.PoolSeq.PoolSNP.001.50.10Nov2020.ann.vcf.gz | head -n 5000 > first_5000.vcf

# compress the new file (5 MiB)
bgzip first_5000.vcf

# use "view" to convert compressed vcf-file to bcf-file (4.9 MiB)
bcftools view first_5000.vcf.gz -O b -o first_5000.bcf.gz       # -O b sets output file format to bcf


# create list with sample names from Linvilla
bcftools query first_5000.bcf.gz -l | grep ^PA_li_ > samples_Linvilla.txt

# create new bcf-file with only samples from Linvilla
bcftools view -Ou -S samples_Linvilla.txt first_5000.bcf.gz -O b -o first_5000_Linvilla.bcf.gz


# create index of bcf-file
bcftools index first_5000_Linvilla.bcf.gz

# save header for easier understanding of file
bcftools view -h first_5000_Linvilla.bcf.gz > header_first_5000_Linvilla.txt

# exclude indels (although I think there where already no indels anymore)
bcftools view -Ou -I first_5000_Linvilla.bcf.gz -O b -o first_5000_Linvilla_without_indels.bcf.gz

# exclude SNPs from "messy" chromosomes and only biallelic SNPs with at least one sample having total depth of "."
bcftools view -Ou -e'CHROM="X" | CHROM="Y" | CHROM="4" | FORMAT/DP="."' -m2 -M2 -v snps first_5000_Linvilla.bcf.gz -O b -o first_5000_Linvilla_filtered.bcf.gz


# get chromosome, position, alternative counts (AD), reference counts (RD) and total depths (DP) for all samples
bcftools query -f'%CHROM\t%POS\t[%AD\t]\t[%RD\t]\t[%DP\t]\n' first_5000_Linvilla_filtered.bcf.gz | head -10

