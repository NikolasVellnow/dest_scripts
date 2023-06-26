# download bcf file with firefox from https://dest.bio/data-files/SNP-tables#bcf-files


# create index of bcf-file
bcftools index -f dest.PoolSeq.PoolSNP.001.50.10Nov2020.header.bcf

# make subset
bcftools view -Ou -r 2L:1-20000 dest.PoolSeq.PoolSNP.001.50.10Nov2020.header.bcf -Ob -o subset.bcf.gz


# create index of subset bcf-file
bcftools index -f subset.bcf.gz


# count number of variants
bcftools plugin counts subset.bcf.gz

# Number of samples: 246
# Number of SNPs:    207
# Number of INDELs:  0
# Number of MNPs:    0
# Number of others:  0
# Number of sites:   207


# create list with sample names from Linvilla
bcftools query dest.PoolSeq.PoolSNP.001.50.10Nov2020.header.bcf -l | grep ^PA_li_ > samples_Linvilla.txt


# create new bcf-file with only samples from Linvilla (gives error when output-type is b=compressed => try u)
bcftools view -S samples_Linvilla.txt subset.bcf.gz -Ob -o linvilla_subset.bcf.gz

# create index of subset bcf-file
bcftools index -f linvilla_subset.bcf.gz

# save header for easier understanding of file
bcftools view -h linvilla_subset.bcf.gz > header_linvilla_subset.txt

# count number of variants
bcftools plugin counts linvilla_subset.bcf.gz

# exclude indels (although I think there where already no indels anymore)
bcftools view -Ob -I linvilla_subset.bcf.gz -Ob -o linvilla_subset.bcf.gz

# exclude SNPs from "messy" chromosomes and only biallelic SNPs with at least one sample having total depth of "."
bcftools view -Ob -e'CHROM="X" | CHROM="Y" | CHROM="4" | FORMAT/DP="."' -m2 -M2 -v snps linvilla_subset.bcf.gz -Ob -o linvilla_subset_filtered.bcf.gz

# count number of variants again
bcftools plugin counts linvilla_subset_filtered.bcf.gz

# get chromosome, position, alternative counts (AD), reference counts (RD) and total depths (DP) for all samples
bcftools query -f'%CHROM\t%POS\t[%AD\t]\t[%RD\t]\t[%DP\t]\n' linvilla_subset_filtered.bcf.gz | head -10

