
# download bcf file with firefox from https://dest.bio/data-files/SNP-tables#bcf-files

# create index of bcf-file
bcftools index dest.PoolSeq.PoolSNP.001.50.10Nov2020.header.bcf

# count number of variants
bcftools plugin counts dest.PoolSeq.PoolSNP.001.50.10Nov2020.header.bcf 

# Number of samples: 246
# Number of SNPs:    4042456
# Number of INDELs:  0
# Number of MNPs:    0
# Number of others:  0
# Number of sites:   4042456


# save header for easier understanding of file
bcftools view -h dest.PoolSeq.PoolSNP.001.50.10Nov2020.header.bcf > header_dest.PoolSeq.PoolSNP.001.50.10Nov2020.header.txt

# create list with sample names from Linvilla
bcftools query dest.PoolSeq.PoolSNP.001.50.10Nov2020.header.bcf -l | grep ^PA_li_ > samples_Linvilla.txt


# create new bcf-file with only samples from Linvilla
bcftools view -S samples_Linvilla.txt dest.PoolSeq.PoolSNP.001.50.10Nov2020.header.bcf -Ob -o linvilla.bcf


bcftools plugin counts linvilla.bcf

# Number of samples: 15
# Number of SNPs:    4042456
# Number of INDELs:  0
# Number of MNPs:    0
# Number of others:  0
# Number of sites:   4042456



# exclude indels (if there were still some indels)
# bcftools view -Ou -I linvilla.bcf -O b -o linvilla.bcf

# exclude SNPs from "messy" chromosomes and only biallelic SNPs with at least one sample having total depth of "."
bcftools view -Ou -e'CHROM="X" | CHROM="Y" | CHROM="4" | FORMAT/DP="."' -m2 -M2 -v snps linvilla.bcf -O b -o linvilla_filtered.bcf

bcftools plugins counts linvilla_filtered.bcf

# Number of samples: 15
# Number of SNPs:    1686856
# Number of INDELs:  0
# Number of MNPs:    0
# Number of others:  0
# Number of sites:   1686856


# get chromosome, position, alternative counts (AD), reference counts (RD) and total depths (DP) for all samples
bcftools query -f'%CHROM\t%POS\t[%AD\t]\t[%RD\t]\t[%DP\t]\n' linvilla_filtered.bcf | head -10

