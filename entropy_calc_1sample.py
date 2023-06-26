from argparse import ArgumentParser
from decimal import * # needed to avoid weird rounding errros/bias
import subprocess
import sys

getcontext().prec = 5 # setting decimal precision to 5 decimals

def create_snp_dict(line: str, snp_keys: list) -> dict:
    snp_data = line.split('\t')
    snp_dict = {k:v for (k,v) in zip(snp_keys,snp_data)}
    return snp_dict

def create_snps_list(output: str, snp_keys: list):
    snps_list = [create_snp_dict(line, snp_keys) for line in output.split('\n') if len(line)>1]
    return snps_list

def calculate_entropy_1loc(snp: dict) -> Decimal:
    x_i = Decimal(snp.get('alt_counts'))/Decimal(snp.get('tot_counts'))
    y_i = Decimal(snp.get('ref_counts'))/Decimal(snp.get('tot_counts'))
    entropy_1loc = 0
    if(x_i > Decimal(0.0) and x_i < Decimal(1.0)):
        entropy_1loc = -x_i * x_i.ln() - y_i * y_i.ln()
    return entropy_1loc


def calculate_entropy(snps_list: list) -> Decimal:
    l = len(snps_list)
    entropy = Decimal(0.0)
    sum_value = sum(map(calculate_entropy_1loc, snps_list))
    entropy = sum_value/l
    return entropy


parser = ArgumentParser()
parser.add_argument("-i", "--input_file", required=True)
#parser.add_argument("-o", "--output_file", required=True)
args = parser.parse_args()


bcf_file_str = args.input_file


# build bcftools command to get SNP data

#
command = "bcftools query -f'%CHROM\\t%POS\\t[%AD]\\t[%RD]\\t[%DP]\\n' " + bcf_file_str
#bcftools query -f'%CHROM\t%POS\t[%AD]\t[%RD]\t[%DP]\n' first_5000_PA_li_09_fall.bcf.gz | head -10

# save output stream from slim run to work with
output_stream = subprocess.run(command, capture_output=True, shell=True, text=True)
output = output_stream.stdout

snp_keys = ['chrom', 'pos', 'alt_counts', 'ref_counts', 'tot_counts']

snps_list = create_snps_list(output, snp_keys)

print(sys.getsizeof(snps_list))

entropy_0 = calculate_entropy(snps_list=snps_list)
print(entropy_0)



