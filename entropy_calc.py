from argparse import ArgumentParser
import subprocess


parser = ArgumentParser()
parser.add_argument("-i", "--input_file", required=True)
#parser.add_argument("-o", "--output_file", required=True)
args = parser.parse_args()


bcf_file_str = args.input_file


# build bcftools command to get SNP data

#command = "bcftools query -f'%CHROM\\t%POS\\t[%AD\\t]\\t[%RD\\t]\\t[%DP\\t]\\n' " + bcf_file_str + " | head -10"


#
command = "bcftools query -f'%CHROM\\t%POS\\t[%AD]\\t[%RD]\\t[%DP]\\n' " + bcf_file_str + " | head -10"
#bcftools query -f'%CHROM\t%POS\t[%AD]\t[%RD]\t[%DP]\n' first_5000_PA_li_09_fall.bcf.gz | head -10
print(command)

# save output stream from slim run to work with
output_stream = subprocess.run(command, capture_output=True, shell=True, text=True)
output = output_stream.stdout


for line in output:
    print(line)
