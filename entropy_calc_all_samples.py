"""
Script takes bcf-file as input, calculates allel frequencies and allelic entropies,
and filters out SNP positions that are not polymorphic at the first sampling time point.
The script then outputs a csv-file witht the processed data.

"""
import math
import subprocess
import time
from argparse import ArgumentParser
from itertools import chain

import pandas as pd

NUM_COL: int = 5


def calculate_entropy_1loc(x_i: float) -> float:
    """ Calculate entropy for allele frequency at one locus """
    y_i = 1-x_i
    entropy_1loc = 0.0
    if 0.0 < x_i < 1.0:
        entropy_1loc = -x_i * math.log(x_i) - y_i * math.log(y_i)
    return entropy_1loc


def calculate_entropy(entropy_list: list) -> float:
    """ Calculate pseudo entropy over number of num_l loci """
    num_l = len(entropy_list)
    entropy = 0.0
    sum_value = sum(entropy_list)
    entropy = sum_value/num_l
    return entropy


def create_snp_dict(count_tuple: tuple, snp_location: list,
                    sample_index: int, sample_list: list) -> dict:
    """ Create dictionary with data for individual snp """
    snp_dict = {}
    snp_dict.setdefault('sample', str(sample_list[sample_index]))
    snp_dict.setdefault('chrom', str(snp_location[0]))
    snp_dict.setdefault('pos', str(snp_location[1]))
    snp_dict.setdefault('alt_counts', int(count_tuple[0]))
    snp_dict.setdefault('ref_counts', int(count_tuple[1]))
    snp_dict.setdefault('tot_counts', int(count_tuple[2]))
    return snp_dict


def create_snp_dict_list(current_line: str, sample_list: list) -> list:
    """ Create list of snp-dictionaries for all samples """
    snp_dict_list = []
    num_samples = len(sample_list)
    splitted_line = current_line.split('\t')
    # remove weird empty space
    splitted_line = [x for x in splitted_line if x != '']
    snp_location = splitted_line[:2]        # chromosome and position of snp
    remaining_line = splitted_line[2:]
    count_tuples = [(alt, ref, tot) for (alt, ref, tot) in zip(
        remaining_line[0:num_samples],
        remaining_line[num_samples:num_samples*2],
        remaining_line[num_samples*2:num_samples*3])]

    snp_dict_list = [create_snp_dict(entry, snp_location, sample_index, sample_list)
                     for sample_index, entry in enumerate(count_tuples)]

    return snp_dict_list


def convert_to_dataframe(all_snp_dicts: list) -> pd.DataFrame:
    """ Convert list of list of dictionaries to pandas dataframe """
    total_snp_dict_list = list(
        chain.from_iterable(all_snp_dicts))
    dataframe = pd.DataFrame.from_records(total_snp_dict_list)
    return dataframe


def convert_samples_to_timepoints(sample: str) -> int:
    """ Convert sample names into time point intergers"""
    switch_dict = {
        'PA_li_09_spring': 0,
        'PA_li_09_fall': 1,
        'PA_li_10_spring': 2,
        'PA_li_10_fall': 3,
        'PA_li_11_spring': 4,
        'PA_li_11_fall': 5,
        'PA_li_11_frost': 6,
        'PA_li_12_spring': 7,
        'PA_li_12_sum': 8,
        'PA_li_12_fall': 9,
        'PA_li_13_spring': 10,
        'PA_li_14_spring': 11,
        'PA_li_14_fall': 12,
        'PA_li_15_spring': 13,
        'PA_li_15_fall': 14
    }
    return switch_dict.get(sample, -1)


PARSER = ArgumentParser()
PARSER.add_argument("-i", "--input_file", required=True)
PARSER.add_argument("-o", "--output_file", required=False)
ARGS = PARSER.parse_args()

BCF_FILE_STR = ARGS.input_file


SAMPLES_COMMAND = "bcftools query -l " + BCF_FILE_STR
# save output stream from bcftools to work with
SAMPLES_OUTPUT_STREAM = subprocess.run(
    SAMPLES_COMMAND, capture_output=True, shell=True, text=True, check=True)
SAMPLES_OUTPUT = SAMPLES_OUTPUT_STREAM.stdout

# save sample names from output
sample_names = [line for line in SAMPLES_OUTPUT.split('\n') if line != '']


# build bcftools command to get SNP data
QUERY = "-f'%CHROM\\t%POS\\t[%AD\\t]\\t[%RD\\t]\t[%DP\\t]\\n' "
SNP_DATA_COMMAND = "bcftools query " + QUERY + BCF_FILE_STR


# save output stream from bcftools to work with
OUTPUT_STREAM = subprocess.run(
    SNP_DATA_COMMAND, capture_output=True, shell=True, text=True, check=True)
OUTPUT = OUTPUT_STREAM.stdout

DATA_COLS = ['chrom', 'pos', 'sample',
             'alt_counts', 'ref_counts', 'tot_counts']


total_snp_dict_list_of_lists = [create_snp_dict_list(line, sample_names)
                                for line in OUTPUT.split('\n') if line != '']


data = convert_to_dataframe(total_snp_dict_list_of_lists)


start_time = time.time()

# add time point column
data['time_point'] = data['sample'].map(convert_samples_to_timepoints)

# add column to uniquely identify each position on each chromosome
data['chrom_pos'] = data['chrom'].map(str) + '_' + data['pos'].map(str)
data['chrom_pos'] = data['chrom_pos'].map(hash)

# Change datatype to save memory
data['sample'] = data['sample'].astype('category')
data['chrom'] = data['chrom'].astype('category')
data['pos'] = data['pos'].astype('int64')

data['time_point'] = data['time_point'].astype('int8')
data['alt_counts'] = data['alt_counts'].astype('int16')
data['ref_counts'] = data['ref_counts'].astype('int16')
data['tot_counts'] = data['tot_counts'].astype('int16')

# check if hashing unique position IDs worked
num_IDs = len(data['chrom_pos'])
num_unique_IDs = len(set(data['chrom_pos']))

if num_unique_IDs != len(set(data['chrom_pos'])):
    print("!!!! Warning: Hashing unique position IDs did not work!!!!")


# Group data by time point
grouped_data = data.groupby('time_point')

# get positions of snps where allele frequency is > 0.0 at first time point
first_sample = grouped_data.get_group(0)
filtered_positions = first_sample.loc[first_sample['alt_counts'] > 0, [
    'chrom_pos']].sort_values('chrom_pos')

filtered_positions_list = list(filtered_positions['chrom_pos'])

filtered_data = data[data['chrom_pos'].isin(filtered_positions_list)].copy()


# Delete original data to save memory (is this the right way?)
del data
del grouped_data
del first_sample
del filtered_positions
del filtered_positions_list

# calculate allele frequencies of alternative allele
filtered_data['freq'] = filtered_data.loc[:, 'alt_counts'] / \
    filtered_data.loc[:, 'tot_counts']
# calculate entropies for individual loci
filtered_data['entropy'] = filtered_data['freq'].map(calculate_entropy_1loc)


print("Filtered data: ")
print(filtered_data.head())
print(filtered_data.describe())

filtered_data.to_csv('linvilla_filtered_data.csv', index=False)

# subset for coding
subset_data = filtered_data.loc[filtered_data['pos'] < 40000]
subset_data.to_csv('linvilla_subset_data.csv', index=False)

end_time = time.time()

print("elapsed time:", end_time-start_time)
