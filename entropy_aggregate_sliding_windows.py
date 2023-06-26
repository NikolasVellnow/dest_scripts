import math
import time
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import pandas as pd


def calculate_entropy(entropy_list) -> float:
    """ Calculate pseudo entropy over number of l loci """
    l = len(entropy_list)
    entropy = 0.0
    sum_value = math.fsum(entropy_list)
    entropy = sum_value/l
    return entropy


def get_win_lim_fix(genome_positions, win_size: int) -> list:
    """
    Takes list or dataframe column of genome positions and fixed window size and
    returns a list of tuples with window start and stop values
    """
    window_limits = range(min(genome_positions), max(
        genome_positions) + win_size, win_size)
    windows = list(zip(window_limits, window_limits[1:]))
    return windows


def get_overlapping_win_lim_fix(genome_positions, win_size: int) -> list:
    """
    Takes list or dataframe column of genome positions and window size and
    returns a list of tuples with window start and stop values. Windows have overlap of 50%
    """
    window_limits = range(min(genome_positions), max(
        genome_positions) + win_size, win_size//2)
    windows = list(zip(window_limits, window_limits[2:]))
    return windows


start_time = time.time()

PARSER = ArgumentParser()
PARSER.add_argument("-i", "--input_file", required=True)
PARSER.add_argument("-ws", "--window_size",
                    required=False, default=5000, type=int)
PARSER.add_argument("-th", "--threshold",
                    required=False, default=0, type=int)
PARSER.add_argument("-o", "--output_file", required=False)
ARGS = PARSER.parse_args()

data_file = ARGS.input_file
window_size = ARGS.window_size
threshold = ARGS.threshold

df = pd.read_csv(data_file)

# increase memory efficiency
df['sample'] = df['sample'].astype('category')
df['chrom'] = df['chrom'].astype('category')
df['pos'] = df['pos'].astype('int64')
df['time_point'] = df['time_point'].astype('int8')
df['alt_counts'] = df['alt_counts'].astype('int16')
df['ref_counts'] = df['ref_counts'].astype('int16')
df['tot_counts'] = df['tot_counts'].astype('int16')


# get chromosome names
chrom_names = list(df['chrom'].unique())

# Group data by time point
grouped_data = df.groupby('chrom')

# loop over chromosomes
for chromosome in chrom_names:
    chrom_data = grouped_data.get_group(chromosome).copy()
    chrom_data['window'] = -1
    chrom_data['window'] = chrom_data['window'].astype('int32')
    chrom_data['window_start'] = -1
    chrom_data['window_start'] = chrom_data['window_start'].astype('int64')
    chrom_data['window_stop'] = -1
    chrom_data['window_stop'] = chrom_data['window_stop'].astype('int64')

    # get positions of SNPs
    positions = chrom_data['pos'].unique()

    # fixed window size
    win_lims = get_win_lim_fix(positions, window_size)

    # loop through windows
    for index, window in enumerate(win_lims):
        print('Chromosome: ', chromosome, window)
        # create mask for only SNPs within window
        win_data = chrom_data[chrom_data['pos'].between(window[0], window[1])]
        win_mask = chrom_data['pos'].between(window[0], window[1])
        # add columns to identify windows
        chrom_data.loc[win_mask, 'window_start'] = window[0]
        chrom_data.loc[win_mask, 'window_stop'] = window[1]
        chrom_data.loc[win_mask, 'window'] = index

    # aggregate data per window and time point
    chrom_data_per_time_point = chrom_data.groupby(
        ['window', 'time_point']).agg(
            entropy=('entropy', calculate_entropy),
            num_snps=('pos', 'count'),
            window_start=('window_start', 'first'),
            window_stop=('window_stop', 'last')
    ).reset_index()

    # increase memory efficiency
    chrom_data_per_time_point['window'] = chrom_data_per_time_point['window'].astype(
        'int32')
    chrom_data_per_time_point['time_point'] = chrom_data_per_time_point['time_point'].astype(
        'int8')
    chrom_data_per_time_point['num_snps'] = chrom_data_per_time_point['num_snps'].astype(
        'int16')
    # add column with entropy at time point 0 for each window
    chrom_data_per_time_point['t0_entropy'] = chrom_data_per_time_point.groupby('window')['entropy'].transform(
        'first')
    chrom_data_per_time_point['window_start'] = chrom_data_per_time_point['window_start'].astype(
        'int32')
    chrom_data_per_time_point['window_stop'] = chrom_data_per_time_point['window_stop'].astype(
        'int32')

    # calculate delta-entropy for each time point
    chrom_data_per_time_point['delta_entropy'] = chrom_data_per_time_point['entropy'] - \
        chrom_data_per_time_point['t0_entropy']

    # exclude windows that have less SNPs htan threshold
    if threshold != 0:
        chrom_data_per_time_point = chrom_data_per_time_point[chrom_data_per_time_point['num_snps'] > threshold]

        file_str = 'Entropy_Linvilla_agg' + \
        f'_chromosome{chromosome}_fixed_win_{window_size}_threshold{threshold}'.format(
            chromosome, window_size, threshold)
        chrom_data_per_time_point.to_csv(file_str + '.csv', index=False)
    
    else:
        file_str = 'Entropy_Linvilla_agg' + \
            f'_chromosome{chromosome}_fixed_win_{window_size}'.format(
                chromosome, window_size)
        chrom_data_per_time_point.to_csv(file_str + '.csv', index=False)


end_time = time.time()

print('Elapsed time: ', (end_time-start_time), 's')
