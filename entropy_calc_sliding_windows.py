import math
import time
from argparse import ArgumentParser
from itertools import starmap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def calculate_entropy_1loc(x_i: float) -> float:
    """ Calculate entropy for allele frequency at one locus """
    y_i = 1-x_i
    entropy_1loc = 0.0
    if (0.0 < x_i < 1.0):
        entropy_1loc = -x_i * math.log(x_i) - y_i * math.log(y_i)
    return entropy_1loc


def calculate_entropy(entropy_list: list) -> float:
    """ Calculate pseudo entropy over number of l loci """
    l = len(entropy_list)
    entropy = 0.0
    sum_value = sum(entropy_list)
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
    returns a list of tuples with window start and stop values
    """
    window_limits = range(min(genome_positions), max(
        genome_positions) + win_size, win_size//2)
    windows = list(zip(window_limits, window_limits[2:]))
    return windows


def plot_delta_h_trajectories_per_chrom(chromosome_data, window_lims: list, win_size: int, chromosome_name: str) -> None:
    """
    Loop through windows, calculate a entropy time series per window
    and then add each delta-entropy time series to plot
    """
    empty_windows_counter = 0
    # loop through windows
    for window in window_lims:
        win_data = chromosome_data[chromosome_data['pos'].between(
            window[0], window[1])]
        # if there are no SNPs in the window continue with next window
        if win_data.empty:
            empty_windows_counter += 1
            continue

        entropy_series = win_data.groupby('time_point')['entropy'].aggregate(
            calculate_entropy)

        # entropy_series = data_by_time_point['entropy'].aggregate(
        # calculate_entropy)

        entropy_t0 = entropy_series[0]
        delta_entropy_series = entropy_t0 - entropy_series

        # entropy_series.plot(alpha=0.3)
        delta_entropy_series.plot()

    print(f'{empty_windows_counter} windows had no SNPs'.format(
        empty_windows_counter))

    title_str = r'$\Delta$-Pseudo-Entropy (Linvilla)' + \
        f' for chromosome {chromosome_name} \n and fixed window size of {win_size} bases'.format(
            chromosome_name, win_size)
    plt.ylabel(r'$\Delta$-Pseudo-Entropy')
    plt.axis([0, 14, -0.7, 0.71])
    plt.yticks(np.arange(-0.7, 0.71, 0.1))
    plt.title(title_str)
    plt.savefig(title_str, dpi=200, format='png')
    # plt.show()


start_time = time.time()

PARSER = ArgumentParser()
PARSER.add_argument("-i", "--input_file", required=True)
PARSER.add_argument("-o", "--output_file", required=False)
ARGS = PARSER.parse_args()

data_file = ARGS.input_file

df = pd.read_csv(data_file)


# increase memory efficiency
df['chrom'] = df['chrom'].astype('category')
df['pos'] = df['pos'].astype('int64')
df['time_point'] = df['time_point'].astype('int8')
df['alt_counts'] = df['alt_counts'].astype('int16')
df['ref_counts'] = df['ref_counts'].astype('int16')
df['tot_counts'] = df['tot_counts'].astype('int16')

# print(df.info(memory_usage='deep'))
# print(df.describe())
# print(df['chrom'].unique())

# plt.scatter(df['pos'], df['entropy'], s=0.01, alpha=0.7)
# plt.show()

# get chromosome names
chrom_names = list(df['chrom'].unique())

# Group data by time point
grouped_data = df.groupby('chrom')

window_size = 50000

# loop over chromosomes
for chromosome in chrom_names:
    print('Chromosome: ', chromosome)
    chrom_data = grouped_data.get_group(chromosome)

    positions = chrom_data['pos'].unique()

    start_win_lim = time.time()
    # fixed window size
    win_lims = get_win_lim_fix(positions, window_size)
    end_win_lim = time.time()

    print("win_lim: ", (end_win_lim - start_win_lim))

    # fixed overlapping window size
    # overlap_win_lims = get_overlapping_win_lim_fix(positions, window_size)

    start_plot_per_chrom = time.time()
    # plot delta-entropy trajectories
    plot_delta_h_trajectories_per_chrom(
        chrom_data, win_lims, window_size, chromosome)

    end_plot_per_chrom = time.time()

    print("plot_per_chrom: ", (end_plot_per_chrom - start_plot_per_chrom))

end_time = time.time()

print('Elapsed time: ', (end_time-start_time), 's')
# df.groupby('chrom').plot(x='pos', y='entropy',
# kind='scatter', s=2, alpha=0.5, xlim=(0, max(df['pos'])), ylim=(0.0, 0.7))
# plt.vlines(x=win_lims, ymin='0.0', ymax='0.7')
# plt.vlines(x=win_lims, ymin='0.0', ymax='0.7')
# plt.vlines(x=win_lims, ymin='0.0', ymax='0.7')
# plt.show()
