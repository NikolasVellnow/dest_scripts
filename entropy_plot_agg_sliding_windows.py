"""
Plots line plots and violine plots of entropy trajectories
"""
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PARSER = ArgumentParser()
PARSER.add_argument("-i", "--input_file", required=True)
ARGS = PARSER.parse_args()

data_file = ARGS.input_file

df = pd.read_csv(data_file)

title_str = data_file[:-4]


windows = df.groupby('window')
for window in windows:
    plt.plot(window[1].time_point, window[1].delta_entropy)

# dummy = [plt.plot(window[1].time_point, window[1].delta_entropy)
    # for window in windows]

plt.figure(1)
plt.ylabel(r'$\Delta$-Pseudo-Entropy')
plt.xlabel('Time points')
plt.axis([-0.6, 14.6, -0.7, 0.5])
plt.yticks(np.arange(-0.7, 0.5, 0.1))
plt.title(title_str)
plt.grid(True, color='grey', linestyle='dotted', axis='y')
plt.savefig(title_str + '_line_plot.png', dpi=300, format='png')

series = df.groupby('time_point')

trajectories = [elem[1].delta_entropy for elem in series]

plt.figure(2)
bla = plt.violinplot(trajectories,
                     df.time_point.unique(),
                     showmedians=True,
                     showextrema=False,
                     widths=0.8)


for pc in bla['bodies']:
    pc.set_facecolor('0.4')
    pc.set_edgecolor('0.2')
    pc.set_alpha(0.7)

plt.ylabel(r'$\Delta$-Pseudo-Entropy')
plt.xlabel('Time points')
plt.axis([-0.6, 14.6, -0.7, 0.5])
plt.yticks(np.arange(-0.7, 0.5, 0.1))
plt.title(title_str)
plt.grid(True, color='grey', linestyle='dotted', axis='y')
plt.savefig(title_str + '_violine_plot.png', dpi=300, format='png')
