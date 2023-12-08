import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import pandas as pd
import os

DIR = 'C:/Users/Ravil/source/repos/matmod_lab3/matmod_lab3/'
subdir = 'task5_b_arr_2/'
# subdir = 'task_a_arr/'
fullpath = DIR+subdir


def order(s):
    return int(s.split('_')[0])


ex = False
filelist = os.listdir(fullpath)
filelist = sorted(filelist, key=order)
for i in range(1, len(filelist) // 3 + 1):

    with open(f'{fullpath}{i}_pars.out') as fin:
        _, _a, _b = map(float, fin.readline().split())
    results_df = pd.read_csv(f'{fullpath}{i}_res.out', index_col=False)
    checkp_cnt = (results_df.shape[1]-1)//2

    max_abs = max([results_df[f'v_{i}'].abs().max()
                   for i in range(checkp_cnt)])

    # q = np.zeros((checkp_cnt, results_df[f'v_{0}'].to_numpy().shape[0]))
    # for j in range(0, checkp_cnt, 1):
    #     max_abs = max(max_abs, results_df[f'v_{j}'].abs().max())
    #     q[i] = results_df[f'v_{j}'].to_numpy()

    st = 1000

    fig, ax = plt.subplots()
    line, = ax.plot(results_df[f'v_{st + 0}'].to_numpy())
    ax.set_title(f"alpha = {_a}, beta = {_b}")
    ax.set_xlim([-1, 501])
    ax.set_ylim([-max_abs, max_abs])

    def update(frame):
        # if frame % 10 == 0:
        #     print(frame)
        line.set_ydata(results_df[f'v_{st + frame}'].to_numpy())
        # line.set_ydata(q[frame, :])
        return line,
    anim = animation.FuncAnimation(fig=fig, func=update, frames=checkp_cnt-st,
                                   interval=30, blit=True, repeat=False)
    plt.show()
    plt.close()
