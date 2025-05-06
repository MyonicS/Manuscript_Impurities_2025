from matplotlib import pyplot as plt


def set_plot_params():
    font = {'family' : 'Arial',
            'weight' : 'normal',
            'size'   : 7}
    plt.rc('font', **font)
    plt.rcParams['figure.figsize'] = [8.8/2.54,6.22/2.54]
