from matplotlib import pyplot as plt
import matplotlib.pyplot as plt
import platform

def set_plot_params():
    system = platform.system()
    
    if system == 'Windows':
        font_family = 'Arial'
    elif system == 'Linux':
        font_family = 'Liberation Sans'
    else:
        font_family = 'sans-serif'  # Fallback
    
    font = {
        'family': font_family,
        'weight': 'normal',
        'size': 7
    }
    
    plt.rc('font', **font)
    plt.rcParams['figure.figsize'] = [8.8 / 2.54, 6.22 / 2.54]