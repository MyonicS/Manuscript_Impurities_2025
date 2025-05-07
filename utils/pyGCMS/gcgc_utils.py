
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import scipy.integrate as integrate
import matplotlib.colors as colors





def parse_chromatogram(path): #from a CSV file
    df = pd.read_csv(path, sep = ',',  skiprows = 1)
    df.columns = ['Time(ms)', 'Time(min)', 'unknown','Absolute Intensity']
    df['Ret.Time'] = df['Time(min)']
    df['Ret.Time']-=df['Ret.Time'][0]
    df['Ret.Time[s]']=df['Ret.Time']*60
    df.drop(['Time(ms)','Time(min)','unknown'], axis=1, inplace=True)
    return df


def split_solvent(df, solvent_time):
    #set absoluyte intensity to 0 when time is < solvent_time
    df.loc[df['Ret.Time[s]'] <= solvent_time, 'Absolute Intensity'] = 0
    return df

def add_split(df,split_time,sampling_interval):#split time in s, sampling interval in ms
    rows_splitting = split_time/sampling_interval*1000
    df['split_no_fromindex'] = df.index//rows_splitting 
    while len(df[df['split_no_fromindex']==df['split_no_fromindex'].max()]) < len(df[df['split_no_fromindex']==df['split_no_fromindex'].min()]):
        df.loc[len(df)+1] = df.iloc[-1]
    return df

def min_correct(df): # 'global minimum' as baseline
    df['Absolute Intensity'] = df['Absolute Intensity'] - df['Absolute Intensity'].min()
    return df

def baseline_stridewise(df_array): #minimum in each stride as baseline
    df_array = df_array - df_array.min(axis=0)
    return df_array

def convert_to2D(df,split_time):
    df_short = df[['split_no_fromindex','Absolute Intensity']]
    array_list = []
    
    for i in range(0,int(df_short['split_no_fromindex'].max())):
        array_list.append(df_short[df_short['split_no_fromindex']==i]['Absolute Intensity'].values)

    #turn arraylist into an 2D array
    array = np.zeros((len(array_list),len(array_list[0])))
    for i in range(0,len(array_list)):
        array[i,:] = array_list[i]

    index_list_retention_time = []
    for i in range(len(array_list)):
        index_list_retention_time.append(i*split_time)

    columns_splittime = []
    for i in range(len(array_list[0])):
        columns_splittime.append(round(i*split_time/len(array_list[0]),3))

    df_array = pd.DataFrame(array, index = index_list_retention_time, columns = columns_splittime)
    df_array= df_array.T
    df_array = df_array.iloc[::-1]
    return df_array

def shift_phase(df_array, shift):
    df_array_shifted = np.roll(df_array, shift, axis=0)
    return df_array_shifted

import tifffile





def normalize_array(df_array):
    arrray_for_integral = np.array(df_array)
    #integrate over rows
    row_integrated_pre = integrate.trapezoid(arrray_for_integral, axis=0)
    #integrating the new array
    integral_non_norm = integrate.trapezoid(row_integrated_pre, axis=0)
    df_array_norm = df_array/integral_non_norm
    return df_array_norm
    
def integrate_masked(df_norm_array, maskpath):
    mask =  tifffile.imread(maskpath)/255 # if the mask is binary no need to divide by 255
    df_norm_array_masked  = df_norm_array*mask

    array_norm_mask_diarom = np.array(df_norm_array_masked)
    #integrate over rows
    row_integrated = integrate.trapezoid(array_norm_mask_diarom, axis=0)
    #integrating the new array
    column_integrated = integrate.trapezoid(row_integrated, axis=0)
    return column_integrated
    


def mask_integrate(df_array_norm, mask_dir):
    mask_list = glob.glob(mask_dir + '*.tif')
    #get the mask names
    mask_names = []
    for i in range(len(mask_list)):
        mask_name = mask_list[i].split('\\')[-1].split('.')[0]
        #split off the 'Mask_' part
        mask_names.append(mask_name.split('Mask_')[-1])


    integral_list = []
    for i in range(len(mask_list)):
        integral_list.append(integrate_masked(df_array_norm, mask_list[i]))

    # make a dataframe with the integral values and the mask_names as column names
    df_integral = pd.DataFrame([integral_list], columns = mask_names)
    df_integral['unassigned'] = 1-df_integral[mask_names].sum(axis=1)
    return df_integral

def process_chromatogram(filepath, split_time, sampling_interval, mask_dir,shift=0, solvent_time=0):
    df = parse_chromatogram(filepath)
    df = split_solvent(df, solvent_time)
    df = add_split(df,split_time,sampling_interval)
    df_array = convert_to2D(df,split_time)
    df_array  = baseline_stridewise(df_array)
    df_array = shift_phase(df_array, shift)
    df_array_norm = normalize_array(df_array)
    df_integral = mask_integrate(df_array_norm, mask_dir)
    return df_integral, df_array_norm

def plot_2Dchromatogram(df_array_norm,maskdir,savedir,plotmask=True,title = '2D Chromatogram',split_time = 20):
    mask_list = glob.glob(maskdir + '*.tif')
    plt.figure(figsize=(8,8/1.615))
    plt.imshow(np.sqrt(df_array_norm), cmap='viridis', interpolation='nearest', extent=[0, 106, 0, split_time], aspect='auto')
    plt.colorbar(label='$\sqrt{\mathrm{intensity}}$')

    colormaplist = [colors.ListedColormap(['none', 'C'+str(i)]) for i in range(6) ]
    annotations = [['Alkanes/Alkenes',28,7],['Monoaromatics',45,12],['Diaromatics',30,19],['Triaromatics',64,4.5],['Pyrenes',93,8]]
    if plotmask:
        for i in range(len(mask_list)):
            mask =  tifffile.imread(mask_list[i])/255 # if the mask is binary no need to divide by 255
            plt.imshow(mask, cmap=colormaplist[i], interpolation='nearest', extent=[0, 106, 0, split_time], aspect='auto', alpha=0.2)
            plt.text(annotations[i][1], annotations[i][2], annotations[i][0], fontsize=8, color='white')

    plt.xlabel('Retention time 1 (min)')
    plt.ylabel('Retention time 2 (s)')
    plt.title(title)
    plt.savefig(savedir, dpi=1000, bbox_inches='tight')

def process_FID(filepath, split_time, sampling_interval,shift=0, solvent_time=0):
    df = parse_chromatogram(filepath)
    df = split_solvent(df, solvent_time)
    df = add_split(df,split_time,sampling_interval)
    df_array = convert_to2D(df,split_time)
    df_array  = baseline_stridewise(df_array)
    df_array = shift_phase(df_array, shift)
    df_array_norm = normalize_array(df_array)
    return df_array_norm







from scipy import sparse
from scipy.sparse.linalg import spsolve
def baseline_als(y, lam, p, niter=50):
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    D = lam * D.dot(D.transpose()) # Precompute this term since it does not depend on `w`
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    for i in range(niter):
        W.setdiag(w) # Do not create a new matrix, just update diagonal values
        Z = W + D
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z