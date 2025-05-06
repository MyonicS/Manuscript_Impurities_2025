
import pandas as pd
import numpy as np
import os
path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Results', 'Results_PyIR.csv')

def write_db(catalyst,property,value,db_path=path):
    df = pd.read_csv(db_path, sep=',',header=0,index_col=0)
    df.loc[catalyst,property] = value
    df.to_csv(db_path, sep=',',header=True,index=True)

def read_db(catalyst,property,db_path=path):
    df = pd.read_csv(db_path, sep=',',header=0,index_col=0)
    return df.loc[catalyst,property]

def read_unit(property,db_path=path):
    df = pd.read_csv(db_path, sep=',',header=0,index_col=0)
    return df.loc['Unit',property]

def get_db(db_path=path):
    return pd.read_csv(db_path, sep=',',header=0,index_col=0)
    
