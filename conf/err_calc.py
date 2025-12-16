#!/bin/usr/env python3
import pandas as pd 
import numpy as np 

df= pd.read_csv ('dihedrals_distances.csv' )

def block_avg(data, block_size):
    data = data[:len(data) // block_size * block_size]
    return [sum(data[i:i+block_size]) / block_size for i in range(0, len(data), block_size)]

gauche = df [df ['Dihedral' ] <= 120]
trans = df [df ['Dihedral' ] > 120]

gauche_rand =  gauche.sample(frac=1).reset_index(drop=True)
trans_rand =  trans.sample(frac=1).reset_index(drop=True)

n_blocks= 100 

gauche_avg= block_avg(gauche_rand['Dihedral'], 100)
trans_avg= block_avg(trans_rand['Dihedral'], 100)

gauche_err = np.std(gauche_avg) / np.sqrt(len(gauche_avg))

trans_err= np.std(trans_avg) / np.sqrt(len(trans_avg))

print ('Std Dev of the mean for gauche is ' + str (gauche_err)) 
print ('Std Dev of the mean for trans is ' + str (trans_err)) 



