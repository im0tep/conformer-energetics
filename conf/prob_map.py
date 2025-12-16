#!/usr/bin/env python3
import argparse
parser=argparse.ArgumentParser() 
parser.add_argument('--pattern', default='ff.*.lammpstrj', help= 'Use only .lammpstrj files matching this pattern (default: \'ff.*lammpstrj\')')
args= parser.parse_args()

import numpy as np
import pandas as pd
from glob import glob 
import matplotlib as mpl
from matplotlib import pyplot as plt 
import os 

def compute_dihedral(a0, a1, a2, a3):
    b1 = a1 - a0
    b2 = a2 - a1
    b3 = a3 - a2

    n1 = np.cross(b1, b2)
    n1 = n1 / np.linalg.norm(n1, axis=1)[:, None]

    n2 = np.cross(b2, b3)
    n2 = n2 / np.linalg.norm(n2, axis=1)[:, None]

    b2n = b2 / np.linalg.norm(b2, axis=1)[:, None]
    m1 = np.cross(n1, b2n)

    x = np.sum(n1 * n2, axis=1)
    y = np.sum(m1 * n2, axis=1)

    return abs(np.degrees(np.arctan2(y, x)))

def compute_distance(p1, p2):
    return np.linalg.norm(p1 - p2, axis=1)

dihedrals_all = []
distances_all = []
transitions= []

files = glob(args.pattern, recursive= True)
for f in files:
    traj = pd.read_csv (f, sep = '\s+', names= ['timestep','ATOMS','type','x','y','z'])
    if len(traj) < 60006:
        print (f, ' trajectory is too short to compute')
        continue #raise(RuntimeError('Too few lines to compute'))
    C1, C2, O1, O2, H1, H2 = ( 
        traj [traj['ATOMS'] ==2],
        traj [traj['ATOMS'] ==1],
        traj [traj['ATOMS'] ==5],
        traj [traj['ATOMS'] ==8],
        traj [traj['ATOMS'] ==13],  
        traj [traj['ATOMS'] ==12]
    ) 
    dihedrals= compute_dihedral(
        O1.iloc[:, -3:].values,
        C1.iloc[:, -3:].values,
        C2.iloc[:, -3:].values,
        O2.iloc[:, -3:].values
    )
    distances= compute_distance(H1.iloc[:, -3:].values , H2.iloc[:, -3:].values)
    cond= dihedrals > 120
    transition = np.count_nonzero(cond[1:] != cond[:-1])
    print ('number of transitions for', f ,' is: ', transition)
    dihedrals_all.append (dihedrals)
    distances_all.append (distances)
    transitions.append (transition)

print ('Std dev in transitions for', len (files),'trajectories is', np.std(transitions)) 

dihedrals_all_flat= np.concatenate(dihedrals_all).ravel()
distances_all_flat= np.concatenate(distances_all).ravel()

basedir = os.getcwd()

df = pd.DataFrame({
    'Dihedral': dihedrals_all_flat,
    'Distance': distances_all_flat
})

df.to_csv(basedir + '/' +'dihedrals_distances.csv', index=False)

gauche = dihedrals_all_flat[ (dihedrals_all_flat < 120  ) ]   
trans = dihedrals_all_flat[ (dihedrals_all_flat >= 120  ) ] 

print ( 'Percentage of gauche conformers in analysed trajectories is ' + str (round (len (gauche) /len (dihedrals_all_flat),2)) +  '%, percentage of trans conformers is ' + str (round (len (trans) /len (dihedrals_all_flat), 2)) + '%' )

def block_avg(data, block_size):
    data = data[:len(data) // block_size * block_size]
    return [np.mean(data[i:i+block_size]) for i in range(0, len(data), block_size)]

gauche_rand = np.random.permutation(gauche)
trans_rand = np.random.permutation(trans)

n_blocks = 100
gauche_avg = block_avg(gauche_rand, n_blocks)
trans_avg = block_avg(trans_rand, n_blocks)

gauche_err = np.std(gauche_avg) / np.sqrt(len(gauche_avg))
trans_err = np.std(trans_avg) / np.sqrt(len(trans_avg))

print ('Std Dev of the mean for gauche is ' + str (gauche_err))
print ('Std Dev of the mean for trans is ' + str (trans_err))

# quick plot but will move plotting to another script 
plt.hist2d(distances_all_flat, dihedrals_all_flat, bins=30, density=True, cmap='turbo')
plt.xlabel('H-H Distance (Å)')
plt.ylabel('OCCO Dihedral (°)')
plt.colorbar(label='Probability Density')
plt.xlim([2,5.5])
plt.savefig (basedir +'/'+'prob_map.pdf', bbox_inches='tight')
