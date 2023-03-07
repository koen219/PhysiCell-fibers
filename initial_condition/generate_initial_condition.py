# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 17:53:43 2023

@author: s169389
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import os

def line(x, a, b):
    return a*x+b


x = np.linspace(0,10, 5)

y = line(x,1,1)
z = np.zeros(len(x))

plt.scatter(x,y)

rows = []

for xi, yi, zi in zip(x,y,z):       
    rows.append([xi, yi, zi, 0])
    
full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)


with open(path+'\intial_cells.csv', 'w', newline='') as f:
    # create the csv writer
    writer = csv.writer(f)

    for row in rows:
        writer.writerow(row)

    
f.close()