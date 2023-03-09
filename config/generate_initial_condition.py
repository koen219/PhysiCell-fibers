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


x = np.random.randint(-50,50, size=10)

y = np.random.randint(-50,50,size=10)
z = np.zeros(len(x))

plt.scatter(x,y)

rows = []

for xi, yi, zi in zip(x,y,z):       
    rows.append([xi, yi, zi, 0])
    
# number_of_pushers = 0 
 
# while number_of_pushers < 100:
    
    # xi = np.random.randint(-500,500)

    # yi = np.random.randint(-500,500)
    
    # if 100**2 < xi**2 + yi**2 < 150**2:
    #     rows.append([xi, yi, 0,1])
    #     number_of_pushers += 1
    
    #     plt.scatter(xi,yi, color='black')
    
full_path = os.path.realpath(__file__)
path, filename = os.path.split(full_path)


with open(path+'\cells.csv', 'w', newline='') as f:
    # create the csv writer
    writer = csv.writer(f)

    for row in rows:
        writer.writerow(row)

    
f.close()