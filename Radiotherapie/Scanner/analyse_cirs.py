#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 13:48:39 2024

@author: quentin
"""

import matplotlib.pyplot as plt
import time
import class_cirs as c

start_time = time.time()

dir_to_save = './Results_CIRS/'

# Creation of the different class's object
thorax_125 = c.cirs_analyse('24-04-23_19-21-58_Thorax125/', 'Thorax 125 kV', dir_to_save)
pelvis_125 = c.cirs_analyse('24-03-26_16-10-02_Pelvis125/' , 'Pelvis 125 kV', dir_to_save, ref=True, reference_coor=thorax_125.pixel_coords)
head_100 = c.cirs_analyse('24-03-26_15-57-49_Tete100/' , 'Head 100 kV', dir_to_save, ref=True, reference_coor=thorax_125.pixel_coords)
sein_125 = c.cirs_analyse('24-04-23_19-09-53_Sein125/', 'Sein 125 kV', dir_to_save, ref=True, reference_coor=thorax_125.pixel_coords)

# Plot of the HU curves 
plt.figure(figsize=(40, 20))
plt.style.use('ggplot')
plt.rc('font', size=40)
plt.plot(thorax_125.df['Density'], thorax_125.df[f'UH {thorax_125.protocol} {thorax_125.tension} kV'], label=f'UH {thorax_125.name}', marker='o', markersize=12, color='navy')
plt.plot(pelvis_125.df['Density'], pelvis_125.df[f'UH {pelvis_125.protocol} {pelvis_125.tension} kV'], label=f'UH {pelvis_125.name}', marker='o', markersize=12, color='red')
plt.plot(head_100.df['Density'], head_100.df[f'UH {head_100.protocol} {head_100.tension} kV'], label=f'UH {head_100.name}', marker='o', markersize=12, color='green')
plt.plot(sein_125.df['Density'], sein_125.df[f'UH {sein_125.protocol} {sein_125.tension} kV'], label=f'UH {sein_125.name}', marker='o', markersize=12, color='orange')
plt.ylabel('Nombre Hounsfield (UH)')
plt.xlabel(r'Masse volumique ($g.cm^{-3}$)')
plt.legend()
plt.show()

# Plot of the HU curves 
plt.figure(figsize=(40, 20))
plt.style.use('ggplot')
plt.rc('font', size=40)
plt.plot(thorax_125.df['Density'], thorax_125.df[f'UH int {thorax_125.protocol} {thorax_125.tension} kV'], label=f'UH interne {thorax_125.protocol}', marker='o', markersize=12, color='navy')
plt.plot(thorax_125.df['Density'], thorax_125.df[f'UH ext {thorax_125.protocol} {thorax_125.tension} kV'], label=f'UH externe {thorax_125.protocol}', marker='o', markersize=12, color='red')
plt.plot(thorax_125.df['Density'], thorax_125.df[f'UH {thorax_125.protocol} {thorax_125.tension} kV'], label=f'UH total {thorax_125.protocol}', marker='o', markersize=12, color='green')
plt.ylabel('Nombre Hounsfield (UH)')
plt.xlabel(r'Masse volumique ($g.cm^{-3}$)')
plt.legend()
plt.show()

duree = time.time() - start_time
print ('\nTotal running time : %5.3g s' % duree)
