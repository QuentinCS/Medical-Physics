#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 13:48:39 2024

@author: quentin
"""
#test
import matplotlib.pyplot as plt
import time
import class_cirs as c

start_time = time.time()

dir_to_save = './Results_CIRS'

# Creation of the different class's object
thorax_125 = c.cirs_analyse('./24-04-23_19-21-58_Thorax125/', 'Thorax_125_kV', dir_to_save)
thorax_rapide_125 = c.cirs_analyse('./24-04-23_19-03-24_ThoraxRapide125/', 'Thorax_rapide_125_kV', dir_to_save, ref=True, reference_coor=thorax_125.pixel_coords)
head_100 = c.cirs_analyse('./24-03-26_15-57-49_Tete100/', 'Head_100_kV', dir_to_save, ref=True, reference_coor=thorax_125.pixel_coords)
head_faible_100 = c.cirs_analyse('./24-03-26_16-08-04_TeteFaibleDose100/', 'Head_faible_100_kV', dir_to_save, ref=True, reference_coor=thorax_125.pixel_coords)
head_petit_100 = c.cirs_analyse('./24-04-23_19-25-58_Tete100PetitFantome/', 'Head_petit_fantome_100_kV', dir_to_save, ref=True, reference_coor=thorax_125.pixel_coords)
sein_125 = c.cirs_analyse('./24-04-23_19-09-53_Sein125/', 'Sein_125_kV', dir_to_save, ref=True, reference_coor=thorax_125.pixel_coords)
pelvis_125 = c.cirs_analyse('./24-03-26_16-10-02_Pelvis125/', 'Pelvis_125_kV', dir_to_save, ref=True, reference_coor=thorax_125.pixel_coords)
pelvis_rapide_125 = c.cirs_analyse('./24-04-23_15-56-44_PelvisRapide125/', 'Pelvis_rapide_125_kV', dir_to_save, ref=True, reference_coor=thorax_125.pixel_coords)
pelvis_140 = c.cirs_analyse('./24-03-26_16-26-21_PelvisLargeRapide140/', 'Pelvis_140_kV', dir_to_save, ref=True, reference_coor=thorax_125.pixel_coords)
dose_restreint_125 = c.cirs_analyse('./24-03-26_16-14-23_ImageDoseRestreinte80/', 'Dose_restreint_80_kV', dir_to_save, ref=True, reference_coor=thorax_125.pixel_coords)


# Plot of the HU curves 
plt.figure(figsize=(40, 20))
plt.style.use('ggplot')
plt.rc('font', size=40)
plt.plot(thorax_125.df['Density'], thorax_125.df[f'UH {thorax_125.protocol} {thorax_125.tension} kV'], label=f'{thorax_125.name}', marker='o', markersize=12, color='navy')
plt.plot(pelvis_125.df['Density'], pelvis_125.df[f'UH {pelvis_125.protocol} {pelvis_125.tension} kV'], label=f'{pelvis_125.name}', marker='o', markersize=12, color='red')
plt.plot(head_100.df['Density'], head_100.df[f'UH {head_100.protocol} {head_100.tension} kV'], label=f'{head_100.name}', marker='o', markersize=12, color='green')
plt.plot(head_faible_100.df['Density'], head_faible_100.df[f'UH {head_faible_100.protocol} {head_faible_100.tension} kV'], label=f'{head_faible_100.name}', marker='o', markersize=12, color='orange')
plt.plot(sein_125.df['Density'], sein_125.df[f'UH {sein_125.protocol} {sein_125.tension} kV'], label=f'{sein_125.name}', marker='o', markersize=12, color='lime')
plt.plot(pelvis_140.df['Density'], pelvis_140.df[f'UH {pelvis_140.protocol} {pelvis_140.tension} kV'], label=f'{pelvis_140.name}', marker='o', markersize=12, color='brown')
plt.plot(pelvis_rapide_125.df['Density'], pelvis_rapide_125.df[f'UH {pelvis_rapide_125.protocol} {pelvis_rapide_125.tension} kV'], label=f'{pelvis_rapide_125.name}', marker='o', markersize=12, color='grey')
plt.plot(thorax_rapide_125.df['Density'], thorax_rapide_125.df[f'UH {thorax_rapide_125.protocol} {thorax_rapide_125.tension} kV'], label=f'{thorax_rapide_125.name}', marker='o', markersize=12, color='pink')
plt.plot(dose_restreint_125.df['Density'], dose_restreint_125.df[f'UH {dose_restreint_125.protocol} {dose_restreint_125.tension} kV'], label=f'{dose_restreint_125.name}', marker='o', markersize=12, color='cyan')
plt.ylabel('Nombre Hounsfield (UH)')
plt.xlabel(r'Masse volumique ($g.cm^{-3}$)')
plt.legend()
plt.show()

# Plot of the HU curves 
plt.figure(figsize=(40, 20))
plt.style.use('ggplot')
plt.rc('font', size=40)
plt.plot(thorax_125.df['Density'], thorax_125.df[f'UH int {thorax_125.protocol} {thorax_125.tension} kV'], label=f'interne {thorax_125.protocol}', marker='o', markersize=12, color='navy')
plt.plot(thorax_125.df['Density'], thorax_125.df[f'UH ext {thorax_125.protocol} {thorax_125.tension} kV'], label=f'externe {thorax_125.protocol}', marker='o', markersize=12, color='red')
plt.plot(thorax_125.df['Density'], thorax_125.df[f'UH {thorax_125.protocol} {thorax_125.tension} kV'], label=f'total {thorax_125.protocol}', marker='o', markersize=12, color='green')
plt.ylabel('Nombre Hounsfield (UH)')
plt.xlabel(r'Masse volumique ($g.cm^{-3}$)')
plt.legend()
plt.show()

duree = time.time() - start_time
print ('\nTotal running time : %5.3g s' % duree)
