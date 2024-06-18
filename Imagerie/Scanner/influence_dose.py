#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 14:41:01 2023

@author: quentin
"""

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import numpy as np
import pandas as pd
import time

# Fonction pour le fit de la dose en fonction des paramètres
def kv_fit(d, a, b):
    return a*np.power(d, b)

def mas_fit(d, a, b):
    return (a*d)+b

start_time = time.time()

######################################################################
############################# Tension ################################
######################################################################

tension_prime_sp = {"Tension (kV)": [80, 100, 120, 135],
           "CTDI mes (mGy)": [1.161, 2.775, 4.958, 7.327],
           "CTDI aff (mGy)": [1 ,3 , 5.4, 7.7],
            }

tension_one_prism = {"Tension (kV)": [80, 100, 120, 135],
           "CTDI mes (mGy)": [1.569, 3.440, 6.010, 8.249],
           "CTDI aff (mGy)": [1.5, 3.3, 5.8, 8],
            }


parameters_kV_prime_sp, cov = curve_fit(kv_fit, tension_prime_sp['Tension (kV)'], tension_prime_sp['CTDI mes (mGy)'])
fit_kv_a_prime_sp = parameters_kV_prime_sp[0]
fit_kv_b_prime_sp = parameters_kV_prime_sp[1]
fit_prime_sp = kv_fit(tension_prime_sp['Tension (kV)'], fit_kv_a_prime_sp, fit_kv_b_prime_sp)
Tension = np.linspace(80, 140, 100)
fit1_prime_sp = kv_fit(Tension, fit_kv_a_prime_sp, fit_kv_b_prime_sp)

parameters_kV_one_prism, cov = curve_fit(kv_fit, tension_one_prism['Tension (kV)'], tension_one_prism['CTDI mes (mGy)'])
fit_kv_a_one_prism = parameters_kV_one_prism[0]
fit_kv_b_one_prism = parameters_kV_one_prism[1]
fit_one_prism = kv_fit(tension_one_prism['Tension (kV)'], fit_kv_a_one_prism, fit_kv_b_one_prism)
Tension = np.linspace(80, 140, 100)
fit1_one_prism = kv_fit(Tension, fit_kv_a_one_prism, fit_kv_b_one_prism)


plt.figure(figsize=(20, 10))
#plt.scatter(tension_prime_sp["Tension (kV)"], tension_prime_sp['CTDI aff (mGy)'], label='Affiché', color='red')
plt.scatter(tension_prime_sp["Tension (kV)"], tension_prime_sp['CTDI mes (mGy)'], label='Aquilion Prime SP', color='navy')
plt.plot(Tension, fit1_prime_sp, color='blue')#, label='Fit: CTDI = %.2e $kV^{%.3f}$\n $r²$ = %.4f'%(fit_kv_a_prime_sp, fit_kv_b_prime_sp, r2_score(tension_prime_sp['CTDI mes (mGy)'], fit_prime_sp)))
plt.text(120, 3.5, 'Fit: CTDI = %.2e $kV^{%.3f}$\n $r²$ = %.4f'%(fit_kv_a_prime_sp, fit_kv_b_prime_sp, r2_score(tension_prime_sp['CTDI mes (mGy)'], fit_prime_sp)), fontsize=14, color='blue')
plt.scatter(tension_one_prism["Tension (kV)"], tension_one_prism['CTDI mes (mGy)'], label='Aquilion One Prism', color='green')
plt.plot(Tension, fit1_one_prism, color='green')#, label='Fit: CTDI = %.2e $kV^{%.3f}$\n $r²$ = %.4f'%(fit_kv_a_one_prism, fit_kv_b_one_prism, r2_score(tension_one_prism['CTDI mes 1 (mGy)'], fit_one_prism)))
plt.text(110, 7, 'Fit: CTDI = %.2e $kV^{%.3f}$\n $r²$ = %.4f'%(fit_kv_a_one_prism, fit_kv_b_one_prism, r2_score(tension_one_prism['CTDI mes (mGy)'], fit_one_prism)), fontsize=14, color='green')
plt.text(125, 1, f'Canon,\n Axial\n 100 mAs\n Fantôme:32 cm\n 8mm x 4', fontsize=14)
plt.xlabel('Tension (kV)', fontsize=20)
plt.ylabel('CTDI vol (mGy)', fontsize=20)
plt.legend(fontsize=20)


######################################################################
############################## Charge ################################
######################################################################

charge = {"Charge (mAs)": [50, 100, 150, 200, 250, 300, 350, 400, 450, 500],
           "CTDI mes (mGy)": [1.423, 2.768, 4.167, 5.681, 7.157, 8.228, 9.562, 11.013, 13.723, 15.075],
           "CTDI aff (mGy)": [1.5, 3, 4.5, 6, 7.6, 9.1, 10.6, 12.1, 14.4, 16],
            }

parameters_mas, cov = curve_fit(mas_fit, charge['Charge (mAs)'], charge['CTDI mes (mGy)'])
fit_mas_a = parameters_mas[0]
fit_mas_b = parameters_mas[1]
fit_mas = mas_fit(np.array(charge['Charge (mAs)']), fit_mas_a, fit_mas_b)
Charge = np.linspace(50, 500, 100)
fit_mas1 = mas_fit(Charge, fit_mas_a, fit_mas_b)


plt.figure(figsize=(20, 10))
plt.scatter(charge["Charge (mAs)"], charge['CTDI mes (mGy)'], label='Mesuré', color='navy')
plt.scatter(charge["Charge (mAs)"], charge['CTDI aff (mGy)'], label='Affiché', color='red')
plt.plot(Charge, fit_mas1, color='blue', label='Fit: CTDI = %.2e mAs + %.3f\n $r²$ = %.4f'%(fit_mas_a, fit_mas_b, r2_score(charge['CTDI mes (mGy)'], fit_mas)))
plt.text(400, 1, f'Canon,\n Acquilion Prime SP\n Axial\n 100 kV\n Fantôme:32 cm\n 8mm x 4', fontsize=14)
plt.xlabel('Charge (mAs)', fontsize=20)
plt.ylabel('CTDI vol (mGy)', fontsize=20)
plt.legend(fontsize=20)

######################################################################
########################### Collimation ##############################
######################################################################

collimation = {"Collimation (mm)": ['0.5 x 4', '1 x 4', '2 x 4', '3 x 4', '4 x 4', '5 x 4', '8 x 4'],
           "CTDI mes (mGy)": [9.137, 5.643, 4.035, 3.406, 3.132, 3.192, 2.775],
           "CTDI aff (mGy)": [10.1, 6.3, 4.5, 3.9, 3.5, 3.6, 3],
            }

plt.figure(figsize=(20, 10))
plt.scatter(collimation["Collimation (mm)"], collimation['CTDI mes (mGy)'], label='Mesuré', color='navy')
plt.scatter(collimation["Collimation (mm)"], collimation['CTDI aff (mGy)'], label='Affiché', color='red')
plt.text(0.2, 3, f'Canon,\n Acquilion Prime SP\n Axial\n 100 kV\n 100 mAs\n Fantôme:32 cm', fontsize=14)
plt.xlabel('Collimation (mm)', fontsize=20)
plt.ylabel('CTDI vol (mGy)', fontsize=20)
plt.legend(fontsize=20)


######################################################################
############################### Pitch ################################
######################################################################

pitch = {"Pitch (-)": [0.625, 0.825, 1.250, 1.500],
         "Pitch 2 (-)": [0.625, 0.825, 1.155, 1.250, 1.375, 1.500, 1.575],
           "CTDI aff (mGy)": [18.1, 13.7, 9.1, 7.5],
           "CTDI mes (mGy)": [17.565, 13.293, 8.774, 7.169],
           "mAs effectives": [160, 121, 86, 80, 72, 66, 63],
            }

plt.figure(figsize=(20, 10))
plt.scatter(pitch["Pitch (-)"], pitch['CTDI mes (mGy)'], label='Mesuré (hélicoidal)', color='navy')
plt.scatter(pitch["Pitch (-)"], pitch['CTDI aff (mGy)'], label='Affiché (axial)', color='red')
plt.text(0.65, 1, f'Canon,\n Acquilion Prime SP\n Hélicoidal\n 100 kV\n Fantôme:32 cm\n 1.0mm x 40', fontsize=14)
plt.xlabel('Pitch (-)', fontsize=20)
plt.ylabel('CTDI vol (mGy)', fontsize=20)
plt.ylim(0)
plt.legend(fontsize=20)

parameters_pitch, cov = curve_fit(kv_fit, pitch['Pitch 2 (-)'], pitch['mAs effectives'])
fit_pitch_a = parameters_pitch[0]
fit_pitch_b = parameters_pitch[1]
fit_pitch = kv_fit(np.array(pitch['Pitch 2 (-)']), fit_pitch_a, fit_pitch_b)
Pitch = np.linspace(0.5, 1.6, 100)
fit_pitch1 = kv_fit(Pitch, fit_pitch_a, fit_pitch_b)

p = []
for i in range(len(pitch['Pitch 2 (-)'])):
    p.append(pitch['mAs effectives'][i]*pitch['Pitch 2 (-)'][i])
pitch['Charge (mAs)'] = p


plt.figure(figsize=(20, 10))
plt.scatter(pitch["Pitch 2 (-)"], pitch['mAs effectives'], label='mAs effectives', color='navy')
plt.scatter(pitch["Pitch 2 (-)"], pitch['Charge (mAs)'], label='mAs', color='red')
plt.plot(Pitch, fit_pitch1, color='blue', label='Fit: mAs = %.2f $PF^{%.3f}$\n $r²$ = %.4f'%(fit_pitch_a, fit_pitch_b, r2_score(pitch['mAs effectives'], fit_pitch)))
plt.text(0.5, 60, f'Canon,\n Acquilion Prime SP\n Hélicoidal\n 100 kV\n Fantôme:32 cm\n 1.0mm x 40', fontsize=14)
plt.xlabel('Pitch (-)', fontsize=20)
plt.ylabel('Charge (mAs)', fontsize=20)
plt.legend(fontsize=20)



######################################################################
######################### Dose radiale ###############################
######################################################################

dose = {"Tension (kV)": [80, 100, 120, 135],
           "Rapport péri/centre": [2.0304, 1.9236, 1.8239, 1.9056],
            }

plt.figure(figsize=(20, 10))
plt.scatter(dose["Tension (kV)"], dose['Rapport péri/centre'], label='Fantôme 32', color='navy')
plt.text(80, 1.85, f'Canon,\n Acquilion Prime SP\n Axial\n 100 mAs\n Fantôme:32 cm\n 1.0mm x 40', fontsize=14)
plt.xlabel('Tension (kV)', fontsize=20)
plt.ylabel('Rapport de dose périphérique/centre (-)', fontsize=20)
#plt.ylim(0)
plt.legend(fontsize=20)


######################################################################
######################### Dose hauteur ###############################
######################################################################

dose_h = {"Hauteur iso (cm)": [-2, 0, 2],
           "CTDI (mGy)": [3.339, 3.553, 3.419],
           "Écart rel (%)": [-7.26, -1.3, -5.02],
            }

plt.figure(figsize=(20, 10))
plt.scatter(dose_h["Hauteur iso (cm)"], dose_h['Écart rel (%)'], label='Fantôme 32', color='navy')
plt.text(-2, -4, f'Canon,\n Acquilion Prime SP\n Axial\n 100 kV\n 100 mAs\n Fantôme:32 cm\n 1.0mm x 40', fontsize=14)
plt.xlabel('Décalage en hauteur isocentre (cm)', fontsize=20)
plt.ylabel('Écart relatif CTDI affiché (%)', fontsize=20)
plt.ylim(-10, 0)
#plt.legend(fontsize=20)






duree = time.time() - start_time
print ('\n \nTotal running time : %5.3g s' % duree)
