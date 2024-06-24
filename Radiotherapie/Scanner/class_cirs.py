#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 12:22:49 2024

@author: quentin
"""

from platform import python_version
import matplotlib.pyplot as plt
from datetime import date
from datetime import datetime
import SimpleITK as sitk
import pandas as pd
import pydicom as dcm
import numpy as np
import shutil
import cv2
import os

# Code based on class to analyse the data from CT (or CBCT) acquisition using the CIRS 062 phantom and plot the hounsfield unit- electronic density curve 
# The code use the relative positions of the insert to position the different ROIs based on the CIRS 062 geometry values obtained using ImageJ software
# It assume that the phantom is well positionned (especially the code doesn't correct for rotation of the phantom)
# The HU's mean values are obtained for each ROI using few slices (nb_slice) centered int the middle of the slices (important to have a centered acquisition)
# It is possible to display the positionned ROI on the image to check the correct position of the different ROIs

class cirs_analyse:
    def __init__(self, Dir, name, dir_to_save, nb_slice=4, size_insert=50, ref=False, reference_coor=False, plot=True):
        
        # Variables
        self.Dir = Dir # The directory to get the CT scans 
        self.name = name # Name can be used for plotting 
        self.plot = plot # Allow to plot directly when analyse the data
        self.dir_to_save = dir_to_save # Directory to save results and details
        self.ref = ref  # Allows to use reference pixel (manual or from another class object), need to fill reference_coords variable
        self.nb_slice = nb_slice # Number of slices to calculate mean ROIs
        self.reference_coor = reference_coor # Value of pixel position (manual or from another class object)
        self.size_insert = size_insert # diameter of insert in mm
        self.verbose = False 
        self.window_name = 'image'
        
        # Global variables to get the reference pixel position 
        self.pixel_coords = None
        self.pixel_value = None
        
        # Dataframe to get parameters values (Dicom and analysis parameters)
        self.data_info = None
        
        # Check if directory already exist, erase and create-the directory
        # WARNING remove directory if already exist
        if os.path.exists(self.dir_to_save):
            print('Directory already exist, data will be saved in former directory')
        else:
            print('Creating directory')
            os.makedirs(self.dir_to_save)
        if os.path.exists(f'{self.dir_to_save}/{self.name}'):
            shutil.rmtree(f'{self.dir_to_save}/{self.name}')
        os.makedirs(f'{self.dir_to_save}/{self.name}')
        if os.path.exists(f'{self.dir_to_save}/{self.name}/Slices'):
            shutil.rmtree(f'{self.dir_to_save}/{self.name}/Slices')
        os.makedirs(f'{self.dir_to_save}/{self.name}/Slices/')
        
        # Creation of a dataframe for the value and the analysis of the ROIs
        self.df = pd.DataFrame(list(zip(material, density)), columns = ['Material', 'Density'])
        # Analyse the images 
        self.analyse()
        if self.plot == True:
            self.plot_hu_ed()
        
        # Assemble data and save it
        self.resume_data()
        self.save()
               
    ############################### Function ###############################################      
            
    # Function to select a pixel using mouse (reference for the position of the differents ROIs)       
    def select_pixel(self, event, x, y, flags, param):
        if event == cv2.EVENT_LBUTTONDOWN:
            self.pixel_coords = (x, y)
            self.pixel_value = self.image[y, x]
            print(f"Selected Pixel: ({x}, {y})")
            print(f"Pixel value (Grayscale): {self.pixel_value}")
            self.is_pixel_selected = True
            cv2.destroyWindow(self.window_name)  # Close the windows after selecting pixel
             
    # Function to create a circular ROI (mask) from the center of the image or a position and a radius 
    def create_circular_mask(self, h, w, center=None, radius=None):
        if center is None: # use the middle of the image
            center = (int(w/2), int(h/2))
        if radius is None: # use the smallest distance between the center and image walls
            radius = min(center[0], center[1], w-center[0], h-center[1])

        Y, X = np.ogrid[:h, :w]
        dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

        mask = dist_from_center <= radius
        return mask

    # Main function to get the data from the Dicom and get different Tag Dicom values
    def analyse(self):        
        
        self.slice = 0
        files = []
        for dirName, subdirList, self.fileList in os.walk(self.Dir):
            files.append(self.fileList)
        self.fileList = [f for f in self.fileList if f.startswith('CT')]

        for i in range(0, len(self.fileList)):

            # Upload Dicom file and obtain the values of Dicom tags
            dicom_file = self.Dir + self.fileList[i]
            dicom_data = dcm.read_file(dicom_file)
            self.date = date(int(dicom_data[0x00080012].value[0:4]), int(dicom_data[0x00080012].value[4:6]), int(dicom_data[0x00080012].value[6:8]))
            self.time = datetime.strptime(dicom_data[0x00080013].value, '%H%M%S').time()
            self.device = dicom_data[0x00081090].value
            self.hospital = dicom_data[0x00080080].value
            self.position = dicom_data[0x00185100].value
            self.protocol = dicom_data[0x0008103E].value
            self.recons = dicom_data[0x00204000].value
            self.tension = int(dicom_data[0x00180060].value)
            self.charge = int(dicom_data[0x00181152].value)
            self.ctdi = dicom_data[0x00189345].value
            self.slice_thickness = dicom_data[0x00180050].value           
            self.pixel_size = float(dicom_data[0x00280030].value[0])
            self.size_row = dicom_data[0x00280010].value
            self.size_column = dicom_data[0x00280011].value
            self.roi_size = int(1*0.2*self.size_insert/self.pixel_size)
            self.central_slice = (int(len(self.fileList)/2))
            self.central_slice = 35
            serie = dicom_data[0x00200013].value
            
            # Loop on the number of slices needed for the computation, center the slices on the middle instance (median slice correct if the phantom is well centered)
            if serie >= (int(self.central_slice) - int(self.nb_slice/2)) and serie < (int(self.central_slice) + int(self.nb_slice/2)):
                # Extract image from Dicom
                self.image = dicom_data.pixel_array
                img_itk = sitk.ReadImage(self.Dir + self.fileList[i])
                image_np = sitk.GetArrayFromImage(img_itk)
         
                if self.slice == 0 and self.ref == False:
                    # Image normalisation to display image (optional)
                    image_visu = cv2.normalize(self.image, None, 0, 255, cv2.NORM_MINMAX)
                    image_visu = np.uint8(image_visu)
                    
                    # Create a windows and call the function to select pixel
                    cv2.namedWindow('image_visu')
                    cv2.setMouseCallback('image_visu', self.select_pixel)
                    print("Select center pixel of the central ROI and press 'ESC' when finish")

                    while True:
                        cv2.imshow('image_visu', image_visu)
                        if cv2.waitKey(1) & 0xFF == 27:  # Press 'ESC' to close window
                             break
                #  Use manual reference coordinate or from another object 
                if self.ref == True:
                    self.pixel_coords = self.reference_coor
                
                # Destroy all the remaining windows if necessary
                cv2.destroyAllWindows()
                
                # Print pixel's coordinate values and pixel value (useless)
                if self.verbose == True:
                    if self.pixel_coords is not None:
                        print(f"The selected pixel's coordinates: {self.pixel_coords}")
                        print(f"The selected pixel's value: {self.pixel_value}")
                    else:
                        print("No selected pixel")
                
                cv2.destroyAllWindows() # Remove all windows that can still be present 

                # Creation of all the ROI using the reference pixel's coordinates and the geometry's data of the CIRS 062 Phantom
                # Use the pixel spacing to adapt to different FOV or pixel size 
                roi_water = self.create_circular_mask(image_np.shape[1], image_np.shape[2], self.pixel_coords, self.roi_size)
                roi_lung_exhale_int = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_lung_exhale_int[0]/self.pixel_size, self.pixel_coords[1] + dist_lung_exhale_int[1]/self.pixel_size], self.roi_size)
                roi_lung_inhale_int = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_lung_inhale_int[0]/self.pixel_size, self.pixel_coords[1] + dist_lung_inhale_int[1]/self.pixel_size], self.roi_size)
                roi_lung_inhale_ext = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_lung_inhale_ext[0]/self.pixel_size, self.pixel_coords[1] + dist_lung_inhale_ext[1]/self.pixel_size], self.roi_size)
                roi_lung_exhale_ext = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_lung_exhale_ext[0]/self.pixel_size, self.pixel_coords[1] + dist_lung_exhale_ext[1]/self.pixel_size], self.roi_size)
                
                roi_bone200_int = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_bone200_int[0]/self.pixel_size, self.pixel_coords[1] + dist_bone200_int[1]/self.pixel_size], self.roi_size)
                roi_bone800_int = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_bone800_int[0]/self.pixel_size, self.pixel_coords[1] + dist_bone800_int[1]/self.pixel_size], self.roi_size)
                roi_bone200_ext = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_bone200_ext[0]/self.pixel_size, self.pixel_coords[1] + dist_bone200_ext[1]/self.pixel_size], self.roi_size)
                roi_bone800_ext = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_bone800_ext[0]/self.pixel_size, self.pixel_coords[1] + dist_bone800_ext[1]/self.pixel_size], self.roi_size)
                
                roi_muscle_int = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_muscle_int[0]/self.pixel_size, self.pixel_coords[1] + dist_muscle_int[1]/self.pixel_size], self.roi_size)
                roi_liver_int = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_liver_int[0]/self.pixel_size, self.pixel_coords[1] + dist_liver_int[1]/self.pixel_size], self.roi_size)
                roi_muscle_ext = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_muscle_ext[0]/self.pixel_size, self.pixel_coords[1] + dist_muscle_ext[1]/self.pixel_size], self.roi_size)
                roi_liver_ext = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_liver_ext[0]/self.pixel_size, self.pixel_coords[1] + dist_liver_ext[1]/self.pixel_size], self.roi_size)
                
                roi_adip_int = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_adip_int[0]/self.pixel_size, self.pixel_coords[1] + dist_adip_int[1]/self.pixel_size], self.roi_size)
                roi_breast_int = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_breast_int[0]/self.pixel_size, self.pixel_coords[1] + dist_breast_int[1]/self.pixel_size], self.roi_size)
                roi_adip_ext = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_adip_ext[0]/self.pixel_size, self.pixel_coords[1] + dist_adip_ext[1]/self.pixel_size], self.roi_size)
                roi_breast_ext = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_breast_ext[0]/self.pixel_size, self.pixel_coords[1] + dist_breast_ext[1]/self.pixel_size], self.roi_size)
                
                #roi_pmma = self.create_circular_mask(image_np.shape[1], image_np.shape[2], [self.pixel_coords[0] + dist_pmma[0]/self.pixel_size, self.pixel_coords[1] + dist_pmma[1]/self.pixel_size], self.roi_size)
                
                # Extraction of the different ROI using mask (zero value outside the ROIs)
                img_roi_water = image_np[0].copy()
                img_roi_water[~roi_water] = 0
                img_roi_lung_exhale_int = image_np[0].copy()
                img_roi_lung_exhale_int[~roi_lung_exhale_int] = 0
                img_roi_lung_exhale_ext = image_np[0].copy()
                img_roi_lung_exhale_ext[~roi_lung_exhale_ext] = 0
                img_roi_lung_inhale_int = image_np[0].copy()
                img_roi_lung_inhale_int[~roi_lung_inhale_int] = 0
                img_roi_lung_inhale_ext = image_np[0].copy()
                img_roi_lung_inhale_ext[~roi_lung_inhale_ext] = 0 
                
                img_roi_bone200_ext = image_np[0].copy()
                img_roi_bone200_ext[~roi_bone200_ext] = 0
                img_roi_bone800_ext = image_np[0].copy()
                img_roi_bone800_ext[~roi_bone800_ext] = 0
                img_roi_bone200_int = image_np[0].copy()
                img_roi_bone200_int[~roi_bone200_int] = 0
                img_roi_bone800_int = image_np[0].copy()
                img_roi_bone800_int[~roi_bone800_int] = 0
                
                img_roi_muscle_int = image_np[0].copy()
                img_roi_muscle_int[~roi_muscle_int] = 0
                img_roi_liver_int = image_np[0].copy()
                img_roi_liver_int[~roi_liver_int] = 0
                img_roi_muscle_ext = image_np[0].copy()
                img_roi_muscle_ext[~roi_muscle_ext] = 0
                img_roi_liver_ext = image_np[0].copy()
                img_roi_liver_ext[~roi_liver_ext] = 0
                
                img_roi_adip_int = image_np[0].copy()
                img_roi_adip_int[~roi_adip_int] = 0
                img_roi_breast_int = image_np[0].copy()
                img_roi_breast_int[~roi_breast_int] = 0
                img_roi_adip_ext = image_np[0].copy()
                img_roi_adip_ext[~roi_adip_ext] = 0
                img_roi_breast_ext = image_np[0].copy()
                img_roi_breast_ext[~roi_breast_ext] = 0
         
                #img_roi_pmma = image_np[0].copy()
                #img_roi_pmma[~roi_pmma] = 0

                # Plot the image and the different superimposed to check the correct position (by default activate on the first slice)
                img_roi_insert = img_roi_water + img_roi_lung_exhale_ext + img_roi_lung_exhale_int + img_roi_lung_inhale_ext + img_roi_lung_inhale_int + img_roi_bone200_ext + img_roi_bone200_int + img_roi_bone800_ext + img_roi_bone800_int + img_roi_muscle_ext + img_roi_muscle_int + img_roi_liver_ext + img_roi_liver_int + img_roi_breast_ext + img_roi_breast_int + img_roi_adip_ext + img_roi_adip_int# + img_roi_pmma
                alpha = 0.3
                o = cv2.addWeighted(image_np[0], alpha, img_roi_insert, 1 - alpha, 0)
                plt.figure(figsize=(40, 25))
                plt.imshow(o, cmap='gray')
                plt.title(f'{self.name}')
                plt.colorbar()
                plt.savefig(f'{self.dir_to_save}/{self.name}/verif_center_ROI.png')
                if self.slice != 0:
                    plt.close()
                    plt.savefig(f'{self.dir_to_save}/{self.name}/Slices/verif_ROI_slice_{serie}.png')

                # Save the value of mean UH for the different ROIs
                # Do for the inner ROIs and the outer ROIs separately
                self.df.loc[:,f'Coupe int {self.slice}'] = [np.sum(img_roi_bone800_int)/(np.pi*self.roi_size**2), np.sum(img_roi_bone200_int)/(np.pi*self.roi_size**2), np.sum(img_roi_muscle_int)/(np.pi*self.roi_size**2), np.sum(img_roi_liver_int)/(np.pi*self.roi_size**2), np.sum(img_roi_breast_int)/(np.pi*self.roi_size**2), np.sum(img_roi_adip_int)/(np.pi*self.roi_size**2), np.sum(img_roi_lung_inhale_int)/(np.pi*self.roi_size**2), np.sum(img_roi_lung_exhale_int)/(np.pi*self.roi_size**2), np.sum(img_roi_water)/(np.pi*self.roi_size**2)]#, np.sum(img_roi_pmma)/(np.pi*self.roi_size**2)]
                self.df.loc[:,f'Coupe ext {self.slice}'] = [np.sum(img_roi_bone800_ext)/(np.pi*self.roi_size**2), np.sum(img_roi_bone200_ext)/(np.pi*self.roi_size**2), np.sum(img_roi_muscle_ext)/(np.pi*self.roi_size**2), np.sum(img_roi_liver_ext)/(np.pi*self.roi_size**2), np.sum(img_roi_breast_ext)/(np.pi*self.roi_size**2), np.sum(img_roi_adip_ext)/(np.pi*self.roi_size**2), np.sum(img_roi_lung_inhale_ext)/(np.pi*self.roi_size**2), np.sum(img_roi_lung_exhale_ext)/(np.pi*self.roi_size**2), np.sum(img_roi_water)/(np.pi*self.roi_size**2)]#, np.sum(img_roi_pmma)/(np.pi*self.roi_size**2)]
                self.df.loc[:,f'Coupe STD int {self.slice}'] = [self.stdev_roi(img_roi_bone800_int), self.stdev_roi(img_roi_bone200_int), self.stdev_roi(img_roi_muscle_int), self.stdev_roi(img_roi_liver_int), self.stdev_roi(img_roi_breast_int), self.stdev_roi(img_roi_adip_int), self.stdev_roi(img_roi_lung_inhale_int), self.stdev_roi(img_roi_lung_exhale_int), self.stdev_roi(img_roi_water)]#, self.stdev_roi(img_roi_pmma)]
                self.df.loc[:,f'Coupe STD ext {self.slice}'] = [self.stdev_roi(img_roi_bone800_ext), self.stdev_roi(img_roi_bone200_ext), self.stdev_roi(img_roi_muscle_ext), self.stdev_roi(img_roi_liver_ext), self.stdev_roi(img_roi_breast_ext), self.stdev_roi(img_roi_adip_ext), self.stdev_roi(img_roi_lung_inhale_ext), self.stdev_roi(img_roi_lung_exhale_ext), self.stdev_roi(img_roi_water)]#, self.stdev_roi(img_roi_pmma)]
                self.slice += 1

        # Compute the mean on all the selected slice for the mean HU values and save in the dataframe
        self.df.loc[:,f'UH int {self.protocol} {self.tension} kV'] = sum([self.df[f'Coupe int {i}'] for i in range(self.nb_slice)]) / self.nb_slice
        self.df.loc[:,f'UH ext {self.protocol} {self.tension} kV'] = sum([self.df[f'Coupe ext {i}'] for i in range(self.nb_slice)]) / self.nb_slice
        self.df.loc[:,f'UH {self.protocol} {self.tension} kV'] = (self.df[f'UH int {self.protocol} {self.tension} kV'] + self.df[f'UH ext {self.protocol} {self.tension} kV'])/2
        self.df.loc[:,f'STD int {self.protocol} {self.tension} kV'] = sum([self.df[f'Coupe STD int {i}'] for i in range(self.nb_slice)]) /self.nb_slice
        self.df.loc[:,f'STD ext {self.protocol} {self.tension} kV'] = sum([self.df[f'Coupe STD ext {i}'] for i in range(self.nb_slice)]) /self.nb_slice
        self.df.loc[:,f'STD {self.protocol} {self.tension} kV'] = (self.df[f'STD int {self.protocol} {self.tension} kV'] + self.df[f'STD ext {self.protocol} {self.tension} kV'])/2
        self.df = self.df.sort_values(by=['Density'], ascending=False)

    # Function for plot the HU curve as function of the electronic densities (disable by default)
    def plot_hu_ed(self): 
        plt.figure(figsize=(40, 20))
        plt.style.use('ggplot')
        plt.rc('font', size=40)
        plt.plot(self.df['Density'][:], self.df[f'UH {self.protocol} {self.tension} kV'], label=f'UH {self.protocol} {self.tension} kV', marker='o', markersize=12, color='navy')
        plt.text(1.2, -600, f'{self.device}\n {self.tension} kV\n {self.charge} mAs\n {self.position}\n Slice: {self.slice_thickness} mm\n Pixel: {self.pixel_size:.3f} mm\n {self.size_row} x {self.size_column}', fontsize=30, color='gray')
        plt.ylabel('Hounsfield unit (HU)')
        plt.xlabel(r'Density ($g.cm^{-3}$)')
        plt.legend()
        plt.savefig(f'{self.dir_to_save}/{self.name}/HU.png')
        
    # Function to calculate the noise in the ROIs     
    def stdev_roi(self, img):
        array = np.array(img)
        image_arr = array[array != 0] # Remove zeros values (outside roi)
        return np.std(image_arr)
    
    # Function to assemble data from Dicom and analysis to save into dataframe
    def resume_data(self):
        data = {
                'Parameter': ['Version', 'Device', 'Hospital', 'Date', 'Time', 'Protocol', 'Tension (kV)', 'Charge (mAs)','CTDI (mGy)', 'Slice thickness (mm)', 'Pixel size (mm)', 'Matrix', 'Patient position', 'Reconstruction', 'Nb slices used', 'Central Slice'],
                'Value': [f'Python {python_version()}', self.device, self.hospital, self.date, self.time, self.protocol, self.tension, self.charge, self.ctdi, self.slice_thickness, self.pixel_size, f'{self.size_row} x {self.size_column}', self.position, self.recons, self.nb_slice, self.central_slice]
                }
        self.data_info = pd.DataFrame(data, columns=['Parameter', 'Value'])
        
    # Function to save to excel file results and data parameters
    def save(self):
        with pd.ExcelWriter(f'{self.dir_to_save}/{self.name}/Data CIRS {self.name}.xlsx') as file:
            self.data_info.to_excel(file, sheet_name='Parameters')
            self.df.to_excel(file, sheet_name='Results', columns=['Material', 'Density', f'UH {self.protocol} {self.tension} kV', f'STD {self.protocol} {self.tension} kV'])
                   
################################# Data ###############################################

# CIRS's Insert data for the different materials 

material = ['Bone 800', 'Bone 200', 'Muscle', 'Liver', 'Breast', 'Adipose', 'Lung inhale', 'Lung exhale', 'Water']#, 'PMMA']
materiaux = ['Os 800', 'Os 200', 'Muscle', 'Foie', 'Sein', 'Tissus adipeux', 'Poumon inspi', 'Poumon expi', 'Eau']#, 'PMMA']
density = [1.53, 1.16, 1.06, 1.07, 0.99, 0.96, 0.2, 0.5, 1]#, 1.029]

# Geometric distance of the different ROI from the central ROI
# The value can be different in case of different insert 
# The ROIs are positionned: (clockwise from the top)
# Head insert: Bone 200 int - Liver int - Lung exhale int - Adipose int - Bone 800 int - Muscle int - Lung inhale int - Breast int
# Body insert: Bone 800 ext - Muscle ext - Lung inhale ext - Breast ext - Bone 200 ext - Liver ext - Lung exhale ext - Adipose ext

# Distances in mm 
R1 = 62
R2 = 115
R3 = 60
R4 = 111

dist_lung_exhale_int = [R1, 0]
dist_lung_inhale_ext = [R2, 0]
dist_lung_exhale_ext = [-R2, 0]
dist_lung_inhale_int = [-R1, 0]

dist_bone200_int = [0, -R3]
dist_bone800_ext = [0, -R4-5]
dist_bone200_ext = [0, R4]
dist_bone800_int = [0, R3]
   
dist_liver_int = [R1*np.cos(np.radians(45)), -R3*np.sin(np.radians(45))]
dist_liver_ext = [-R2*np.cos(np.radians(45)), R4*np.sin(np.radians(45))]
dist_muscle_ext = [R2*np.cos(np.radians(45)), -R4*np.sin(np.radians(45))-3]
dist_muscle_int = [-R1*np.cos(np.radians(45)), R3*np.sin(np.radians(45))]  
    
dist_adip_int = [R1*np.cos(np.radians(45)), R3*np.sin(np.radians(45))]
dist_adip_ext = [-R2*np.cos(np.radians(45)), -R4*np.sin(np.radians(45))]
dist_breast_ext = [R2*np.cos(np.radians(45)), R4*np.sin(np.radians(45))]
dist_breast_int = [-R1*np.cos(np.radians(45)), -R3*np.sin(np.radians(45))]

dist_pmma = [R2*np.cos(np.radians(25)), R2*np.sin(np.radians(25))]
