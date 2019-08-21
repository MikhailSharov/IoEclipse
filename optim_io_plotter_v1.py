# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 08:10:37 2019

@author: msharov
"""

#importing any neccessary functions (some are not required and are from previous revs)
import numpy as np
from astropy import units as u
from astropy.nddata import CCDData
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.time import Time
from scipy.ndimage.interpolation import shift
from astropy.modeling.models import Lorentz1D
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
import scipy.optimize

def f(factor, frame):
    #defining any lists and variables
    x = []
    y = []
    multfactor = factor
    scatter = []
    science = []
    counts = [100884,108497,114845,122934,134215,145975]
    
    #importing all files including scatter and science and stored in defined lists
    for i in range(0,6):  
        scatter_file = fits.open('fullspecJupiter_center_Spectrum.0001.ec.fits')            #importing scatter frames into list scatter
        scatter_data = scatter_file[0].data
        scatter.append(scatter_data)
        scatter_file.close()
    
    for i in range(1,7):
        image_file = fits.open('fullspecIo_eclipsed.000'+str(i)+'.ec.fits')     #importing science frames into list science
        image_data = image_file[0].data
        science.append(image_data)
        image_file.close()
    
    #unnormalizing the data and subtracting scatter
    for i in range(0,6):
        science[i] = science[i]*counts[i]    

    for i,j in zip(range(0,6), [1.5,1.75,2,2.25,2.5,2.75]):
        scatter[i] = shift(scatter[i], j)

    box_kernel = Box1DKernel(2)
    for i in range(0,6):
        scatter[i] = convolve(scatter[i], box_kernel)                                                   #unnormalization of science to counts
    
    y.append([science[frame][j]-(multfactor*(scatter[frame][j])) for j in range(57733,57938)])         #subraction for science frames 1-4

    for i in range(57733,57938):
        x.append(3359.29052734 + i*0.043806734019)                                              #x axis in angstroms
    
    summ = sum(map(abs, y[0]))
    return(summ)

multfactor = []    
for i in range(0,6):
    res = scipy.optimize.minimize(f, x0 = 10000, args = i, method = 'SLSQP', options = {'maxiter':10000, 'ftol': 1e-04})
    multfactor.append(res.x[0])



#defining any lists and variables
xsci = []
xsca = []
x = []
y = []
y1 = []
#multfactor = [101749,109474.1,114828.6,123332,134674.2,146222.7]                              #minimization factors
multfactor = [100402.03106618,108053.45309066,113607.95109659,121980.62509868,133171.22724197,144453.97572916]  
scatter = []
science = []
counts = [100884,108497,114845,122934,134215,145975]                                        #the counts from the box on science frams for unnormalization
centercount = 204525  

for i in range(1,7):
    image_file = fits.open('fullspecIo_eclipsed.000'+str(i)+'.ec.fits')                     #importing science frames into list science
    image_data = image_file[0].data
    science.append(image_data)
    image_file.close()

for i in range(0,6):  
    scatter_file = fits.open('fullspecJupiter_center_Spectrum.0001.ec.fits')            #importing scatter frames into list scatter
    scatter_data = scatter_file[0].data
    scatter.append(scatter_data)
    scatter_file.close()

for i in range(0,6):
    science[i] = science[i]*counts[i]

for i,j in zip(range(0,6), [1.5,1.75,2,2.25,2.5,2.75]):
    scatter[i] = shift(scatter[i], j)

box_kernel = Box1DKernel(14)
for i in range(0,6):
    scatter[i] = convolve(scatter[i], box_kernel)
    
for i in range(57724,57953):
    x.append(3359.29052734 + i*0.043806734019)   
for i in range(0,6):
    y.append([science[i][j]-(multfactor[i]*(scatter[i][j])) for j in range(57724,57953)])

for i in range(0,6):
    y1.append([(science[i][j])+(10000*(i+1)) for j in range(57724,57953)])
    
for i in range(0,6):
    y1.append([(scatter[i][j])*multfactor[i]+(10000*(i+1)) for j in range(57724,57953)])

plt.clf()
plt.subplot(2, 1, 1)
plt.title('6363 $\AA$ Science Frames with Jupiter Disk Frames (Wavelength - Total Count)', fontsize = 13)
plt.ylabel('Total Pixel Count', fontsize = 13)
plt.tick_params(labelsize = 11)
plt.plot(x,y1[0], c='#EF0000')
plt.plot(x,y1[1], c='#EFB200')
plt.plot(x,y1[2], c='#78EE00')
plt.plot(x,y1[3], c='#00ED39')
plt.plot(x,y1[4], c='#00ECE9')
plt.plot(x,y1[5], c='#003DEB')
plt.plot(x,y1[6], c='#000000')
plt.plot(x,y1[7], c='#000000')
plt.plot(x,y1[8], c='#000000')
plt.plot(x,y1[9], c='#000000')
plt.plot(x,y1[10], c='#000000')
plt.plot(x,y1[11], c='#000000')



#conversion from counts to raleighs/angstrom
sensitivity = (5500000.0*1.0)/204525                #calculation of sensitivity given 5.5MR/A with 1 second exposure and 204525 counts

for i in range(0, len(y)):
    y[i] = [(x / 179.99) * sensitivity for x in y[i]]  #conversion from counts to R/A using counts --> counts/s * sensitivty --> R/A (179.99s exposure time)
    

plt.subplot(2, 1, 2)
plt.title('6363 $\AA$ Subtracted Science Frames (Wavelength - Flux)', fontsize = 13)
plt.ylabel('Flux [R/$\AA$]', fontsize = 13)
plt.xlabel('Wavelength [$\AA$]', fontsize = 13)
plt.tick_params(labelsize = 11)
plt.plot(x,y[0], c='#EF0000')
plt.plot(x,y[1], c='#EFB200')
plt.plot(x,y[2], c='#78EE00')
plt.plot(x,y[3], c='#00ED39')
plt.plot(x,y[4], c='#00ECE9')
plt.plot(x,y[5], c='#003DEB')
plt.plot(x, np.zeros(len(x)), 'k--')


