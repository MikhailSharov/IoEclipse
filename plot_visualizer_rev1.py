# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 08:33:01 2019

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

#defining any lists and variables
x = []
y = []
y1 = []
multfactor = [101536.4,109094.5,114662.1,122910.9,134650,146947]                              #minimization factors
scatter = []
science = []
counts = [100884,108497,114845,122934,134215,145975]                                        #the counts from the box on science frams for unnormalization
centercount = 204525                                                                        #the counts from the box on Jupiter Center Spectrum for unnormalization

#importing all files including scatter and science and stored in defined lists
for i in range(1,4):
    scatter_file = fits.open('fullspecJovian_Scatter.000' + str(i) + '.ec.fits')            #importing scatter frames into list scatter
    scatter_file.info()
    scatter_data = scatter_file[0].data
    print(type(scatter_data))
    print(scatter_data.shape)
    scatter.append(scatter_data)
    scatter_file.close()

for i in range(1,7):
    image_file = fits.open('fullspecIo_eclipsed.000'+str(i)+'.ec.fits')                     #importing science frames into list science
    image_file.info()
    image_data = image_file[0].data
    print(type(image_data))
    print(image_data.shape)
    science.append(image_data)
    image_file.close()
    
#unnormalizing the data and subtracting scatter
for i in range(0,6):
    science[i] = science[i]*counts[i]                                                       #unnormalization of science to counts
y.append([science[0][j]-(multfactor[0]*(scatter[1][j] + (0.0925*scatter[2][j]))/1.0925) for j in range(57733,57938)])
y.append([science[1][j]-(multfactor[1]*(scatter[1][j] + (0.131*scatter[0][j]))/1.131) for j in range(57733,57938)])
y.append([science[2][j]-(multfactor[2]*(scatter[1][j] + (0.40*scatter[0][j]))/1.40) for j in range(57733,57938)])
y.append([science[3][j]-(multfactor[3]*(scatter[1][j] + (0.79*scatter[0][j]))/1.79) for j in range(57733,57938)])
y.append([science[4][j]-(multfactor[4]*(scatter[0][j] + (0.6368*scatter[1][j]))/1.6368) for j in range(57733,57938)])
y.append([science[5][j]-(multfactor[5]*(scatter[0][j] + (0.3038*scatter[1][j]))/1.308) for j in range(57733,57938)])     #subraction of scatter from science frames
                                                                                            #only one scatter frame was used, as I found this to give the best results so far
    
#first plot with normalized data with offset
for i in range(0,6):
    y1.append([(science[i][j])+(10000*(i+1)) for j in range(57733,57938)])
for i in range(57733,57938):
    x.append(3359.29052734 + i*0.043806734019)                                              #x axis in angstroms
for i in range(0,3):
    y1.append([(scatter[i][j]*multfactor[0]) for j in range(57733,57938)])

    
    
plt.clf()
plt.subplot(2, 1, 1)
plt.plot(x,y1[0], c='#FDD321')
plt.plot(x,y1[1], c='#F4C91C')
plt.plot(x,y1[2], c='#ED9E16')
plt.plot(x,y1[3], c='#E57311')
plt.plot(x,y1[4], c='#DD490C')
plt.plot(x,y1[5], c='#EA061B')
plt.plot(x,y1[6], c='#686868')
plt.plot(x,y1[7], c='#000000')
plt.plot(x,y1[8], c='#343434')

#conversion from counts to raleighs/angstrom
sensitivity = (5500000.0*1.0)/204525                #calculation of sensitivity given 5.5MR/A with 1 second exposure and 204525 counts

for i in range(0, len(y)):
    y[i] = [(x / 179.99) * sensitivity for x in y[i]]  #conversion from counts to R/A using counts --> counts/s * sensitivty --> R/A (179.99s exposure time)

#second plot with residual in units of R/A
plt.subplot(2, 1, 2)
plt.plot(x,y[0], c='#FDD321')
plt.plot(x,y[1], c='#F4C91C')
plt.plot(x,y[2], c='#ED9E16')
plt.plot(x,y[3], c='#E57311')
plt.plot(x,y[4], c='#DD490C')
plt.plot(x,y[5], c='#EA061B')
plt.plot(x, np.zeros(len(x)), 'k--')

plt.show()
