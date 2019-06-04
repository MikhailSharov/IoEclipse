# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 20:30:37 2019

@author: msharov
"""

import numpy as np
from astropy import units as u
from astropy.nddata import CCDData
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.time import Time

multfactor = 0.9991
y = []
x = []
y1 = []
x1 = []
y2= []
x2 = []
x3 = []
y3 = []
x4 = []
y4 = []

image_file1 = fits.open('fullspecJovian_Scatter.0002.ec.fits')
image_file1.info()
image_data1 = image_file1[0].data
print(type(image_data1))
print(image_data1.shape)
image_file1.close()

image_file2 = fits.open('fullspecJovian_Scatter.0001.ec.fits')
image_file2.info()
image_data2 = image_file2[0].data
print(type(image_data2))
print(image_data2.shape)
image_file2.close()

image_file3 = fits.open('fullspecJovian_Scatter.0003.ec.fits')
image_file3.info()
image_data3 = image_file3[0].data
print(type(image_data3))
print(image_data3.shape)
image_file3.close()

"""
this group of code plots the raw data
"""
for i in range(1,5):
    image_file = fits.open('fullspecIo_eclipsed.000' + str(i) + '.ec.fits')
    image_file.info()
    image_data = image_file[0].data
    print(type(image_data))
    print(image_data.shape)
    image_file.close()
    
    
    for k in range(len(image_data)):
        image_data[k] = image_data[k]

    
    y.append([image_data[j] for j in range(57724,57953)])
   
for i in range(57724,57953):
    x.append(3359.29052734 + i*0.043806734019)

for i in range(5,7):
    image_file = fits.open('fullspecIo_eclipsed.000' + str(i) + '.ec.fits')
    image_file.info()
    image_data = image_file[0].data
    print(type(image_data))
    print(image_data.shape)
    image_file.close()
    
    
    for k in range(len(image_data)):
        image_data[k] = image_data[k]

    
    y.append([image_data[j] for j in range(57724,57953)])
"""
plt.clf()
plt.subplot(5, 1, 1)
plt.plot(x,y[0], c='#FDD321')
plt.plot(x,y[1], c='#F4C91C')
plt.plot(x,y[2], c='#ED9E16')
plt.plot(x,y[3], c='#E57311')
plt.plot(x,y[4], c='#DD490C')
plt.plot(x,y[5], c='#D51F07')
"""
"""
The follwing group of code plots the Jupiter Scatter
"""
for k in range(len(image_data2)):
    image_data2[k] = image_data2[k]

y1.append([image_data2[j]*multfactor for j in range(57724,57953)])
for i in range(57724,57953):
    x1.append(3359.29052734 + i*0.043806734019)

for k in range(len(image_data1)):
    image_data1[k] = image_data1[k]

y2.append([image_data1[j]*multfactor for j in range(57724,57953)])
for i in range(57724,57953):
    x2.append(3359.29052734 + i*0.043806734019)

counts = [100884,108497,114845,122934,134215,145975]
"""
plt.subplot(5, 1, 2)
plt.plot(x1,y1[0], c='black')
plt.plot(x2,y2[0], c='grey')
"""
"""
The following group is plotting the two on the same plot
"""
plt.subplot(2, 1, 1)
plt.plot(x,y[0], c='#FDD321')
plt.plot(x,y[1], c='#F4C91C')
plt.plot(x,y[2], c='#ED9E16')
plt.plot(x,y[3], c='#E57311')
plt.plot(x,y[4], c='#DD490C')
plt.plot(x,y[5], c='#EA061B')
plt.plot(x1,y1[0], c='black')
plt.plot(x2,y2[0], c='grey')

"""
The following group is plotting the subtraction of the scatter from the science
"""
for i in range(1,5):
    image_file = fits.open('fullspecIo_eclipsed.000' + str(i) + '.ec.fits')
    image_file.info()
    image_data = image_file[0].data
    print(type(image_data))
    print(image_data.shape)
    image_file.close()
    
    
    for k in range(len(image_data)):
        image_data[k] = image_data[k]-(image_data1[k]*multfactor)

    
    counts = [100884,108497,114845,122934,134215,145975]
    
    y3.append([image_data[j]*counts[i-1] for j in range(57724,57953)])
   
for i in range(57724,57953):
    x3.append(3359.29052734 + i*0.043806734019)

for i in range(5,7):
    image_file = fits.open('fullspecIo_eclipsed.000' + str(i) + '.ec.fits')
    image_file.info()
    image_data = image_file[0].data
    print(type(image_data))
    print(image_data.shape)
    image_file.close()    
    
    for k in range(len(image_data)):
        image_data[k] = image_data[k]-((image_data1[k]+image_data2[k])*multfactor/2)

    
    counts = [100884,108497,114845,122934,134215,145975]
    
    y3.append([image_data[j]*counts[i-1] for j in range(57724,57953)])

    
plt.subplot(2, 1, 2)
plt.plot(x3,y3[0], c='#FDD321')
plt.plot(x3,y3[1], c='#F4C91C')
plt.plot(x3,y3[2], c='#ED9E16')
plt.plot(x3,y3[3], c='#E57311')
plt.plot(x3,y3[4], c='#DD490C')
plt.plot(x3,y3[5], c='#EA061B')
plt.plot(x3,np.zeros(len(x3)), 'k--')

for i in range(1,5):
    image_file = fits.open('fullspecIo_eclipsed.000' + str(i) + '.ec.fits')
    image_file.info()
    image_data = image_file[0].data
    print(type(image_data))
    print(image_data.shape)
    image_file.close()
    
    
    for k in range(len(image_data)):
        image_data[k] = image_data[k]-(image_data1[k]*multfactor)
        if image_data[k] < 0:
            image_data[k] = 0

    
    counts = [100884,108497,114845,122934,134215,145975]
    
    y4.append([image_data[j]*counts[i-1] for j in range(57724,57953)])
   
for i in range(57724,57953):
    x4.append(3359.29052734 + i*0.043806734019)

for i in range(5,7):
    image_file = fits.open('fullspecIo_eclipsed.000' + str(i) + '.ec.fits')
    image_file.info()
    image_data = image_file[0].data
    print(type(image_data))
    print(image_data.shape)
    image_file.close()    
    
    for k in range(len(image_data)):
        image_data[k] = image_data[k]-((image_data1[k]+image_data2[k])*multfactor/2)
        if image_data[k] < 0:
            image_data[k] = 0

    
    counts = [100884,108497,114845,122934,134215,145975]
    
    y4.append([image_data[j]*counts[i-1] for j in range(57724,57953)])
""" 
plt.subplot(5, 1, 5)
plt.plot(x4,y4[0], c='red')
plt.plot(x4,y4[1], c='gold')
plt.plot(x4,y4[2], c='green')
plt.plot(x4,y4[3], c='blue')
plt.plot(x4,y4[4], c='purple')
plt.plot(x4,y4[5], c='magenta')
"""
plt.show()