#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 18:58:50 2024

@author: lorenzopizzuti

Module containing all the cluster info, cosmological parameters

NOTE: For the gas mass profile of  MACS 1720, the double beta model is not working (odd values of parameters). 
Use the single beta-model instead! 

"""

import numpy as np
import scipy.integrate as ints


#Cosmological parameters:
G=4.302e-9
H0=70
Om=0.3
Ol=0.7



# Name of the cluster according to Lensing data
nameclust=np.array(['a209',
       'a383',
      'a611', 
       'a2261',
       'macsj0329', 
       'macsj0416',
       'macsj0429', 
       'macsj0647',
       'macsj0717', 
       'macsj0744', 
       'macsj1115', 
       'macsj1149', 
       'macsj1206', 
      'macsj1720', 
       'macsj1931',
       'ms2137', 
       'rxj1347',
       'rxj1532', 
       'rxj2129',
       'rxj2248'],dtype=str)


# r200 values of the cluster (given by lensing)
r200=np.array([2.323925,
 1.7268267857142858,
 2.120977142857143,
 2.4829532142857147,
 1.6335985714285715,
 1.8391732142857145,
 1.7119378571428574,
 1.8240728571428573,
 2.416177142857143,
 1.916089285714286,
 2.1618860714285715,
 2.3142478571428575,
 2.112519285714286,
 1.9882757142857146,
 2.166578571428572,
 1.9948735714285717,
 2.5992032142857147,
 1.4722378571428574,
 1.5643414285714285,
 2.155100714285714],dtype=float)

# Average redshift of the cluster
za=np.array([0.206,0.187,0.288,0.224,0.450,
             0.396,0.399,0.584,0.548,0.686,0.352,0.544,0.440,
    0.391,0.352,0.313,0.451,0.345,0.234,0.348],dtype=float)


#Cluster gas data according to Salzano and checked on the Chandra Data Archive
#NOTE: For MACS 1720 the double beta model is not working (odd values of parameters). 
#Use the single beta-model instead! 

#******************************************************************************
# 0: name of cluster 
# 1: number of cluster data points used for the best fit. I include this for comparison reasons. 
#### Some clusters have very few data, so it might be useful to use this as a test.
# 2: ne0 - electron density 
# 3: r0 - characteristic radius in beta model
# 4: alpha - exponent in beta model
# 5: re0  - characteristic radius  in beta model
# 6: beta0 - exponent in beta model
# 7: ne1  - similar to ne0, but for double beta model
# 8: re1 - similar to re0, but for double beta model
# 9: beta1 - similar to beta0, but for double beta model
# 10: relaxation state
# Note: For the clusters which are fit with a single beta model the entries 7,8,9 are = 0.


cluster_data = [['A209', 4.0, 0.0733712, 47585.8, 0.507418, 299.051, 0.474847, 0, 0, 0, 1],
['A611', 16.0, 0.0677311, 36559.1, 0.633827, 152.511, 0.408396, 0, 0, 0, 0],
['A2261', 14.0, 8.58559, 100478., 3.18934e-7, 28.314, 0.429423, 1.7609, 179.855, 1.,0],
#['MACS0329', 24, 44.39367934,   13.33108055, -1.20122628, 19.1043168, 1.2308636, 
# 1.78694419, -209.24418617, 0.72303934, 1],
['MACS0329', 24,2.44484568e+01, 1.64167019e-01, 5.85505620e-08, 4.00077116e+01,
        1.21601137e+00, 2.38352733e+00, 1.78760282e+02, 6.98601473e-01, 1],
['MACS0416', 3.0, 1.31078, 152.779, 0.0102307, 414.69, 0.954723, 0, 0, 0, 1],
['MACS0429', 5.0, 1.17371971e+01, 1.40725014e+00, 2.91809574e-07, 6.00658535e+02,
       3.16261348e+02, 1.23365937e+01, 4.27381188e+01, 5.36021032e-01, 0],
['MACS0647', 3.0, 9.40547, 18.4831, 0.670707, 497.919, 0.990862, 0, 0, 0, 1],
['MACS0744', 5.0, 3.47974, 103.019, 0.947664, 289.22, 0.355337, 0, 0, 0, 0],
['MACS1115', 10.0, 1.90507435e+01, 2.87921965e+03, 4.74349188e-30, 4.00717121e+03,
       8.84343950e+03, 6.93694732e+00, 1.08621075e+02, 6.68815776e-01, 0],
['MACS1149', 3.0,  0.124954, 59731.0, 0.401771, 564.733, 0.571188, 0, 0, 0, 1],
['MACS1206', 10.0,  0.083374, 35192.4,  0.646376, 251.218, 0.448791, 0, 0, 0, 0],
['MACS1720', 6.0, 7.86793, 35.3772, 1.0, 288.243, 0.35858, 420260., 25.809, 0.856979, 0],
['MS2137', 24, 7.19048317e+00, 6.66017747e+01, 6.61309517e-01, 4.36666437e+03,
       4.58800822e+03, 6.25193014e+00, 1.03932983e+02, 7.81272900e-01, 1],
['RXJ1347' , 24.0, 3.32492587e+00, 4.20004428e+01, 7.72379499e-20, 3.16571584e+02,
       2.66378335e+00, 5.02110825e+01, 1.50379603e+01, 4.65907658e-01, 1],
['RXJ2129', 24.0, 4.33718813e+00, 3.70246256e+01, 9.83175299e-01, 2.45123867e+03,
       2.07090002e+03, 5.95066775e+00, 7.91614638e+01, 6.04192549e-01, 0],
['RXJ2248', 18.0, 0.179542, 47882.7, 0.517691, 298.564, 0.642956, 0, 0, 0, 1]]

cluster_data = np.array(cluster_data)

# Defines the beta model function. 
# Selects whether it is standard beta model, or double beta model based on the best fit parameters.
def rho_fit(x, ne0, r0, alpha, re0, beta0, ne1, re1, beta1):
    if ne1 == 0:
        return  ((ne0*10.0**5.0)*(x/r0)**(-alpha))*(1.0 + (x/re0)**2.0)**(-1.5*beta0) 
    else:
        return  ((ne0*10.0**5.0)*(x/r0)**(-alpha))*(1.0 + (x/re0)**2.0)**(-1.5*beta0) + (ne1*10**5)*(1 + (x/re1)**2)**(-1.5*beta1)



def rhoint(x,ne0, r0, alpha, re0, beta0, ne1, re1, beta1):
    return rho_fit(x,ne0, r0, alpha, re0, beta0, ne1, re1, beta1)*x**2


def mgas(x,ne0, r0, alpha, re0, beta0, ne1, re1, beta1):
    """
    
    Parameters
    ----------
    x : FLOAT (1D array)
        radius in kpc.
    ne0 : TYPE
        DESCRIPTION.
    r0 : TYPE
        DESCRIPTION.
    alpha : TYPE
        DESCRIPTION.
    re0 : TYPE
        DESCRIPTION.
    beta0 : TYPE
        DESCRIPTION.
    ne1 : TYPE
        DESCRIPTION.
    re1 : TYPE
        DESCRIPTION.
    beta1 : TYPE
        DESCRIPTION.

    Returns
    -------
    FLOAT (1D array)
        The value of the mass profile ad x.

    """
    
    return 4*np.pi*ints.quad(rhoint,0,x, args=(ne0, r0, alpha, re0, beta0, \
                                               ne1, re1, beta1), \
                        epsabs=1.49e-4, epsrel=1.49e-8)[0]
mgas=np.vectorize(mgas)



#Cluster BCG masses ===========================================================

cluster_BCG_mass = [['MACS0329', 3.69,0.21,0.79,0.06],
['MACS1115', 2.19,0.17,0.51,0.04],
['MS2137', 8.38,0.42,0,0],
['RXJ1347',3.68,0.85,1.30,0.10],
['RXJ2129',2.21,0.32,1.16,0.09],
['A209', 2.0,0.16,0,0.0],
['A611', 3.45,0.42,0.31,0.03],
['A2261', 1.74,0.18,4.14,0.33],
['MACS0416', 2.89,0.99,1.50,0.12],
['MACS0429', 4.68,0.23,1.61,0.13],
['MACS0647', 5.13,1.79,10.43,0.83],
['MACS0744', 8.20,0.43,0.65,0.05],
['MACS1149', 3.05,0.58,1.10,0.09],
['MACS1206', 4.93,1.39,1.10,0.09],
['MACS1720',3.83,0.94,0.78,0.06],
['RXJ2248',5.33,0.61,1.31,0.10]]


#Cluster temperature: Name, T, error(68%) =====================================
cluster_Temperature = np.array([['A209', 7.3, 0.5], 
['A611', 7.90, 0.35],
['A2261', 7.6, 0.3],
['MACS0329', 8.0, 0.5],
['MACS0416', 7.5, 0.8],
['MACS0429', 6.00, 0.44],
['MACS0647', 13.3, 1.8],
['MACS0744', 8.9, 0.8],
['MACS1115', 8.0, 0.4], 
['MACS1149', 8.7, 0.9],
['MACS1206', 10.8, 0.6],
['MACS1720', 6.6, 0.4],
['MS2137', 5.9, 0.3],
['RXJ1347', 15.5, 0.6],  
['RXJ2129', 5.8, 0.4],
['RXJ2248', 12.4, 0.6]], dtype = None)





