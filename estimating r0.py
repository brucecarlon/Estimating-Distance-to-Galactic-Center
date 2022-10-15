# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 17:53:18 2020

@author: user
"""
import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.mlab as mlab
import numpy.random as random
from  astropy.stats import jackknife_resampling
from astropy.stats import jackknife_stats
from scipy.optimize import curve_fit #uses LM algorithm through least sq
from statistics import mean
import math
f = open('LMC_t2cephs.dat.txt','r')
header = f.readline()

#lmc_v = []
#lmc_p = []
#lmc_logp =[] 
#udata = []

i = 0
n = 187
lmc_v = np.zeros(n)
lmc_i = np.zeros(n)
lmc_p = np.zeros(n)
lmc_logp = np.zeros(n)
WI = np.zeros(n)
v_0 = np.zeros(n)

#tdata,ydata = [],[]
for line in f : 
    line = line.strip()
    columns = line.split()
    lmc_i[i] = float(columns[1])
    lmc_v[i] = float(columns[2])
    lmc_p[i] = float(columns[3])
    lmc_logp[i] = float(columns[4])
    WI[i] = float(columns[5])
    v_0[i] = float(columns[6])
    i = i + 1
   
vi_0 = np.zeros(n)
m = -2.62
c = 17.33
R = 2.55
EBI = 0.09
Av = 0.005
for i in range(n):
    vi_0[i] = 1/R*(v_0[i] - m*lmc_logp[i] - c)
    
print('mean =',np.mean(vi_0))
####    LOGP vs MAGV (LMC)   
plt.scatter(lmc_logp,lmc_v)
plt.xlabel('log(P)')
plt.ylabel('V mag')
plt.savefig('pl')
plt.show()


                                ########################################
ydata = WI
xdata = lmc_logp                               
N = len(xdata)
# calculate the sums needed for the least squares fit

sum_xy = 0.0
sum_x = 0.0
sum_y = 0.0
sum_xx = 0.0
sum_dd = 0.0

for j in range(N):

    sum_xy += xdata[j]*ydata[j]
    sum_x += xdata[j]
    sum_y += ydata[j]
    sum_xx += xdata[j]*xdata[j]


# do the linear least squares fit
# (compare these to equations on page 110 of the measurement manual)

mfit = (N*sum_xy - sum_x * sum_y) / (N*sum_xx - sum_x * sum_x)
cfit = (sum_xx*sum_y - sum_xy*sum_x) / (N*sum_xx - sum_x * sum_x)

for j in range(N):
    d = ydata[j] - (mfit * xdata[j] + cfit)
    sum_dd += d*d;

umfit = math.sqrt((sum_dd)/(N*sum_xx - sum_x * sum_x)*N/(N-2))
ucfit = math.sqrt((sum_dd * sum_xx)/(N*(N*sum_xx - sum_x * sum_x))*N/(N-2))

# print the results

print("m = %6.2f +/- %5.2f"%(mfit,umfit))
print("c = %6.2f +/- %5.2f"%(cfit,ucfit))


## plot the model prediction with the best-fit parameters 
## (you will need matplotlib and numpy installed to do this)
#
import numpy as np
import matplotlib.pyplot as plt
#
xdata_array = np.array(xdata)
ydata_array = np.array(ydata)

plt.scatter(lmc_logp,WI, label = 'data' )
plt.xlabel('WI')
plt.ylabel('logp')

#plt.plot(xdata_array, ydata_array, 'bs', label = "Data")
#
yfit = mfit*xdata_array + cfit
plt.plot(xdata_array, yfit, '-r', label = "Line of Best Fit")
#
#plt.title("The result of the model fit to the linear data")
plt.ylabel("WI")
plt.xlabel("log(P)")
plt.legend()
plt.savefig('plc')
plt.show()
                                
                                
                                #######################################
    
