import pandas as pd
import numpy as np
import scipy
import scipy.optimize as opt

# dataPath='data/10102022_Liroc_PZC8_Vth460_noProbePA_allMasked_ch58_50ohm_noExtPZ_EPIC2x2x3_PbF2_black_FBK_NUVHDRH_UHDDA_vs_HF_TAC2x2x3_LYSOCeCa_4FP_BRCM_Nr3/'
# filename='staircase_board6_liroc_HV45_PZC8.csv'
from matplotlib import pyplot as plt
import matplotlib
from scipy.optimize import curve_fit
from math import *
matplotlib.use('Qt4Agg')
font = {'size'   : 22}
matplotlib.rc('font', **font)

dataPath='data/LIROC_UHD-DE/'
filename='Staircase_pzc8_HV45.csv'
channel=58
df=pd.read_csv(dataPath+filename, sep=' ', usecols=[0,channel+1], names=['threshold','hits'])
# print(df)

df['derivative']=np.diff(df['hits'],prepend=df['hits'][0])


def modelGauss(x, mu, sigma, amp):
    return amp * (1 / (sigma * (np.sqrt(2 * np.pi)))) * (np.exp((-1.0 / 2.0) * (((x - mu) / sigma) ** 2)))

def myerf(x,A,mu,sigma):
   return A*(1+ scipy.special.erf( -(x - mu)/(sqrt(2)*sigma)))

mask=(df['threshold']>420) & (df['threshold']<450)
x=df.loc[mask,'threshold'].to_numpy()
y=np.log(df.loc[mask,'hits'].to_numpy())
ymin = min(y)
print(ymin)
y = y - ymin
ymax=max(y)
print(ymax)
y=y/ymax


# print(y)# p0=[0.5, 435, -0.5]
# popt, pcov = curve_fit(myerf, x, y) #, p0=p0)
# print("noise ", abs(popt[2]))
# print("threshold ", popt[1])
# yfit=myerf(x,*popt)
# yfit=myerf(x,popt[0], popt[1], popt[2])


# sigmoid function
def f(x, c, d):
    # return ((-1 / (1 + np.exp(-c * (x - d)))) + 1)
    return (1 / (1 + np.exp(-c * (x - d))))

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)

p0 = [max(y), np.median(x),1,min(y)] # this is an mandatory initial guess

popt, pcov = curve_fit(sigmoid, x,y,p0, method='dogbox')

#  np.log(df.loc[mask,'hits'].to_numpy()
# # Fit sigmoid
# p0 = [5,435] # this is an mandatory initial guess
# b0 = ([0.1,415],[20,455]) # bounds
#
# popt, pcov  = opt.curve_fit(f,x, y, method='dogbox') #p0=p0,bounds=b0)

noise=abs(popt[0])
print("noise ", noise)
print("threshold ", popt[1])
yfit=sigmoid(x,*popt)

# fitParam, cov = opt.curve_fit(modelGauss, x, y,p0=[edges[np.argmax(count)], .01, max(count)])
#
# meanHisto = fitParam[0]
# sigmaHisto = fitParam[1]
# ampHisto = fitParam[2]

fig,ax = plt.subplots(1,1,figsize=(10,10))
plt.plot(x,y, marker='+', label='raw data')
plt.plot(x,yfit,label='fit noise= %.2f dacu'%noise, marker='+')
plt.legend()
plt.xlabel("Global threshold [dacu]")
plt.ylabel("Normalized trigger efficiency ")
plt.savefig(dataPath+filename+".png")
plt.close()


# gaussian fit noise
mask=(df['threshold']>477) & (df['threshold']<490)
x2=df.loc[mask,'threshold'].to_numpy()
y2=np.log(df.loc[mask,'hits'].to_numpy())
def modelGauss(x, mu, sigma, amp):
    return amp * (1 / (sigma * (np.sqrt(2 * np.pi)))) * (np.exp((-1.0 / 2.0) * (((x - mu) / sigma) ** 2)))

p0=[x2[np.argmax(y2)], 1, max(y2)]
popt2, pcov2 = curve_fit(modelGauss, x2, y2,p0) #, method='dogbox')
xfit2 = np.arange(min(x2), max(x2), .1)
yfit2 = np.exp(modelGauss(xfit2, *popt2))
sigma=popt2[1]

fig2,ax2 = plt.subplots(1,1,figsize=(20,10))
ax2.set_yscale('log')
plt.scatter(df['threshold'], df['hits'], marker='+')
plt.plot(x,np.exp(yfit*ymax+ymin),label='fit noise= %.2f dacu'%noise, linestyle='--', color='tab:red')
plt.plot(xfit2,yfit2,label='fit gaussian sigma=%.2f dacu'%sigma, linestyle='--', color='tab:orange')
plt.legend()
plt.xlabel("Global threshold [dacu]")
plt.ylabel("Raw trigger efficiency (log scale)")
plt.savefig(dataPath+filename+"_full.png")
plt.xlim([415,495])
plt.ylim([0.5e4,1e9])
plt.savefig(dataPath+filename+"_zoom.png")
plt.show()
