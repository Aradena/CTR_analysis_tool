import os

import pandas as pd
import numpy as np
import scipy
import scipy.optimize as opt

# dataPath='data/10102022_Liroc_PZC8_Vth460_noProbePA_allMasked_ch58_50ohm_noExtPZ_EPIC2x2x3_PbF2_black_FBK_NUVHDRH_UHDDA_vs_HF_TAC2x2x3_LYSOCeCa_4FP_BRCM_Nr3/'
# filename='staircase_board6_liroc_HV45_PZC8.csv'
from matplotlib import pyplot as plt
import matplotlib
from scipy.interpolate import splrep, splev
from scipy.optimize import curve_fit
from math import *
matplotlib.use('Qt4Agg')
font = {'size'   : 22}
matplotlib.rc('font', **font)


def modelGauss(x, mu, sigma, amp):
    return amp * (1 / (sigma * (np.sqrt(2 * np.pi)))) * (np.exp((-1.0 / 2.0) * (((x - mu) / sigma) ** 2)))


def myerf(x, A, mu, sigma):
    return A * (1 + scipy.special.erf(-(x - mu) / (sqrt(2) * sigma)))

# sigmoid function
def f(x, c, d):
    # return ((-1 / (1 + np.exp(-c * (x - d)))) + 1)
    return (1 / (1 + np.exp(-c * (x - d))))

def sigmoid(x, L, x0, k, b):
    y = L / (1 + np.exp(-k * (x - x0))) + b
    return (y)

def fitStairCase(dataPath, filename, channel, method, leftEdge, rightEdge):

    df=pd.read_csv(dataPath+filename, sep=' ', usecols=[0,channel+1], names=['threshold','hits'])
    # print(df)
    if not os.path.exists(dataPath+"/"+method):
        os.mkdir(dataPath+"/"+method)


    # # df['derivative']=np.diff(df['hits'],prepend=df['hits'][0])
    # fig2,ax2 = plt.subplots(1,1,figsize=(20,10))
    # ax2.set_yscale('log')
    # plt.scatter(df['threshold'], df['hits'], marker='+', label='raw efficiency')
    # plt.legend()
    # plt.xlabel("Global threshold [dacu]")
    # plt.ylabel("Raw trigger efficiency (log scale)")
    # plt.show()


    mask=(df['threshold']>leftEdge) & (df['threshold']<rightEdge)
    x=df.loc[mask,'threshold'].to_numpy()
    y=np.log(df.loc[mask,'hits'].to_numpy())
    ymin = min(y)
    # print(ymin)
    y = y - ymin
    ymax=max(y)
    # print(ymax)
    y=y/ymax

    # using diff method
    # maskDeriv = (df['threshold'] > leftEdge) & (df['threshold'] < rightEdge+1)
    # yderiv = np.diff(np.log(df.loc[maskDeriv, 'hits'].to_numpy()))
    # derivee centrale + moyennage sur bords
    yderiv = np.gradient(y, edge_order=2)
    p0 = [x[np.argmax(yderiv)], 1, max(y)]


    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    xfit = np.linspace(min(x), max(x), 100)

    if method == "erf":
        p0 = [0.5, np.median(x), -0.5]
        popt, pcov = curve_fit(myerf, x, y, p0, method='dogbox')
        yfit = myerf(xfit, *popt)
        noise = abs(popt[2])
        threshold=popt[1]
        amplitude=popt[0]
        plt.plot(xfit, yfit, label='Erf fit\n $\sigma$=%.2f dacu\n $\mu$=%.1f dacu\n Ampl=%.1f' % (noise,threshold,amplitude), linestyle='--')

    elif method == "sigmoid":
        p0 = [max(y), np.median(x),1,min(y)] # this is an mandatory initial guess
        # b0 = ([0.1,415],[20,455]) # bounds
        popt, pcov  = opt.curve_fit(f,x, y, method='dogbox') #p0=p0,bounds=b0)
        noise = abs(popt[0])
        yfit = f(xfit, *popt)
        plt.plot(xfit, yfit, label='Sigmoid fit with sigma= %.2f dacu' % noise, marker='+')

    elif method == "spline":
        spl = splrep(x, y, k=1, s=0)
        yfit = splev(xfit, spl)
        plt.plot(xfit, yfit, label='B-spline fit', marker='+')

        knots, BspliceCoef, _ = spl
        # indexes= np.where((knots>431) & (knots<437))[0]
        thres50=x[np.argmax(yderiv)]
        indexes= np.where((knots>thres50-5) & (knots<thres50+5))[0]
        for slopeIndex in indexes:
            # print(BspliceCoef[slopeIndex])
            ax.text(knots[slopeIndex], y[slopeIndex], "%.2f"%BspliceCoef[slopeIndex])

    elif method == "derivative":

        # print(p0)
        popt, cov = opt.curve_fit(modelGauss, x, yderiv, p0=p0)
        yfit = modelGauss(xfit, *popt)
        plt.scatter(x, yderiv, marker='+', label='derivative')
        noise, threshold, amplitude = popt[1],popt[0],popt[2]
        plt.plot(xfit, yfit, linestyle='--', label='Gaussian fit\n $\sigma$=%.2f dacu\n $\mu$=%.1f dacu\n Ampl=%.1f' % (noise,threshold,amplitude))


    #  np.log(df.loc[mask,'hits'].to_numpy()

    # print("noise ", noise)
    # print("threshold ", popt[1])

    # fitParam, cov = opt.curve_fit(modelGauss, x, y,p0=[edges[np.argmax(count)], .01, max(count)])
    #
    # meanHisto = fitParam[0]
    # sigmaHisto = fitParam[1]
    # ampHisto = fitParam[2]


    plt.plot(x,y, marker='+', label='raw data')
    plt.legend()
    plt.xlabel("Global threshold [dacu]")
    plt.ylabel("Normalized trigger efficiency ")
    plt.tight_layout()
    plt.savefig(dataPath+"/"+method+"/"+filename+".png")
    plt.close()


    # gaussian fit noise
    mask=(df['threshold']>leftEdgeGaussian) & (df['threshold']<rightEdgeGaussian)
    x2=df.loc[mask,'threshold'].to_numpy()
    y2=np.log(df.loc[mask,'hits'].to_numpy())

    p0=[x2[np.argmax(y2)], 1, max(y2)]
    popt2, pcov2 = curve_fit(modelGauss, x2, y2,p0) #, method='dogbox')
    xfit2 = np.arange(min(x2), max(x2), .1)
    yfit2 = np.exp(modelGauss(xfit2, *popt2))
    sigma=popt2[1]

    fig2,ax2 = plt.subplots(1,1,figsize=(20,10))
    ax2.set_yscale('log')
    plt.scatter(df['threshold'], df['hits'], marker='+')
    plt.plot(xfit,np.exp(yfit*ymax+ymin),label='fit', linestyle='--', color='tab:red')
    plt.plot(xfit2,yfit2,label='fit gaussian sigma=%.2f dacu'%sigma, linestyle='--', color='tab:orange')
    plt.legend()
    plt.xlabel("Global threshold [dacu]")
    plt.ylabel("Raw trigger efficiency (log scale)")
    plt.tight_layout()
    plt.savefig(dataPath+filename+"_full.png")
    plt.xlim([leftEdge-10,rightEdgeGaussian+10])
    plt.ylim([0.5e4,1e9])
    plt.savefig(dataPath+"/"+method+"/"+filename+"_zoom.png")
    # plt.show()

if __name__ == '__main__':
    dataPath = 'data/LIROC_UHD-DE/'
    channel = 58
    # filename = 'Staircase_pzc8_HV45.csv'
    filename='Staircase_pzc8_HV50.csv'

    leftEdge, rightEdge = 415,450 # 45V
    # 50V
    leftEdge, rightEdge = 390,445
    leftEdgeGaussian, rightEdgeGaussian = 475,490

    fitStairCase(dataPath, filename, channel, "spline",leftEdge, rightEdge)
    fitStairCase(dataPath, filename, channel, "sigmoid",leftEdge, rightEdge)
    fitStairCase(dataPath, filename, channel, "erf",leftEdge, rightEdge)
    fitStairCase(dataPath, filename, channel, "derivative", leftEdge, rightEdge)
