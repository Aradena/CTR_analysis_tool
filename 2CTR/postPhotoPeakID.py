import scipy.optimize as opt
import os
import matplotlib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
matplotlib.use('Qt5Agg')
font = {'size'   : 22}
matplotlib.rc('font', **font)
from scipy.stats import norm
def modelGauss(x, amp, mu, sigma):
    return amp * (1 / (sigma * (np.sqrt(2 * np.pi)))) * (np.exp((-1.0 / 2.0) * (((x - mu) / sigma) ** 2)))


dataPath= '14102022_Radio_vs_Radio_BiasScan_allMasked_ch62_10ohm_comp0_gain2_PMI1X025_vs_PMI1X026_LYSOCeCa_BRCM_NUVMT'
dataPath='15102022_Radio_vs_Radio_allMasked_ch62_10ohm_PMI3X001_vs_PMI3X002_BGO_3x3x20_BRCM_NUVMT'

maskFiles = [filename for filename in os.listdir('outputBokeh/'+dataPath+'/') if filename.startswith("F3") and filename.endswith(".txt.csv")]
# maskFiles= maskFiles[-2:-1]
# print(maskFiles)

# maskFiles= maskFiles[-2:0]
# maskFiles= maskFiles[0:1]
number=5
maskFiles= maskFiles[number:number+1]
for maskFile in maskFiles:
    print(maskFile)
    df_mask = pd.read_csv('outputBokeh/'+dataPath+'/'+maskFile, usecols=['ID', 'deltaT', 'cluster'])

    maskFile=maskFile.replace(".csv", "")
    riseTimeA_file=maskFile.replace("F3", "F4")
    riseTimeB_file=maskFile.replace("F3", "F5")

    csvSeparator = '\t'  # ',' #
    riseTimeA_df = pd.read_csv('data/' + dataPath + '/' + riseTimeA_file, skiprowupdats=4, encoding="ISO-8859-1",
                               sep=csvSeparator, na_values="-nan(ind)")
    riseTimeA_df = riseTimeA_df.dropna().rename(columns={'Time': 'ID', 'Ampl': 'riseTimeA'})
    riseTimeA_df['riseTimeA'] *= 10 ** 9

    riseTimeB_df = pd.read_csv('data/' + dataPath + '/' + riseTimeB_file, skiprows=4, encoding="ISO-8859-1",
                               sep=csvSeparator, na_values="-nan(ind)")
    riseTimeB_df = riseTimeB_df.dropna().rename(columns={'Time': 'ID', 'Ampl': 'riseTimeB'})
    riseTimeB_df['riseTimeB'] *= 10 ** 9

    df = pd.merge(riseTimeA_df, riseTimeB_df, how='inner', on="ID")
    df = pd.merge(df, df_mask, how='inner', on="ID")

    #
    bins = np.linspace(df['deltaT'].min(), df['deltaT'].max(), 1000)
    hist, bin_edges = np.histogram(df['deltaT'], bins=bins)
    hist = np.append(hist, 0)
    xmax = bin_edges[np.argmax(hist)]
    ymax = max(hist)
    fitParam, cov = opt.curve_fit(modelGauss, bin_edges, hist, p0=[ymax, xmax, 300])
    mean = fitParam[1]
    sigma = fitParam[2]
    xmin = min(bin_edges)
    xmax = max(bin_edges)
    xfit = np.arange(xmin, xmax, 1e-3)
    yfit = modelGauss(xfit, *fitParam)

    fig,ax = plt.subplots(2,2,figsize=(20,10)) #, sharex=True, sharey=True)

    ax[1, 0].hist(df['deltaT'], bins=bins, density=False, label='deltaT')
    ax[1, 0].plot(xfit, yfit, 'r')
    ax[1, 0].legend()

    # select for photopeak
    zoom_sigma=2
    # df=df[(df['cluster']==0)]
    df=df[(df['cluster']==0) & (df['deltaT']<(mean+zoom_sigma*sigma)) & (df['deltaT']>(mean-zoom_sigma*sigma))]
    print(len(df))


    ax[0,0].scatter(df['deltaT'],df['riseTimeA'], marker='+', label='A')
    ax[0,0].scatter(df['deltaT'],df['riseTimeB'], marker='o', label='B')
    ax[0,0].legend()
    ax[0,0].set_xlabel("deltaT [ns]")
    ax[0,0].set_ylabel("Rise Time [ns] ")
    # ax[0, 0].set_xlim([mean-zoom_sigma*sigma, mean+zoom_sigma*sigma])

    bins = np.linspace(df['riseTimeA'].min(), df['riseTimeA'].max(), 1000)
    ax[0, 1].hist(df['riseTimeA'], bins=bins, density=False, label='histo rise time A')
    # ax[0, 1].set_xlim([mean-zoom_sigma*sigma, mean+zoom_sigma*sigma])
    ax[0, 1].legend()
    # ax[1, 0].text(0.1,0.9,"Zoom on mu+-"+str(zoom_sigma)+"*sigma")

    # ax[0, 1].xlabel("deltaT [ns]")
    # ax[0, 1].ylabel("Rise Time [ns] ")

    bins = np.linspace(df['riseTimeB'].min(), df['riseTimeB'].max(), 1000)
    ax[1, 1].hist(df['riseTimeB'], bins=bins, density=False, color='tab:orange', label='histo rise time B')
    ax[1, 1].set_title("Zoom on deltaT mu+-" + str(zoom_sigma) + "*sigma")
    # ax[1, 1].set_xlim([mean-zoom_sigma*sigma, mean+zoom_sigma*sigma])

    ax[1, 1].legend()
    fig.suptitle(maskFile)


    fig.tight_layout()
    # plt.savefig('outputBokeh/'+dataPath+"/"+maskFile+"riseTime.png")
    plt.show()
    print("done")