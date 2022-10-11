import datetime
import os
import re
import webbrowser

import scipy.optimize as opt
from scipy.special import erf

import numpy as np
import pandas as pd
from bokeh.io import export_png, export_svg

from bokeh.layouts import gridplot, column, row
from bokeh.models import BoxSelectTool, LassoSelectTool, ColumnDataSource, Slider, TextInput, Paragraph, Button, Select, \
    Text, Div, Range1d, Legend, WheelZoomTool, CustomJS, Toggle
from bokeh.plotting import curdoc, figure

from matplotlib import pyplot as plt
# debug in pycharm : https://open-biotech.github.io/bio-rtd-docs/getting_started.html

TOOLS="pan,wheel_zoom,box_zoom,box_select,lasso_select,reset"

data = pd.DataFrame(columns=['ID','area_trigger_asic', 'amplitude_ref', 'areaPALG_asic', 'charge_asic', 'deltaT', 'cluster'])
# charge_asic column corresponds to selected column between area_trigger_asic and areaPALG_asic
src = ColumnDataSource(data)
res_dict = {}

squaredFigSize = 600

lirocList = ['HV_SiPM_asic', 'threshold_asic', 'PZC_asic', 'HV_SiPM_ref', 'threshold_ref', 'detectorASICside',
             'detectorREFside', 'noProbePA', 'asicChannel', 'extPZC_circuit', 'externalResAtChannelInput', 'allMasked']
sipmCaracList = ['sipmModel', 'HV_SiPM_HF', 'thres_HF', 'HV_SiPM_ref']


def findNumberIfAny(patternToFind, string, numberBefore=False):
    if re.search(r''.join(patternToFind), string) is None:
        return None
    elif re.search(r''.join(patternToFind)+'lab', string) is not None:
        tmp = str(re.findall(r'(?<=_)[0-9]+(?=.txt)', string)[0])
        print("Automatic labview sweep detected with param: ",tmp)
        return tmp+' sweep'
    else:
        if numberBefore:
            tmp = re.findall(r'[0-9]+(?=' + patternToFind + ')', string)
            if len(tmp)==0:
                return None
            else:
                return int(tmp[0])
        else:
            return int(re.findall(r'(?<='+patternToFind+')[+-]?[0-9]+', string)[0])

def plotMainFig(update):
    p = figure(tools=TOOLS, width=squaredFigSize, height=squaredFigSize, min_border=0, #10, min_border_left=50,
               toolbar_location="above", x_axis_location=None, y_axis_location=None,
               title="Charge ASIC versus charge REF")
    p.background_fill_color = "#fafafa"
    p.select(BoxSelectTool).select_every_mousemove = False
    p.select(LassoSelectTool).select_every_mousemove = False

    r = p.scatter('amplitude_ref', 'charge_asic', source=src, size=3, color="#3A5785", alpha=0.6)
    r.data_source.selected.on_change('indices', update)
    return p

def plotSecondFig(fig, update):
    p = figure(tools=TOOLS, width=squaredFigSize, height=squaredFigSize, min_border=0, y_range=fig.y_range,# min_border_top=0,
               toolbar_location="above", x_axis_location=None, y_axis_location=None,
               title="Charge ASIC versus time difference REF-ASIC")
    p.background_fill_color = "#fafafa"
    p.select(BoxSelectTool).select_every_mousemove = False
    p.select(LassoSelectTool).select_every_mousemove = False

    r = p.scatter('deltaT', 'charge_asic', source=src, size=3, color="#3A5785", alpha=0.6)
    r.data_source.selected.on_change('indices', update)
    return p


LINE_ARGS = dict(color="#3A5785", line_color=None)
colorScheme = plt.rcParams['axes.prop_cycle'].by_key()['color']
listColorScheme = [colorScheme[cluster % 8] for cluster in range(2)]
listLegend = ['selected', 'rest']
binStep = {'charge_ref':1.0 , 'charge_asic':1.0, 'deltaT':1.0} #'area_trigger_asic':1.0, 'areaPALG_asic':1.0}
measures = ['amplitude_ref', 'area_trigger_asic'] #'areaPALG_asic'] #

def plotTimingHisto(fig, globalHistoSrc, globalFitSrc):
    tmpFig = figure(plot_width=fig.width, plot_height=200, toolbar_location='right', min_border=0, #min_border_top=0, , y_axis_label='Count #'
                     x_range=fig.x_range, y_axis_location="left", x_axis_label='Timing [ns]', tools=TOOLS, background_fill_color="#fafafa")


    # tmpFig.output_backend = "svg"

    # tmpFig.quad(source=globalHistoSrc, bottom=0, top='count_0', left='left_0', right='right_0')
    #             # fill_color=colorScheme[1 % 8], fill_alpha=0.7, line_color="white")

    listHistoToStack = ['count_' + str(cluster) for cluster in range(1,-1,-1)]
    listBarCenterToStack = ['center_' + str(cluster) for cluster in range(1,-1,-1)]
    listBarWidthToStack = ['width_' + str(cluster) for cluster in range(1,-1,-1)]

    tmpFig.vbar_stack(listHistoToStack, x=listBarCenterToStack, width=listBarWidthToStack, source=globalHistoSrc,
                      fill_color=listColorScheme, fill_alpha=0.7, line_color=None) #"white",  alpha=0.5)
    # tmpFig.legend.orientation = "horizontal"
    # tmpFig.legend.click_policy = "hide"
    # tmpFig.line(source=globalFitSrc, x='xfit', y='yfit', color='red', line_width=2)

    return tmpFig

def plotTimingHistoWithFit(fig, globalHistoSrc, globalFigSrc):
    cluster = 1 # selected

    tmpFig = figure(plot_width=fig.width, plot_height=400, y_axis_location="left", x_axis_label='Timing [ns]', x_range=fig.x_range, min_border=0, #, y_axis_label='Count #'
                    tools=TOOLS, background_fill_color="#fafafa", active_scroll ="wheel_zoom")
    # tmpFig.toolbar.active_scroll = tmpFig.select_one(WheelZoomTool())
    # tmpFig.add_tools(WheelZoomTool(zoom_on_axis=False))
    tmpFig.output_backend = "svg"


    tmpFig.vbar(top='count_' + str(cluster), x='center_' + str(cluster), width='width_' + str(cluster),
                source=globalHistoSrc) # , fill_alpha=0.7, line_color="white"
    h0 = tmpFig.line(source=globalFigSrc, x='xfit', y='yfit', color='red', line_width=2)
    tmpFig.add_layout(Legend(items=[("selected", [h0]),]))

    return tmpFig


def plotGlobalRefHisto(fig, globalHistoSrc):
    tmpFig = figure(plot_width=fig.width, plot_height=200, toolbar_location=None, min_border=0, #min_border_top=0,
                    x_range=fig.x_range, y_axis_location="left", x_axis_label='Amplitude Reference [mV]',
                    background_fill_color="#fafafa")
    # tmpFig.output_backend = "svg"

    listHistoToStack = ['count_' + str(cluster) for cluster in range(1,-1,-1)]
    listBarCenterToStack = ['center_' + str(cluster) for cluster in range(1,-1,-1)]
    listBarWidthToStack = ['width_' + str(cluster) for cluster in range(1,-1,-1)]

    tmpFig.vbar_stack(listHistoToStack, x=listBarCenterToStack, width=listBarWidthToStack, source=globalHistoSrc,
                      fill_color=listColorScheme, fill_alpha=0.7, line_color=None) #, line_color="white", alpha=0.5)
    # tmpFig.legend.orientation = "horizontal"

    # color="white", line_color="#3A5785")
    return tmpFig

def plotGlobalASICHisto(fig, globalHistoSrc):
    tmpFig = figure(plot_width=200, plot_height=fig.height, toolbar_location=None, y_range = fig.y_range, min_border=0,# min_border_right = 0, y_axis_label = 'ASIC TOT'
                background_fill_color="#fafafa")
    # tmpFig.output_backend = "svg"

    listHistoToStack = ['count_' + str(cluster) for cluster in range(1,-1,-1)]
    listBarCenterToStack = ['center_' + str(cluster) for cluster in range(1,-1,-1)]
    listBarWidthToStack = ['width_' + str(cluster) for cluster in range(1,-1,-1)]

    tmpFig.hbar_stack(listHistoToStack, y=listBarCenterToStack, height=listBarWidthToStack,
                      legend_label= listLegend, source=globalHistoSrc,
                      fill_color = listColorScheme, fill_alpha=0.7, line_color=None) #, line_color="white", alpha=0.5)

    return tmpFig

def makeGlobalChargeHisto(): # For all clusters
    # return "all cluster in one" ColumnDataSource, only one for using in vbar_stack (stacked histograms)
    global data, binChargeAsic

    histo_df_list = []

    for charge in measures:
        histo_dict = {}
        if charge in data:
            if data[charge].dropna().empty:
                print("No data for "+charge)
                # histo_dict.update({'count'+suffix: None, 'left'+suffix: None, 'right'+suffix: None})
            else:
                if charge == 'amplitude_ref':
                    chargeType = 'charge_ref'
                else:
                    chargeType = 'charge_asic'
                bins = np.arange(data[charge].min(), data[charge].max() + .1, binStep[chargeType])
                need2refreshBinVar = False
                if len(bins)<10 or len(bins)>5000:
                    bins = "auto" # let numpy determine the bin size
                    need2refreshBinVar = True
                # print(charge, data[charge].min(), data[charge].max(), binStep[charge])
                for cluster in range(2):
                    clusterSelection = (data['cluster'] == cluster)

                    count, edges = np.histogram(data.loc[clusterSelection, charge], density=False, bins=bins)

                    if need2refreshBinVar:
                        newBining = edges[1]-edges[0]
                        binStep[charge] = newBining
                        bins = edges  # update for second round using the same bining as the first
                        binChargeAsic.update(value="%.3f" % newBining)
                        need2refreshBinVar = False

                    suffix = '_' + str(cluster)
                    histo_dict.update({'count'+suffix: count,
                                      'center'+suffix: (edges[:-1]+edges[1:])/2,
                                       'width'+suffix: (edges[1:]-edges[:-1])})
        else:
            print("No data for " + charge)

        histo_df = pd.DataFrame(histo_dict)
        # print(histo_df.describe())
        histo_df_list.append(ColumnDataSource(histo_df))

    # print(data.columns)
    return histo_df_list[0], histo_df_list[1]


def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

def half_max_x(x, y):
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    if len(zero_crossings_i) ==2:
        return lin_interp(x, y, zero_crossings_i[0], half), lin_interp(x, y, zero_crossings_i[1], half)
    else:
        return min(x),max(x)


def makeGlobalTimingHisto():
    histo_dict = {}
    legendLatex=legend=''
    global data
    fit_df = pd.DataFrame(columns=['xfit', 'yfit'])
    if data['deltaT'].dropna().empty:
        print("No deltaT data")
        suffix = '_1'
        # histo_df = pd.DataFrame(columns=['count_0', 'left_0', 'right_0'])
        histo_df = pd.DataFrame(columns=['count_0', 'center_0', 'width_0'])
    # if len(data.loc[data['cluster'] == 1, 'deltaT'])==0:
    #     print("No deltaT data")
    #     suffix = '_1'
    #     histo_df = pd.DataFrame(columns=['count'+suffix, 'left'+suffix, 'right'+suffix])
    else:
        tmp = data[data['cluster'] == 1]
        if len(tmp)>1:
            bins = np.arange(tmp['deltaT'].min(), tmp['deltaT'].max() + .1, binStep['deltaT'])
        else:
            bins = np.arange(data['deltaT'].min(), data['deltaT'].max() + .1, binStep['deltaT'])

        for cluster in range(2):
            clusterSelection = (data['cluster'] == cluster)
            count, edges = np.histogram(data.loc[clusterSelection, 'deltaT'], density=False, bins=bins)

            suffix = '_' + str(cluster)
            histo_dict.update({'count' + suffix: count,#'left'+suffix: edges[:-1], 'right'+suffix: edges[1:]})
                               'center' + suffix: (edges[:-1] + edges[1:]) / 2,
                               'width' + suffix: (edges[1:] - edges[:-1])})

            if cluster == 1 and autoFit : # selected
                entries = len(data.loc[clusterSelection, 'deltaT'])
                if entries > 1:  # nothing selected yet
                    count = np.append(count, 0)
                    if fitMethod == 'Gaussian':
                        try: # Fit gaussian timing difference
                            fitParam, cov = opt.curve_fit(modelGauss, edges, count,
                                                          p0=[edges[np.argmax(count)], .01, max(count)])
                            meanHisto = fitParam[0]
                            sigmaHisto = fitParam[1]
                            ampHisto = fitParam[2]
                            # xmin = min(edges)
                            # xmax = max(edges)
                            xmin = meanHisto - 5*sigmaHisto
                            xmax = meanHisto + 5*sigmaHisto

                            xfit = np.arange(xmin, xmax, 1e-3)
                            yfit = modelGauss(xfit, meanHisto, sigmaHisto, ampHisto)

                            r1, r2 = half_max_x(xfit, yfit)

                            sigmaHisto *= 10 ** 3
                            fwhmHisto = abs(r2 - r1) * 10 ** 3
                            legendLatex = ''.join((
                                '<p>Entries %.0f</p>' % (entries,),
                                r'<p>$$\mu$$ %.3f ns</p>' % (meanHisto,),
                                r'<p>$$\sigma$$ %.0f ps</p>' % (sigmaHisto,),
                                '<p>FWHM %.0f ps</p>' % (fwhmHisto,)
                            ))
                            legend = ''.join((
                                'Entries: %.0f \n' % (entries,),
                                'Mean:%.3f ns\n' % (meanHisto,),
                                'Sigma:%.0f ps\n' % (sigmaHisto,),
                                'FWHM:%.0f ps\n' % (fwhmHisto,)
                            ))

                            fit_df = pd.DataFrame({'xfit': xfit, 'yfit': yfit})

                            res_dict['sigma'] = "%.1f" % sigmaHisto
                            res_dict['entries'] = entries
                            res_dict['fwhm'] = "%.1f" % fwhmHisto
                            # res_dict['CTR_asic'] = "%.1f" % ((2 * (fwhmHisto ** 2) - 65 ** 2) ** 0.5)
                            res_dict['CTR_asic_sigma'] = "%.1f" % (fwhmHisto / (2 * (entries - 1)) ** 0.5)
                        except RuntimeError:
                            print("ERROR fit Gaussian")

                    elif fitMethod=='Gaussian+Exponential':
                        try:  # Fit gaussian timing
                            expectedMean = edges[np.argmax(count)]
                            mask = (edges>(expectedMean-2)) & (edges<(expectedMean+2))
                            edges = edges[mask]
                            count = count[mask]
                            fitParam, cov = opt.curve_fit(modelGaussExp, edges, count,
                                                          p0=[expectedMean, .01, max(count) / 2.718, 50])

                            meanHisto = fitParam[0]
                            sigmaHisto = fitParam[1]
                            ampHisto = fitParam[2]
                            lambdHisto = fitParam[3]

                            xmin = meanHisto - 20 * sigmaHisto
                            xmax = meanHisto + 20 * sigmaHisto

                            xfit = np.arange(xmin, xmax, 1e-3)
                            yfit = modelGaussExp(xfit, meanHisto, sigmaHisto, ampHisto, lambdHisto)

                            r1, r2 = half_max_x(xfit, yfit)

                            sigmaHisto *= 10 ** 3
                            fwhmHisto = abs(r2 - r1) * 10 ** 3
                            legendLatex = ''.join((
                                '<p>Entries %.0f</p>' % (entries,),
                                r'<p>$$\mu$$ %.3f ns</p>' % (meanHisto,),
                                r'<p>$$\sigma$$ %.0f ps</p>' % (sigmaHisto,),
                                '<p>FWHM %.0f ps</p>' % (fwhmHisto,),
                                r'<p>$$\lambda$$ %.0f ps</p>' % (lambdHisto,)

                            ))
                            legend = ''.join((
                                'Entries: %.0f \n' % (entries,),
                                'Mean:%.3f ns\n' % (meanHisto,),
                                'Sigma:%.0f ps\n' % (sigmaHisto,),
                                'FWHM:%.0f ps\n' % (fwhmHisto,),
                                'Lambda:%.0f ps\n' % (lambdHisto,)
                            ))

                            fit_df = pd.DataFrame({'xfit': xfit, 'yfit': yfit})

                            res_dict['sigma'] = "%.1f" % sigmaHisto
                            res_dict['fwhm'] = "%.1f" % fwhmHisto
                            res_dict['lambda'] = "%.1f" % lambdHisto
                            # res_dict['CTR_asic'] = "%.1f" % ((2 * (fwhmHisto ** 2) - 72 ** 2) ** 0.5)
                            res_dict['CTR_asic_sigma'] = "%.1f" % (fwhmHisto / (2 * (entries - 1)) ** 0.5)
                        except RuntimeError:
                            print("ERROR fit Gaussian+Exp")
                    else:
                        print('Undefined fit method')

        histo_df = pd.DataFrame(histo_dict)


    # print(histo_df.describe())
    return ColumnDataSource(histo_df), ColumnDataSource(fit_df), legendLatex, legend


def modelGauss(x, mu, sigma, amp):
    return amp * (1 / (sigma * (np.sqrt(2 * np.pi)))) * (np.exp((-1.0 / 2.0) * (((x - mu) / sigma) ** 2)))

def modelGaussExp(x, mu, sigma, amp, lambd):
    return amp * (lambd / 2 * np.exp(lambd / 2 * (2 * mu + lambd * sigma * sigma - 2 * x)) * (1 - erf((mu + lambd * sigma * sigma - x) / (1.414213562 * sigma))))
    # return amp * (lambd / 2 * np.exp(lambd / 2 * (2 * mu + lambd * sigma * sigma - 2 * np.log(x))) * (1 - erf((mu + lambd * sigma * sigma - np.log(x)) / (1.414213562 * sigma))))

def redrawTiming(withFitZoom=True):
    global globalHistoSrcTiming, globalFitSrc, p3Legend, p3, autoFit

    new_globalHistoSrcTiming, new_globalFitSrc, legendLatex, legend = makeGlobalTimingHisto()
    globalFitSrc.data.update(new_globalFitSrc.data)
    globalHistoSrcTiming.data.update(new_globalHistoSrcTiming.data)

    # Zoom
    if withFitZoom:
        if len(globalFitSrc.data['xfit']) != 0:
            p3.x_range.update(start=globalFitSrc.data['xfit'].min(), end=globalFitSrc.data['xfit'].max())
            fitHeigth = globalFitSrc.data['yfit'].max()
            p3.y_range.update(start=-0.05 * fitHeigth, end=1.2 * fitHeigth)
            p3Legend.text = legendLatex
            p3.legend.update(items=[(legend, p3.legend.items[0].renderers), ])


def redraw(withFitZoom=True):
    global globalHistoSrcRef, globalHistoSrcAsic, globalHistoSrcTiming, globalFitSrc
    new_globalHistoSrcRef, new_globalHistoSrcAsic = makeGlobalChargeHisto()
    globalHistoSrcRef.data.update(new_globalHistoSrcRef.data)
    globalHistoSrcAsic.data.update(new_globalHistoSrcAsic.data)

    redrawTiming(withFitZoom)


    # p2.y_range.update(p.y_range)
    # CustomJS(code="""document.querySelectorAll('.bk-tool-icon-reset[title="Reset"]').forEach(d => d.click())""")
    # CustomJS(args=dict(p=p2), code="""p.reset.emit()""")
    # p2 = plotSecondFig(p, update)

    print("Redrawn")

def update(attr, old, new):
    global data, p3

    indices = src.selected.indices
    clusterSelection = (data['cluster'] == 1)
    data.loc[clusterSelection, 'cluster'] = 0
    data.loc[indices, 'cluster'] = 1

    src.data.update(data)
    redraw()


def loadFiles(dataPath, deltaT_file):
    global data, src, measures, btnSave

    amplitude_ref_file = deltaT_file.replace("F3", "F1")
    # deltaT_file = deltaT_file.replace("F1", "F3")
    # deltaT_file = areaPALG_asic_file.replace("M2", "M3")
    # amplitude_ref_file = areaPALG_asic_file.replace("M2", "M1")

    area_trigger_file = deltaT_file.replace("F3", "F2")
    areaPALG_asic_file = deltaT_file.replace("F3", "F4")

    csvSeparator = '\t' #',' #

    areaPALG_asic_df = pd.read_csv('data/'+dataPath+'/'+areaPALG_asic_file, skiprows=4, encoding="ISO-8859-1", sep=csvSeparator, na_values="-nan(ind)")
    areaPALG_asic_df = areaPALG_asic_df.dropna().rename(columns={'Time':'ID', 'Ampl': 'areaPALG_asic'})
    areaPALG_asic_df['areaPALG_asic'] *= 10**9

    amplitude_ref_df = pd.read_csv('data/'+dataPath+'/'+amplitude_ref_file, skiprows=4, encoding="ISO-8859-1", sep=csvSeparator, na_values="-nan(ind)")
    amplitude_ref_df = amplitude_ref_df.dropna().rename(columns={'Time':'ID', 'Ampl': 'amplitude_ref'})
    amplitude_ref_df['amplitude_ref'] *= 10**9

    deltaT_df = pd.read_csv('data/'+dataPath+'/' + deltaT_file, skiprows=4, encoding="ISO-8859-1", sep=csvSeparator, na_values="-nan(ind)")
    deltaT_df = deltaT_df.dropna().rename(columns={'Time': 'ID', 'Ampl': 'deltaT'})
    deltaT_df['deltaT'] *= -10 ** 9  # fit quality
    # deltaT_df['deltaT'] += abs(deltaT_df['deltaT'].min()) -82 # shift all to avoid error when fitting with gaussian + exponential

    area_trigger_df = pd.read_csv('data/'+dataPath+'/' + area_trigger_file, skiprows=4, encoding="ISO-8859-1", sep=csvSeparator, na_values="-nan(ind)")
    area_trigger_df = area_trigger_df.dropna().rename(columns={'Time': 'ID', 'Ampl': 'area_trigger_asic'})
    area_trigger_df['area_trigger_asic'] *= 10**9

    data = pd.merge(areaPALG_asic_df, amplitude_ref_df, how='inner', on="ID")
    data = pd.merge(data, deltaT_df, how='inner', on="ID")
    data = pd.merge(data, area_trigger_df, how='inner', on="ID")

    data['cluster'] = 0
    chargeAsicCurSelected = measures[1]
    data['charge_asic'] = data[chargeAsicCurSelected]
    res_dict["charge_asic"] = chargeAsicCurSelected

    # src.data.update(data)
    src.data = src.from_df(data)
    # src = ColumnDataSource(data)

    # btnSave.button_type = 'primary'

def save():
    global deltaT_file, btnSave

    outputFolder = './'+outputPath+'/'+dataPath+'/'
    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)

    # filename= re.findall(r'--_(.*?)(?=.txt)', areaPALG_asic_file)[0]
    # filename= re.findall(r'--(.*?)(?=--[0-9]+.txt)', deltaT_file)[0]
    filename=deltaT_file
    globalPNG = outputFolder+filename+'.png'
    # detailSVG = outputFolder+filename+'_zoom.svg'
    detailPNG = outputFolder+filename+'_zoom.png'
    if os.path.isfile(globalPNG): os.remove(globalPNG)
    # if os.path.isfile(detailSVG): os.remove(detailSVG)
    if os.path.isfile(detailPNG): os.remove(detailPNG)
    export_png(doc, filename=globalPNG)
    # export_svg(p3, filename=detailSVG)
    export_png(p3, filename=detailPNG)

    # print(res_dict)
    res_dict["timeSaved"]=datetime.datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
    res_df = pd.DataFrame(res_dict, index=[0],
                          columns=sipmCaracList+['fwhm',  'CTR_asic_sigma', 'sigma', 'lambda',  'entries', 'charge_asic',
                                   'filename', 'folder', 'timeSaved', 'comment'])
    # 'discri_bias', 'curConv', 'fb', 'paOutBias', 'gain_asic','paBW_asic','comp_asic', 'Hyst_asic',  'oscilloscopeReferenceTimingHysteresis', 'CTR_asic',

    # Read if already existing result file
    csv_file = outputFolder+'results.csv'
    if os.path.isfile(csv_file):
        res_df = pd.concat([pd.read_csv(csv_file), res_df])
        os.remove(csv_file)
    res_df.to_csv(csv_file, index=False)

    dataCSV_file = outputFolder + filename + '.csv'
    data[['ID','areaPALG_asic','amplitude_ref','deltaT','area_trigger_asic','cluster']].to_csv(dataCSV_file)

    btnSave.button_type = 'success'
    print('Saved')
    # browser = webbrowser.get('firefox')
    # browser.open(detailSVG)

def extractInformationFromFilename(dataPath, filename):
    global res_dict

    # folder name
    # 20221010_HFsetup_NUVHDRH_UHDDA_vs_Reference
    print(dataPath)
    sipmModel = (re.findall(r'HFsetup(.*)_vs', dataPath)[0])

    # filename
    # C1--_45VNUVHDRHvs37VRef_Thresh100mV_80mV_--00000
    HV_SiPM_HF = findNumberIfAny('--_', filename)
    thres_HF = findNumberIfAny('Thresh', filename)
    HV_SiPM_ref = findNumberIfAny('vs', filename)



    # folder name
    # # 10102022_Liroc_PZC8_Vth460_noProbePA_allMasked_ch58_50ohm_noExtPZ_EPIC2x2x3_PbF2_black_FBK_NUVHDRH_UHDDA_vs_HF_TAC2x2x3_LYSOCeCa_4FP_BRCM_Nr3
    # asicChannel = findNumberIfAny('ch', dataPath)
    # externalResAtChannelInput = findNumberIfAny('ohm', dataPath, numberBefore=True)
    # if re.search(r"noExtPZ", dataPath) is not None:
    #     extPZC_circuit = 'None'
    # elif re.search(r"StefanPZ", dataPath) is not None:
    #     extPZC_circuit = 'StefanPZC'
    # else:
    #     extPZC_circuit = 'Unknown'
    #
    # crystalBrand = ['EPIC', 'TAC']
    # print('Crystal brand looked for in folder name : ', crystalBrand, ' Add if missing any !')
    # for brand in crystalBrand:
    #     detectorASICside = re.findall(r'' + brand + '(.*?)_vs', dataPath)
    #     if detectorASICside:
    #         detectorASICside = brand + detectorASICside[0]
    #         break
    #
    # detectorREFside = (re.findall(r'vs_HF_(.*)', dataPath)[0]).replace('_hysteresis', '')
    #
    # # F1_Liroc_PZC8_Th460_noProbePA_allMasked_vs_HF_37V_80mV_00000
    # HV_SiPM_asic = findNumberIfAny('mV_', filename)
    # threshold_asic = findNumberIfAny('Th', filename)
    #
    # PZC_asic = findNumberIfAny('PZC', filename)
    # # paBW_asic = findNumberIfAny('paBW', areaPALG_asic_file)
    # Hyst_asic = findNumberIfAny('Hyst', filename)
    # noProbePA = False if (re.search(r"noProbePA", dataPath) is None) else True
    # allMasked = False if (re.search(r"allMasked", dataPath) is None) else True
    #
    #
    # HV_SiPM_ref = int(re.findall(r'(?<=HF_)[0-9]+', filename)[0])
    # threshold_ref = int(re.findall(r'(?<=V_)[0-9]+(?=mV)', filename)[0])


    # # Elsaroc specific
    # comp_asic = findNumberIfAny('Comp', areaPALG_asic_file)
    # gain_asic = findNumberIfAny('gain_asic', areaPALG_asic_file)
    # discri_bias = findNumberIfAny('discriBias', areaPALG_asic_file)
    # curConv = findNumberIfAny('curConv', areaPALG_asic_file)
    # fb = findNumberIfAny('fb', areaPALG_asic_file)
    # paOutBias = findNumberIfAny('paOutBias', areaPALG_asic_file)
    #
    # asicPolarityUnticked = False if (re.search(r"PolarityUnticked", dataPath) is None) else True

    # oscilloscopeReferenceTimingHysteresis = False if (re.search(r"hysteresis", dataPath) is None) else True

    # 'discri_bias', 'curConv', 'fb', 'paOutBias',  'comp_asic', 'gain_asic', 'Hyst_asic','paBW_asic', 'oscilloscopeReferenceTimingHysteresis',


    for variableName in sipmCaracList:
        res_dict[variableName] = locals()[variableName]

    res_dict['folder'] = dataPath
    res_dict['filename'] = filename

    print("Folder : "+dataPath)

def chargeRefBinStepSelection(attr, old, new):
    global binStep
    binStep['charge_ref'] = float(new)
    redraw()

def chargeAsicBinStepSelection(attr, old, new):
    global binStep
    binStep['charge_asic'] = float(new)
    redraw()

def TimingBinStepSelection(attr, old, new):
    global binStep
    binStep['deltaT'] = float(new)
    redraw()

def commentToSaveInput(attr, old, new):
    global res_dict
    res_dict["comment"]=new

def gotoNextFile(attr, old, new):
    global i, p, deltaT_file, binTiming
    i = i + 1
    # selectData.
    deltaT_file = new
    print(new)
    # p.title.text = areaPALG_asic_file
    loadFiles(dataPath, deltaT_file)
    extractInformationFromFilename(dataPath, deltaT_file)
    binStep['deltaT'] = 1
    binTiming.update(value="1.0")
    # print(data.describe())
    redraw()
    # print(globalHistoSrcAsic.data)

def changeChargeASIC(attr, old, new):
    global measures, data
    if new == asicChargeList[0]:
        measures[1] = 'area_trigger_asic'
        data['charge_asic']=data['area_trigger_asic']
        # src.data.update(data)
        src.data = src.from_df(data)
        redraw(withFitZoom=False)
        res_dict["charge_asic"] = 'area_trigger_asic'

    elif new== asicChargeList[1]:
        measures[1] = 'areaPALG_asic'
        data['charge_asic']=data['areaPALG_asic']
        src.data = src.from_df(data)
        # src.data.update(data)
        redraw(withFitZoom=False)
        res_dict["charge_asic"] = 'areaPALG_asic'
    else:
        print("Charge asic variable not correctly defined")

def changeFitMethod(attr, old, new):
    global fitMethod
    fitMethod = new
    redrawTiming()

def changeFitAutoMode(attr, old, new):
    global autoFit
    autoFit = new
    redrawTiming()

def modify_doc(doc):
    global binTiming, binChargeAsic, btnSave

    # binChargeAsic = TextInput(value=str(binStep["charge_asic"]), width=70)
    binChargeAsic.on_change("value", chargeAsicBinStepSelection)
    # binChargeAsicLabel = Paragraph(text=" bin charge ASIC", width=70)

    binChargeRef = TextInput(value=str(binStep['charge_ref']), width=70)
    binChargeRef.on_change("value", chargeRefBinStepSelection)
    binChargeRefLabel = Paragraph(text=" bin charge REF", width=70)

    binTiming = TextInput(value=str(binStep['deltaT']), width=70)
    binTiming.on_change("value", TimingBinStepSelection)
    binTimingLabel = Paragraph(text=" bin timing", width=70)

    # btnNextData = Button(label='Next data', button_type='success', width = 100)
    # btnNextData.on_click(gotoNextFile)
    selectData = Select(title="Data in "+str(dataPath), value=fileList[0], options=fileList)
    selectData.on_change("value", gotoNextFile)

    selectASICcharge = Select(title="ASIC charge : ", value=asicChargeList[0], options=asicChargeList, width = 150)
    selectASICcharge.on_change("value", changeChargeASIC)

    selectFitMethod = Select(title="Fit method : ", value=fitMethodList[0], options=fitMethodList, width = 150)
    selectFitMethod.on_change("value", changeFitMethod)

    btnFitAuto = Toggle(label="Auto fit", button_type="warning", width = 100)
    btnFitAuto.on_change('active', changeFitAutoMode)

    btnSave = Button(label='Save', button_type='primary', width = 100)
    btnSave.on_click(save)

    commentToSaveLabel = Paragraph(text="Comment to save", width=150)
    commentToSave = TextInput(value='', width=squaredFigSize)
    commentToSave.on_change("value", commentToSaveInput)

    binSelect= column(selectASICcharge, row(binChargeAsic, binChargeAsicLabel), row(binChargeRef, binChargeRefLabel), row(binTiming, binTimingLabel))
    selectionArea = gridplot([[p, pv, p2],
                              [ph, binSelect, ph2],
                              [column(commentToSaveLabel, commentToSave), column(selectFitMethod, btnFitAuto, btnSave, p3Legend), p3]], merge_tools=False)
    layout = column(selectData, selectionArea)

    doc.add_root(layout)
    doc.title = "Selection Histogram"


i = 0
outputPath = 'outputBokeh'

# dataPath = '05052022_Liroc_ch22_EPIC2x2x3_PbF2_black_BRCM_Nr1_bias+_vs_HF_TAC2x2x3_LYSOCeCa_5FP_BRCM_Nr2_hysteresis'
# dataPath = '04052022_Liroc_ch10_50ohm_noExtPZ_EPIC2x2x3_PbF2_black_HBK_Nr79457_vs_HF_TAC2x2x3_LYSOCeCa_5FP_BRCM_Nr2_hysteresis'

# dataPath = '04052022_Liroc_ch2_50ohm_noExtPZ_EPIC2x2x3_PbF2_black_BRCM_Nr1_vs_HF_TAC2x2x3_LYSOCeCa_5FP_BRCM_Nr2_hysteresis' # crap
# dataPath = '04052022_Liroc_ch2_100ohm_noExtPZ_EPIC2x2x3_PbF2_black_BRCM_Nr1_vs_HF_TAC2x2x3_LYSOCeCa_5FP_BRCM_Nr2_hysteresis' # crap

# dataPath = '03052022_Liroc_ch8_StefanPZ_TAC2x2x3_LYSOCeCa_5FP_BRCM_Nr1_vs_HF_TAC2x2x3_LYSOCeCa_5FP_BRCM_Nr2_hysteresis'

# dataPath = '01052022_Radioroc_ch2_TAC2x2x3_LYSO_CeCa_5FP_BRCM_Nr1_vs_HF_TAC2x2x3_LYSO_CeCa_5FP_BRCM_Nr2'

# dataPath = '05052022_Elsaroc_ch6_EPIC2x2x3_PbF2_black_HBK_Nr79457_bias+_vs_HF_TAC2x2x3_LYSOCeCa_5FP_BRCM_Nr2_hysteresis'
# dataPath = '05052022_Elsaroc_ch6_EPIC2x2x3_PbF2_black_HBK_Nr79457_bias-_vs_HF_TAC2x2x3_LYSOCeCa_5FP_BRCM_Nr2_hysteresis'
# dataPath = '05052022_Elsaroc_ch6_EPIC2x2x3_PbF2_black_BRCM_Nr1_bias-_vs_HF_TAC2x2x3_LYSOCeCa_5FP_BRCM_Nr2_hysteresis'
# dataPath = '04052022_Elsaroc_ch6_EPIC2x2x3_PbF2_black_HBK_Nr79457_vs_HF_TAC2x2x3_LYSOCeCa_5FP_BRCM_Nr2_hysteresis_scan'
# dataPath = '1'
# dataPath = '03052022_Radioroc_ch0_EPIC2x2x3_BGO_5FP_FBKLF2M1_Nr28_vs_HF_TAC2x2x3_LYSOCeCa_5FP_BRCM_Nr2_hysteresis'
dataPath = '20221010_Reference/20221010_2_BRD_Nr3vsNr4_30um_TAC_PMI1X050vs051_LYSOCeCa_2x2x3_4FP_Teflon_Melt1.582_16C_5V_57mA_BGA2851Bal10nFPZ10pF680WHF1nF'
# dataPath = '10102022_Liroc_PZC8_Vth460_noProbePA_allMasked_ch58_50ohm_noExtPZ_EPIC2x2x3_PbF2_black_FBK_NUVHDRH_UHDDA_vs_HF_TAC2x2x3_LYSOCeCa_4FP_BRCM_Nr3'
# dataPath = '20221010_HFsetup_NUVHDRH_UHDDA_vs_Reference'

# selectionner les points à fitter plutôt que fit auto sur tout l'interval = rapidité
# corriger bug charge asic quand area pa ne change pas
# reset comment & zoom & selection when changing file

fileList = [filename for filename in os.listdir('data/'+dataPath+'/') if filename.startswith("F3")] # M2
deltaT_file = fileList[i]
asicChargeList = ['area trigger', 'TOT']
fitMethod = 'Gaussian'
fitMethodList = ['Gaussian', 'Gaussian+Exponential']
autoFit = False

loadFiles(dataPath, deltaT_file)
extractInformationFromFilename(dataPath, deltaT_file)

binChargeAsic = TextInput(value=str(binStep["charge_asic"]), width=70)
binChargeAsicLabel = Paragraph(text=" bin charge ASIC", width=70)

globalHistoSrcRef, globalHistoSrcAsic = makeGlobalChargeHisto()
globalHistoSrcTiming, globalFitSrc, _, _ = makeGlobalTimingHisto()

p = plotMainFig(update)
p2 = plotSecondFig(p, update)
ph = plotGlobalRefHisto(p, globalHistoSrcRef)
pv = plotGlobalASICHisto(p, globalHistoSrcAsic)
ph2 = plotTimingHisto(p2, globalHistoSrcTiming, globalFitSrc)
p3 = plotTimingHistoWithFit(p2, globalHistoSrcTiming, globalFitSrc)
p3Legend = Div(text='')

doc = curdoc()
modify_doc(doc)

# curdoc().add_root(layout)
# curdoc().title = "Selection Histogram"


