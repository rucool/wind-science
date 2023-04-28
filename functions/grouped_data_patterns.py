#!/usr/bin/env python

"""
Author: Laura Nazzaro on 1/12/2023 
Last modified: Laura Nazzaro 3/22/2023
Generate grouped summary statistics and boxplots of windspeed, power, and capacity factor
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
import calendar
plt.rcParams.update({'font.size': 12})
pd.set_option('display.width', 320, "display.max_columns", 15)  # for display in pycharm console



# boxplots: speed, capacity factor w/ power (https://www.python-graph-gallery.com/line-chart-dual-y-axis-with-matplotlib)
# heatmap: way to show amount of data at top and bottom (https://stackoverflow.com/questions/71417866/get-information-from-plt-hexbin)

def main(df, dv='speed', group=None, max_power=15000, boxplotFile=None, heatmapFile=None, csvFile=None, ttl=None, max_speed=40, all_black=False):
    # define units
    unit_labels = dict(speed='m/s',power='kW',capacity_factor=' ')
    if max_power>50000:
        max_power = max_power/1000
        df['power'] = df['power']/1000
        unit_labels['power'] = 'MW'
        if max_power>1500:
            max_power = max_power/1000
            df['power'] = df['power']/1000
            unit_labels['power'] = 'GW'
    season_order = dict(winter=1,spring=2,summer=3,fall=4)
    heat_cmap = mpl.colormaps['Reds']
    heat_cmap.set_under('white')
    boxplot_colormap=cm.tab20
    boxplot_colororder=np.append(np.arange(0,20,2),np.arange(1,20,2))
    boxplot_colorlist=[boxplot_colormap(n) for n in boxplot_colororder]
    boxplot_colorlist_8=['#7b85d4','#f37738','#83c995','#d7369e','#c4c9d8','#859795','#e9d043','#ad5b50']
    boxplot_hatches=['none','/','x','.','///','xxx','...']
    season_order = ['winter', 'spring', 'summer', 'fall']
    month_order = [12,1,2,3,4,5,6,7,8,9,10,11]   
    if len(np.unique(df['month']))<7:
        month_order = np.unique(df['month'])
        if 12 in month_order:
            month_order=np.append(12,month_order[:-1])
        month_order=list(month_order)
    t0=min(df['time']).strftime('%Y-%m-%d %H:%M')
    t1=max(df['time']).strftime('%Y-%m-%d %H:%M')
    if not ttl:
        ttl = f'{t0} to {t1}'   

    shortnames={'upwelling': {'upwelling NA': 'upw-NA', 'upwelling absent': 'upw-N', 'upwelling present': 'upw-Y'},
        'turbines': {'wind farm': 'WF', 'control': 'ctrl'},
        'season': {'winter': 'Win', 'spring': 'Spr', 'summer': 'Sum', 'fall': 'Aut'}}
    shortlabels=False
    if len(group)>1:
        shortlabels=True

    if (boxplotFile or csvFile) and group:
        df['order_var'] = np.nan
        df['xlabels'] = np.nan
        ov=1
        xt=[]
        xtl=[]

        if not isinstance(group,list):
            group = [group]

        if len(group) > 3:
            print('Can only accept three group types. Exiting.')
            return 0
    
    if csvFile and group:
        all_vars = group.copy()
        all_vars.append(dv)
        stats = df[all_vars].groupby(group).describe()
        if dv=='power':
            sums = df[all_vars].groupby(group).sum()
            stats['total'] = sums['power']
        stats.to_csv(csvFile)
    
    df['doy'] = df['time'].dt.dayofyear
    df.loc[df['month']==12, 'doy'] -= 365
    df.loc[df['month']==12, 'year'] += 1

    if boxplotFile and group:
        order_vars = []
        for g in group:
            if g=='month':
                order_vars.append(month_order)
            elif g=='season':
                order_vars.append(season_order)
            else:
                order_vars.append(list(np.unique(df[g])))

        if len(group) > 1 and len(order_vars[1]) <= 8:
            boxplot_colorlist=boxplot_colorlist_8
        
        ovdf=pd.DataFrame()
        for a in order_vars[0]:
            ia = df[group[0]]==a
            alabel = str(a)
            if group[0]=='month':
                alabel = calendar.month_abbr[a]
            if shortlabels and group[0] in shortnames.keys():
                alabel = shortnames[group[0]][a]
            xl = alabel
            if len(group) > 1:
                for b in order_vars[1]:
                    blabel = str(b)
                    if group[1]=='month':
                        blabel = calendar.month_abbr[b]
                    if shortlabels and group[1] in shortnames.keys():
                        blabel = shortnames[group[1]][b]
                    if b==order_vars[1][0] and a != order_vars[0][0]:
                        xl=' '
                        ovdf = pd.concat([ovdf,pd.DataFrame({'order_var': [ov], 'xlabels': [xl], dv: [np.nan], 'color': [' '], 'hatch': [' ']})], axis=0, ignore_index=True)
                        xt=np.append(xt,ov)
                        xtl=np.append(xtl,xl)
                        ov+=1
                    ib = df[group[1]]==b
                    xl = ', '.join([alabel,blabel])
                    if len(group) > 2:
                        for c in order_vars[2]:
                            clabel = str(c)
                            if group[2]=='month':
                                clabel = calendar.month_abbr[c]
                            if shortlabels and group[2] in shortnames.keys():
                                clabel = shortnames[group[2]][c]
                            if c==order_vars[2][0] and a != order_vars[0][0]:
                                xl=' '
                                ovdf = pd.concat([ovdf,pd.DataFrame({'order_var': [ov], 'xlabels': [xl], dv: [np.nan], 'color': [' '], 'hatch': [' ']})], axis=0, ignore_index=True)
                                xt=np.append(xt,ov)
                                xtl=np.append(xtl,xl)
                                ov+=1
                            ic = df[group[2]]==c
                            xl = ', '.join([alabel,blabel,clabel])
                            i = np.logical_and(ia,np.logical_and(ib,ic))
                            df['order_var'][i]=ov
                            df['xlabels'][i]=xl
                            ovdf = pd.concat([ovdf,pd.DataFrame({'order_var': [ov], 'xlabels': [xl], dv: [np.nan], 'color': [boxplot_colorlist[order_vars[1].index(b)]], 'hatch': [boxplot_hatches[order_vars[2].index(c)]]})], axis=0, ignore_index=True)
                            xt=np.append(xt,ov)
                            xtl=np.append(xtl,xl)
                            ov+=1
                    else:
                        i = np.logical_and(ia,ib)
                        df['order_var'][i]=ov
                        df['xlabels'][i]=xl
                        ovdf = pd.concat([ovdf,pd.DataFrame({'order_var': [ov], 'xlabels': [xl], dv: [np.nan], 'color': [boxplot_colorlist[order_vars[1].index(b)]], 'hatch': ['none']})], axis=0, ignore_index=True)
                        xt=np.append(xt,ov)
                        xtl=np.append(xtl,xl)
                        ov+=1
            else:
                i = ia
                df['order_var'][i]=ov
                df['xlabels'][i]=xl
                ovdf = pd.concat([ovdf,pd.DataFrame({'order_var': [ov], 'xlabels': [xl], dv: [np.nan]})], axis=0, ignore_index=True)
                xt=np.append(xt,ov)
                xtl=np.append(xtl,xl)
                ov+=1
        df = pd.concat([df,ovdf[['order_var','xlabels',dv]]],axis=0,ignore_index=True)

        plt.rcParams.update({'font.size': 12})
        fig, ax = plt.subplots(figsize=(12,8))
        if all_black or len(group)==1:
            df.boxplot(column=dv, by='order_var', ax=ax, notch=True, showmeans=True,
                                boxprops=dict(linewidth=1.25),
                                whiskerprops=dict(linewidth=1.25),
                                medianprops=dict(linewidth=5,color='black'),
                                flierprops=dict(marker='x'),
                                meanprops=dict(marker='^', markeredgecolor='black', markerfacecolor='black', markersize=15))
        else:
            box = df.boxplot(column=dv, by='order_var', ax=ax, patch_artist=True, return_type='both', \
                                notch=True, showmeans=True,
                                boxprops=dict(linewidth=1.25),
                                whiskerprops=dict(linewidth=1.25),
                                medianprops=dict(linewidth=5,color='black'),
                                flierprops=dict(marker='x'),
                                meanprops=dict(marker='^', markeredgecolor='black', markerfacecolor='black', markersize=15))
            #for boxid, color in zip(box['boxes'], ovdf['color']):
            for row_key, (ax,row) in box.iteritems():
                for i,boxid in enumerate(row['boxes']):
                    if ovdf['color'][i] not in ['',' ','none']:
                        boxid.set_facecolor(ovdf['color'][i])
                    if ovdf['hatch'][i] not in ['',' ','none']:
                        boxid.set_hatch(ovdf['hatch'][i])
        plt.xticks(xt,xtl)
        plt.xticks(rotation=45, ha='right')
        plt.xlabel(None)
        plt.ylabel(dv + ' (' + unit_labels[dv] + ')')
        axis_extend = 0.005
        if dv=='speed':
            plt.ylim([0,max_speed])
        elif dv=='power':
            plt.ylim([-max_power*axis_extend,max_power*(1+axis_extend)])
            axb=ax.twinx()
            plt.ylim([-axis_extend,1+axis_extend])
            plt.ylabel('capacity_factor')
        elif dv=='capacity_factor':
            plt.ylabel(dv)
            plt.ylim([-axis_extend,1+axis_extend])
            axb=ax.twinx()
            plt.ylim([-max_power*axis_extend,max_power*(1+axis_extend)])
            plt.ylabel('power (' + unit_labels['power'] + ')')
        #plt.suptitle(ttl)
        plt.title('')
        plt.suptitle('')
        ax.set_title(ttl)
        plt.savefig(boxplotFile, dpi=300)
        plt.close()

    if heatmapFile:
        plt.rcParams.update({'font.size': 24}) 
        doy_ticks = pd.DataFrame(np.append(12,range(1,12)),columns=['month'])
        doy_ticks['day'] = 15
        doy_ticks['year'] = 2011
        doy_ticks['time'] = pd.to_datetime(doy_ticks)
        doy_ticks['doy'] = doy_ticks['time'].dt.dayofyear
        doy_ticks.loc[doy_ticks['month']==12, 'doy'] -= 365
        doy_ticks['labels']=['D','J','F','M','A','M','J','J','A','S','O','N']
        # doy_ticks['labels']=['Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov']
        season_breaks=pd.DataFrame({'date': ['2011-03-01','2011-06-01','2011-09-01']})
        season_breaks['date']=pd.to_datetime(season_breaks['date'])
        season_breaks['doy']=season_breaks['date'].dt.dayofyear

        if dv=='speed':
            vm = len(df)/500
        else:
            vm=len(df)/1000
        fig, ax = plt.subplots(figsize=(12,12))
        hexcoll = plt.hexbin(df['doy'], df[dv], cmap=heat_cmap, 
                        gridsize=52, vmin=1, vmax=vm)
        hexx = hexcoll.get_offsets()[:,0]
        hexy = hexcoll.get_offsets()[:,1]
        hexc = hexcoll.get_array()
        i = np.where(hexc>vm)
        #min_i = np.where(np.logical_and(hexy==np.min(hexy),hexc>vm))
        #max_i = np.where(np.logical_and(hexy==np.max(hexy),hexc>vm))
        base_marker_size = 8
        axis_extend = 0.03
        plt.scatter(hexx[i], hexy[i], s=base_marker_size*(hexc[i]/vm)**2, c='blue', alpha=0.25)
        plt.scatter(hexx[i], hexy[i], s=base_marker_size, c='lightblue')
        #plt.scatter(hexx[min_i], hexy[min_i], s=base_marker_size*hexc[min_i]/vm, c='none', edgecolors='blue')
        #plt.scatter(hexx[min_i], hexy[min_i], s=base_marker_size, c='lightblue')
        #plt.scatter(hexx[max_i], hexy[max_i], s=base_marker_size*hexc[max_i]/vm, c='none', edgecolors='blue')
        #plt.scatter(hexx[max_i], hexy[max_i], s=base_marker_size, c='lightblue')
        plt.xticks(doy_ticks['doy'],doy_ticks['labels'])
        plt.xlabel(None)
        plt.ylabel(dv + ' (' + unit_labels[dv] + ')')
        if dv=='speed':
            plt.ylim([0,max_speed])
            for x in season_breaks['doy']:
                plt.plot([x,x],[0,max_speed],c='black')
        elif dv=='power':
            plt.ylim([-max_power*axis_extend,max_power*(1+axis_extend)])
            for x in season_breaks['doy']:
                plt.plot([x,x],[-max_power*axis_extend,max_power*(1+axis_extend)],c='black')
            axb=ax.twinx()
            plt.ylim([-axis_extend,1+axis_extend])
            plt.ylabel('capacity_factor')
        elif dv=='capacity_factor':
            plt.ylabel(dv)
            plt.ylim([-axis_extend,1+axis_extend])
            for x in season_breaks['doy']:
                plt.plot([x,x],[-axis_extend,1+axis_extend],c='black')
            axb=ax.twinx()
            axb.ylim([-max_power*axis_extend,max_power*(1+axis_extend)])
            axb.ylabel('power (' + unit_labels['power'] + ')')
        plt.title(ttl)
        plt.savefig(heatmapFile, dpi=300)
        plt.close()