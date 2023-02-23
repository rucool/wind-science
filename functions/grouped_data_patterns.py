#!/usr/bin/env python

"""
Author: Laura Nazzaro on 1/12/2023 
Last modified: Laura Nazzaro 1/12/2023
Generate grouped summary statistics and boxplots of windspeed, power, and capacity factor
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
pd.set_option('display.width', 320, "display.max_columns", 15)  # for display in pycharm console



# boxplots: speed, capacity factor w/ power (https://www.python-graph-gallery.com/line-chart-dual-y-axis-with-matplotlib)
# heatmap: way to show amount of data at top and bottom (https://stackoverflow.com/questions/71417866/get-information-from-plt-hexbin)

def main(df, dv='speed', group=None, max_power=15000, boxplotFile=None, heatmapFile=None, csvFile=None, ttl=None, max_speed=40):
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
    season_order = ['winter', 'spring', 'summer', 'fall']
    month_order = [12,1,2,3,4,5,6,7,8,9,10,11]   
    t0=min(df['time']).strftime('%Y-%m-%d %H:%M')
    t1=max(df['time']).strftime('%Y-%m-%d %H:%M')
    if not ttl:
        ttl = f'{t0} to {t1}'   

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
        
        ovdf=pd.DataFrame()
        for a in order_vars[0]:
            ia = df[group[0]]==a
            xl = str(a)
            if len(group) > 1:
                for b in order_vars[1]:
                    if b==order_vars[1][0] and a != order_vars[0][0]:
                        xl=' '
                        ovdf = pd.concat([ovdf,pd.DataFrame({'order_var': [ov], 'xlabels': [xl], dv: [np.nan]})], axis=0, ignore_index=True)
                        xt=np.append(xt,ov)
                        xtl=np.append(xtl,xl)
                        ov+=1
                    ib = df[group[1]]==b
                    xl = ' '.join([str(a),str(b)])
                    if len(group) > 2:
                        for c in order_vars[2]:
                            if c==order_vars[2][0] and a != order_vars[0][0]:
                                xl=' '
                                ovdf = pd.concat([ovdf,pd.DataFrame({'order_var': [ov], 'xlabels': [xl], dv: [np.nan]})], axis=0, ignore_index=True)
                                xt=np.append(xt,ov)
                                xtl=np.append(xtl,xl)
                                ov+=1
                            ic = df[group[2]]==c
                            xl = ' '.join([str(a),str(b),str(c)])
                            i = np.logical_and(ia,np.logical_and(ib,ic))
                            df['order_var'][i]=ov
                            df['xlabels'][i]=xl
                            ovdf = pd.concat([ovdf,pd.DataFrame({'order_var': [ov], 'xlabels': [xl], dv: [np.nan]})], axis=0, ignore_index=True)
                            xt=np.append(xt,ov)
                            xtl=np.append(xtl,xl)
                            ov+=1
                    else:
                        i = np.logical_and(ia,ib)
                        df['order_var'][i]=ov
                        df['xlabels'][i]=xl
                        ovdf = pd.concat([ovdf,pd.DataFrame({'order_var': [ov], 'xlabels': [xl], dv: [np.nan]})], axis=0, ignore_index=True)
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
        df = pd.concat([df,ovdf],axis=0,ignore_index=True)

        fig, ax = plt.subplots(figsize=(12,7))
        df.boxplot(column=dv, by='order_var', ax=ax, notch=True, showmeans=True,
                            boxprops=dict(linewidth=1.25),
                            whiskerprops=dict(linewidth=1.25),
                            medianprops=dict(linewidth=5,color='black'),
                            flierprops=dict(marker='x'),
                            meanprops=dict(marker='^', markeredgecolor='black', markerfacecolor='black', markersize=15))
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
        plt.title(ttl)
        plt.suptitle('')
        ax.set_title('')
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