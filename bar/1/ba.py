#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

def fancybox(ax,loc,heigt,width,color,alpha):
    patches=[]
    fancybox = mpatches.FancyBboxPatch(
        loc, height, width,
        boxstyle=mpatches.BoxStyle("square", pad=0))
    patches.append(fancybox)
    collection = PatchCollection(patches, cmap=plt.cm.hsv, alpha=alpha,edgecolor="none",facecolor=color)
    ax.add_collection(collection)

def plot_bar(root,data,color,name):
    for i in range(len(data)):
        start=0.3
        y=heigt-i*0.1+0.006
        plt.text(start-0.1,y,str(name[i]))
        for j in range(len(data[i])):
            loc=[start,heigt-0.1*i]
            width=float(data[i][j])/sum(data[i])
            fancybox(root,loc,width,0.025,color[j],0.5)
            x=start+0.5*width
            plt.text(x-0.02,y,str(round(width*100))+"%")
            start+=width

if __name__=="__main__":
    plt.figure(figsize=(10, 10))
    root = plt.axes([0, 0, 1, 1])
    name=['123','456']
    data=[[1,3,2,4,5],[1,3,2,4,5]]
    color=['red','blue','green','orange','lime']
    height=0.7
    plot_bar(root,data,color,name)
    root.set_xlim(0, 1.5)
    root.set_ylim(0, 1)
    root.set_axis_off()
    plt.savefig("bar.png",dpi=300)