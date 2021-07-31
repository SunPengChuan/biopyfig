#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import matplotlib.pyplot as plt; plt.rcdefaults()
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import sys 
def fancybox(ax,loc,heigt,width,color,alpha):
    patches=[]
    fancybox = mpatches.FancyBboxPatch(
        loc, heigt, width,
        boxstyle=mpatches.BoxStyle("square", pad=0))
    patches.append(fancybox)
    collection = PatchCollection(patches, cmap=plt.cm.hsv, alpha=alpha,edgecolor="none",facecolor=color)
    ax.add_collection(collection)
	
def plot_bar(root,data,h,start_x,dict,color):
    for i in range(len(data)):
        print(data[i][0])
        k=dict[data[i][0]]
        heigt=h[k[0]]
        x=start_x+float(k[4])*step
        if float(data[i][1])<=40:
            col=color[0]
        elif float(data[i][1])<=50:
            col=color[1]
        elif float(data[i][1])<=60:
            col=color[2]
        elif float(data[i][1])<=70:
            col=color[3]
        else:
            col=color[4]
        plt.plot([x,x],[heigt,heigt+0.025],color=col,alpha=1,lw=0.5)

def read_gff(fn):
    f=open(fn)
    data,dict=[],{}
    for line in f.readlines():
        a=line.strip().split("\t")
        dict[a[1]]=[a[0],a[2],a[3],a[4],a[5]]
        data.append(a)
    return data,dict

def read_lens(fn):
    fp=open(fn)
    data=[]
    for row in fp.readlines():
        r1,r2=row.split()
        data.append([r1,r2])
    return data
def plot_chr(root,lens,gl,name):
    total_lens=max([float(k[1]) for k in lens])
    step=gl/float(total_lens)
    gl_start,n,start_x=0.85,0,0.08
    h={}
    align = dict(family='Times New Roman',style='normal',horizontalalignment="center", verticalalignment="center")
    for i in range(len(lens)):
        width=step*float(lens[i][1])
        heigt=gl_start-0.08*i
        h[lens[i][0]]=heigt
        loc=[start_x,heigt]
        fancybox(root,loc,width,0.025,'blue',0.2)
        plt.text(start_x-0.03,heigt+0.01,name[i],color='black',fontsize = 12,rotation = 0,weight='semibold',**align)        
    return step,h,start_x


if __name__=="__main__":
    plt.figure(figsize=(10, 6))
    root = plt.axes([0, 0, 1, 1])
    gff_1,dict_gff1=read_gff("Ad.order.gff")
    lens_1=read_lens("Ad_chrs.lens")
    gl=0.85
    name=['Ad'+k[0] for k in lens_1]
    step,h,start_x=plot_chr(root,lens_1,gl,name)
    data=read_lens("ad.conv.gc.txt")
    color=['blue','green','orange','red']
    plot_bar(root,data,h,start_x,dict_gff1,color)
    region=['<40','40-50','50-60','>70']
    plt.text(0.84,0.32,'GC(%)',fontsize=8)
    for i in range(len(color)):
        x=0.88
        y=0.3-i*0.04
        plt.plot([x,x+0.02],[y,y],color=color[i],lw=5)
        plt.text(x-0.04,y-0.01,region[i],fontsize=8)
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()
    #plt.show()
    plt.savefig("bar.png",dpi=600)
    plt.savefig("bar.pdf",dpi=600)