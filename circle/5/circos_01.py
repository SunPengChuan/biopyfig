#!/usr/bin/env python
# -*- coding: utf-8 -*- 

"""Read a colinearity file and store the data in a Python object
    Read a Gff file in order to obtain the gene location 
   
    colinearity file :
    the MCScanX out
    0-  0:        Os01g0584200    Os01g0713850      7e-15

    gff format:
    seqid	name	start	end	phase	order
    Os01    Os01g0100100    2449    9297    +       1
"""

import sys
import re
import numpy as np
import matplotlib.pyplot as plt; plt.rcdefaults()

def gff_all(fn):
    f=open(fn)
    data,dict=[],{}
    for line in f.readlines():
        line=line.strip()
        a=line.split("\t")
        a[0]=re.findall(r"\d+", a[0])
        if a[0]:
            a[0]=int(a[0][0])
        else:
            continue
        dict[a[1]]=[a[0],a[2],a[3],a[4],a[5]]
        data.append(a)
    return data,dict

def Bezier3(plist, t):    
    p0, p1, p2 = plist    
    return p0*(1-t)**2+2*p1*t*(1-t)+p2*t**2

def loc_real(rad, r):
    x=float(r)*np.cos(float(rad))
    y=float(r)*np.sin(float(rad))
    return x,y 

def gff_lens(data):
    chr=[k[0] for k in data]
    d = {k:chr.count(k) for k in set(chr)}
    newdata=[]
    chr=[]
    for k in data:
        if k[0] in chr:
            continue
        else:
            chr.append(k[0])
            newdata.append([k[0],d[k[0]]])
    return newdata

def plot_circle(lens,angle_gap,angle,radius,color='black',lw=1,alpha=1):
    start,end=0,0
    for i in range(0,len(lens)):
        end+=angle_gap+angle*(int(lens[i][1]))
        start=end-angle*(int(lens[i][1]))
        t=np.arange(start,end,0.005)
        x,y=(radius)*np.cos(t),(radius)*np.sin(t)
        plt.plot(x,y,linestyle='-',color=color, lw=lw,alpha=alpha)

def read_colinearity(fn):
    f1=open(fn)#
    data=[]
    flag=0
    for line in f1.readlines():
        line=line.strip()
        if re.match(r"^#",line):
            continue
        else:
            a=re.split(r"\t",line)
            data.append([a[1],a[2]])
    f1.close()
    return data

def gene_loction(lens,dict_gff,angle_gap,angle):
    dict1,chr_loc={},{}
    start,end=0,0
    for i in range(0,len(lens)):
        end+=angle_gap+angle*(int(lens[i][1]))
        start=end-angle*(int(lens[i][1]))
        chr_loc[lens[i][0]]=[start,end]
    for d,x in dict_gff.items():
        dict1[d]=float(chr_loc[x[0]][0])+int(x[4])*angle
    return dict1,chr_loc


def colinearity_loction(colinearity,gene_loc):
    colinearity_loc=[]
    for k in colinearity:
            colinearity_loc.append([gene_loc[k[0]],gene_loc[k[1]]])
    return colinearity_loc

def plot_colinearity(data,radius,color='gray',lw=0.02,alpha=1):
    for i in range(0, len(data)):
        ex1x, ex1y = radius*np.cos(float(data[i][0])), radius*np.sin(float(data[i][0]))
        ex2x, ex2y = radius*np.cos(float(data[i][1])), radius*np.sin(float(data[i][1]))
        x = [ex1x, 0.5*(ex1x+ex2x),ex2x]
        y = [ex1y,0.2*(ex1y+ex2y),ex2y]
        step = .01
        t = np.arange(0, 1+step, step)
        xt = Bezier3(x, t)
        yt = Bezier3(y, t)
        plt.plot(xt, yt, color=color, lw=lw,alpha=alpha)

def plot_labels(labels,chrlens_loc,radius,horizontalalignment="center", verticalalignment="center",fontsize = 6, color = 'black',):
    for i in range(len(chrlens_loc)):
        loc=sum(chrlens_loc[i])*0.5
        x,y=loc_real(loc,radius)
        if 1*np.pi<loc<2*np.pi:
            loc+=np.pi        
        plt.text(x,y,labels[i],horizontalalignment=horizontalalignment, verticalalignment=verticalalignment,fontsize = fontsize, color = color ,rotation = loc*180/np.pi-90)



if __name__=="__main__":
    plt.figure(figsize=(12, 12))
    root = plt.axes([0, 0, 1, 1])
    radius=0.6
    data,dict_gff=gff_all(sys.argv[1]) #gff
    lens=gff_lens(data)
    lens=sorted(lens, key=lambda x:int(x[0]))
    labels=['chr'+str(k[0]) for k in lens]
    angle_gap=0.08
    angle=(2*np.pi-(len(lens))*angle_gap)/(int(sum([k[1] for k in lens])))
    gene_loc,chr_loc=gene_loction(lens,dict_gff,angle_gap,angle)
    chrlens_loc=[chr_loc[k[0]] for k in lens]
    plot_circle(lens,angle_gap,angle,radius,lw=3)
    colinearity=read_colinearity(sys.argv[2])  #colinearity
    colinearity_loc=colinearity_loction(colinearity,gene_loc)
    plot_colinearity(colinearity_loc,radius-0.015,lw=0.05)
    plot_labels(labels,chrlens_loc,radius+0.05,fontsize=16)
    root.set_xlim(-1, 1)
    root.set_ylim(-1, 1)
    root.set_axis_off()
    plt.savefig("circos01.png", dpi=600)
