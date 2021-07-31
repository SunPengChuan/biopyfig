#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import matplotlib.pyplot as plt; plt.rcdefaults()
import sys
import re
import numpy as np


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


def plot_chr(data,height,step_x,step_y,start,r):
    data_line,dict_line=[],{}
    for i in range(len(data)):
        x=start+i*float(step_x)
        y=[height,height-float(data[i][1])*step_y]
        data_line.append(x)
        dict_line[data[i][0]]=x
        x1=x+2*r
        plt.plot([x,x],y,linestyle='-',color='black',lw=1.5)
        plt.plot([x1,x1],y,linestyle='-',color='black',lw=1.5)
        halfcicle(x+r,y,r)
    return data_line,dict_line

def halfcicle(x,y,r):
    tx=np.arange(0,np.pi,np.pi/180)
    ty=np.arange(np.pi,2*np.pi,np.pi/180)
    x1=float(x)+r*np.cos(tx)
    y1=float(y[0])+r*np.sin(tx)
    plt.plot(x1,y1,'k-',lw=1.5)
    x2=float(x)+r*np.cos(ty)
    y2=float(y[1])+r*np.sin(ty)
    plt.plot(x2,y2,'k-',lw=1.5)

def plot_name(data,line_1,height):
    align = dict(family='Times New Roman',style='normal',horizontalalignment="center", verticalalignment="center")
    for i in range(len(data)):
        x=line_1[i]+0.006
        plt.text(x,height,data[i],color='black',fontsize = 12,rotation = 0,weight='semibold',**align)

def read_gene(fn):
    f=open(fn)
    data=[]
    for line in f.readlines():
        a=line.strip().split("\t")
        for k in a:
            data.append(k)
    return data

def gene_loc(data,dict,line,height,step,r):
    data_loc,loc=[],[]
    data=list(set(data))
    for k in data:
        a=dict[k]
        x,y=line[a[0]],float(height)-float(a[4])*step
        plt.plot([x,x+2*r],[y,y],lw=0.5,color='black',alpha=0.5)
        loc.append([k,x,y])
    ids=list(set([k[1] for k in loc]))
    for k in ids:
        b=[a for a in loc if a[1] ==k]
        c=sorted(b,key=lambda x:float(x[2]),reverse=True)
        #print(c)       
        data_loc.append(c)
    return data_loc
def loc_new(loc,height,h):
    data=[k for k in loc]
    loc_new=[k for k in loc]
    for j in range(2):
        arr=[]
        print(loc)
        for i in range(len(loc)):
            if i< len(loc)-1 and loc[i]-loc[i+1] <= h+0.001:
                loc[i+1]=loc[i]-h
                arr.append(i)   
                continue
            else:
                arr.append(i)
                if arr[0] == 0 :
                    if i == len(loc)-1 and loc[i-1]-loc[i] <= h:
                        loc[i]=loc[i-1]-h
                    print(arr)
                    ave_old=sum([data[i] for i in arr])/len(arr)
                    ave_new=sum([loc[i] for i in arr])/len(arr)
                    step=ave_old-ave_new
                    if float(height)<float(step)+loc[arr[0]]:
                        step=height-loc[arr[0]]
                    for i in range(len(loc)):
                        if i in arr:
                            loc_new[i]=loc[i]+step
                    arr=[]
                    continue    
                lim1=loc[arr[0]-1]-loc[arr[0]]
                ave_old=sum([data[i] for i in arr])/len(arr)
                ave_new=sum([loc[i] for i in arr])/len(arr)
                step=ave_old-ave_new
                print(lim1,step)
                if float(lim1)<float(step):
                    step=lim1-h
                for i in range(len(loc)):
                    if i in arr:
                        loc_new[i]=loc[i]+step
                arr=[]
        if loc_new == loc :
            break
        else:
            loc=[k for k in loc_new]
    return loc_new
def chr_loc(data,height,h=0.012):
    for k in data:
        loc_chr=[a[2] for a in k]
        #print(loc_chr)
        loc_chr_new=loc_new(loc_chr,height,h)
        for i in range(len(k)):
        	x=float(k[i][1])-0.037
        	plt.plot([float(k[i][1]),float(k[i][1])-0.013],[float(k[i][2]),float(loc_chr_new[i])],lw=0.5,color='black')
        	plt.text(x,float(loc_chr_new[i]),k[i][0],color='blue',fontsize = 5.5,rotation = 0,**align)
        #print(loc_chr,loc_chr_new)

if __name__=="__main__":
    plt.figure(figsize=(12, 6))
    root = plt.axes([0, 0, 1, 1])
    align = dict(family='Times New Roman',style='normal',horizontalalignment="center", verticalalignment="center")
    gff_1,dict_gff_1=read_gff(sys.argv[1])  #gff_1
    lens_1=read_lens(sys.argv[2]) #lens_1
    name=['Os'+k[0] for k in lens_1]
    height,length=0.85,0.6
    start,r=0.1,0.006
    step_x=0.9/len(lens_1)
    step_y=length/max([float(k[1]) for k in lens_1])
    line_1,dict_line_1=plot_chr(lens_1,height,step_x,step_y,start,r)
    plot_name(name,line_1,height+0.05)
    data=read_gene(sys.argv[3])
    data_loc=gene_loc(data,dict_gff_1,dict_line_1,height,step_y,r)
    #print(data_loc)
    chr_loc(data_loc,height)
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()
    plt.savefig("chr_gene.png", dpi=300) 
