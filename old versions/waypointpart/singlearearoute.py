# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from shapely import geometry
def rotate(point):
    x=cos*point[0]-sin*point[1]
    y=cos*point[1]+sin*point[0]
    rpoint=np.array([[x,y]]).reshape(1,2)
    return rpoint
def findcross(pointl,pointr,lineheight):
    x=(lineheight-pointl[1])/(pointr[1]-pointl[1])*(pointr[0]-pointl[0])+pointl[0]
    return x
num=int(input("Enter the number of boundry points:"))
#gcoord=np.zeros((1,0),dtype=[('longitude',float),('latitude',float)])
gcoord=np.zeros((0,2))
for i in range(num):
    longitude_in=float(input("Enter boundry longitude:"))
    latitude_in=float(input("Enter boundry latitude:"))
#    cin=np.array([(longitude_in,latitude_in)],dtype=[('longitude',float),('latitude',float)]).reshape(1,1)
    cin=np.array([[longitude_in,latitude_in]]).reshape(1,2)
    gcoord=np.append(gcoord,cin,axis=0)
#rearrange coordinate(according to x)
gcoord = gcoord[gcoord[:,0].argsort()]
#convert to relative frame
#startx=gcoord.min(0)[0]
#starty=gcoord.min(0)[1]
cgcoord=gcoord-gcoord[0]
#coord=preprocess(gcoord)
#convert to rotated frame
angle=float(input("Enter the angle of the route(degree):"))
sin=math.sin(angle/180*math.pi)
cos=math.cos(angle/180*math.pi)
rcoord=np.zeros((0,2))
for i in range(num):
#    x=math.cos(angle/180*math.pi)*coord[i][0]-math.sin(angle/180*math.pi)*coord[i][1]
#    y=math.cos(angle/180*math.pi)*coord[i][1]+math.sin(angle/180*math.pi)*coord[i][0]
#    rpoint=np.array([[x,y]]).reshape(1,2)
    rpoint=rotate(cgcoord[i]).reshape(1,2)
    rcoord=np.append(rcoord,rpoint,axis=0)
coord=rcoord[rcoord[:,0].argsort()]
#rearrange the sequence of boundry points so that they form a polygon
vertex=np.zeros((0,2))
vertex=np.append(vertex,coord[0].reshape(1,2),axis=0)
vertex=np.append(vertex,coord[num-1].reshape(1,2),axis=0)
for i in range(1,num-1):
    split=(coord[i][0]-coord[0][0])/(coord[num-1][0]-coord[0][0])*(coord[num-1][1]-coord[0][1])+coord[0][1]
    if(coord[i][1]>split):
        up=True
    else:
        up=False
    if(up):
        j=0
        while(vertex[j][0]<=coord[i][0]):
            j=j+1
        vertex=np.insert(vertex, j, values=coord[i], axis=0)
    else:
        j=np.argwhere(vertex[:,0]==coord[num-1][0])[0][0]
        while(j<len(vertex)):
            if(vertex[j][0]>=coord[i][0]):
                j=j+1
            else:
                vertex=np.insert(vertex, j, values=coord[i], axis=0)
                break
        if(j==len(vertex)):
            vertex=np.append(vertex, coord[i].reshape(1,2), axis=0)
pivot=np.argwhere(vertex[:,0]==coord[num-1][0])[0][0]
#define the width,generate cut position
width=float(input("Enter width of the route:"))
cutpos=np.arange(coord[0][0]+width,coord[num-1][0]+width,width)
pointpos=np.arange(coord[0][0]+width/2,coord[num-1][0]+width/2,width)
#find cross point with boundry
allvertex=vertex
upaddnum=0
for i in range(len(cutpos)-1):
    j=0
    k=pivot
    while(cutpos[i]>vertex[j][0]):
        j=j+1
    lefttop=vertex[j-1]
    righttop=vertex[j]
    while(k<len(vertex)):
        if(cutpos[i]<vertex[k][0]):
            k=k+1
        else:
            rightbuttom=vertex[k-1]
            leftbuttom=vertex[k]
            break
    if(k==len(vertex)):
        rightbuttom=vertex[k-1]
        leftbuttom=vertex[0]
    upperpos=[(cutpos[i],(cutpos[i]-lefttop[0])/(righttop[0]-lefttop[0])*(righttop[1]-lefttop[1])+lefttop[1])]
    lowerpos=[(cutpos[i],(cutpos[i]-leftbuttom[0])/(rightbuttom[0]-leftbuttom[0])*(rightbuttom[1]-leftbuttom[1])+leftbuttom[1])]
    allvertex=np.insert(allvertex, upaddnum+j, values=upperpos, axis=0)
    allvertex=np.insert(allvertex, upaddnum+k+1, values=lowerpos, axis=0)
    upaddnum=upaddnum+1
#find range of waypoints ineach cut
pointrange=np.zeros(shape=[0,2])
waypoint=np.zeros(shape=[0,2])
allvertex=np.append(allvertex,allvertex[0].reshape(1,2),axis=0)
allpivot=np.argwhere(allvertex[:,0]==coord[num-1][0])[0][0]
for i in range(len(cutpos)):
    if(i==0):
        leftupindex=0
        leftdownindex=len(allvertex)-1
    else:
        leftupindex=np.argwhere(allvertex[:,0]>=cutpos[i-1])[0][0]
        leftdownindex=np.argwhere(allvertex[allpivot:len(allvertex),0]<cutpos[i-1])[0][0]+allpivot-1
    if(i==len(cutpos)-1):
        rightupindex=allpivot
        rightdownindex=allpivot
    else:
        rightupindex=np.argwhere(allvertex[:,0]>cutpos[i])[0][0]-1
        rightdownindex=np.argwhere(allvertex[allpivot:len(allvertex),0]<cutpos[i])[0][0]+allpivot-1       
    upsearchrange=allvertex[leftupindex:rightupindex+1,:]
    downsearchrange=allvertex[rightdownindex:leftdownindex+1,:]   
    topmaxindex=upsearchrange.argmax(0)[1]
    topmin=upsearchrange.min(0)
    buttommax=downsearchrange.max(0)
    buttomminindex=downsearchrange.argmin(0)[1]
    if(topmin[1]-width/2>buttommax[1]+width/2):
        newpoint=np.array([topmin[1]-width/2,buttommax[1]+width/2])
    else:
        newpoint=np.array([(topmin[1]+buttommax[1])/2])
    newpoint=newpoint.reshape(len(newpoint),1)    
    newpoint=np.insert(newpoint,0,[pointpos[i]],axis=1)
    if(topmaxindex==0):
        topx1=upsearchrange[topmaxindex][0]
        topx2=findcross(upsearchrange[topmaxindex],upsearchrange[topmaxindex+1],upsearchrange[topmaxindex][1]-width/2)
    elif(topmaxindex==len(upsearchrange)-1):
        topx1=findcross(upsearchrange[topmaxindex-1],upsearchrange[topmaxindex],upsearchrange[topmaxindex][1]-width/2)
        topx2=upsearchrange[topmaxindex][0]
    else:
        topx1=findcross(upsearchrange[topmaxindex-1],upsearchrange[topmaxindex],upsearchrange[topmaxindex][1]-width/2)
        topx2=findcross(upsearchrange[topmaxindex],upsearchrange[topmaxindex+1],upsearchrange[topmaxindex][1]-width/2)
    if(topx1<=pointpos[i] and topx2>=pointpos[i]):
        newpoint=np.insert(newpoint,0,[[pointpos[i],upsearchrange[topmaxindex][1]-width/2]],axis=0)
    elif(topx1>pointpos[i]):
        newpoint=np.insert(newpoint,0,[[topx1,upsearchrange[topmaxindex][1]-width/2]],axis=0)
    else:
        newpoint=np.insert(newpoint,0,[[topx2,upsearchrange[topmaxindex][1]-width/2]],axis=0)
        
    if(buttomminindex==len(upsearchrange)-1):
        buttomx1=downsearchrange[buttomminindex][0]
        buttomx2=findcross(downsearchrange[buttomminindex],downsearchrange[buttomminindex-1],downsearchrange[buttomminindex][1]-width/2)
    elif(buttomminindex==0):
        buttomx1=findcross(downsearchrange[buttomminindex+1],downsearchrange[buttomminindex],downsearchrange[buttomminindex][1]-width/2)
        buttomx2=downsearchrange[buttomminindex][0]
    else:
        buttomx1=findcross(downsearchrange[buttomminindex+1],downsearchrange[buttomminindex],downsearchrange[buttomminindex][1]-width/2)
        buttomx2=findcross(downsearchrange[buttomminindex],downsearchrange[buttomminindex-1],downsearchrange[buttomminindex][1]-width/2)
    if(buttomx1<=pointpos[i] and buttomx2>=pointpos[i]):
        newpoint=np.append(newpoint,[[pointpos[i],downsearchrange[buttomminindex][1]-width/2]],axis=0)
    elif(topx1>pointpos[i]):
        newpoint=np.append(newpoint,[[buttomx1,downsearchrange[buttomminindex][1]-width/2]],axis=0)
    else:
        newpoint=np.append(newpoint,[[buttomx2,downsearchrange[buttomminindex][1]-width/2]],axis=0)
    if(i%2==1):
        newpoint=np.flipud(newpoint)
    
    waypoint=np.append(waypoint,newpoint,axis=0)
#rotate back everything
rwaypoint=np.zeros(shape=[len(waypoint),2])
rvertex=np.zeros(shape=[len(vertex),2])
sin=math.sin(-angle/180*math.pi)
cos=math.cos(-angle/180*math.pi)
for i in range(len(waypoint)):
    rwaypoint[i]=rotate(waypoint[i])
for i in range(len(vertex)):
    rvertex[i]=rotate(vertex[i])
rwaypoint=rwaypoint+gcoord[0]
rvertex=rvertex+gcoord[0]
#plot the waypoints
plt.plot(rwaypoint[:,0],rwaypoint[:,1])
poly = geometry.Polygon(rvertex)
x,y = poly.exterior.xy
plt.plot(x,y)
