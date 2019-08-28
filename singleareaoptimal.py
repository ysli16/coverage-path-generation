# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from shapely import geometry
global velocity
global turnpanalty
global num

def rotate(point,sin,cos):
    x=cos*point[0]-sin*point[1]
    y=cos*point[1]+sin*point[0]
    rpoint=np.array([[x,y]]).reshape(1,2)
    return rpoint

def plot(rvertex,rwaypoint):
    plt.plot(rwaypoint[:,0],rwaypoint[:,1])
    poly = geometry.Polygon(rvertex)
    x,y = poly.exterior.xy
    plt.plot(x,y)
    
#estimate time comsumption   
def timeconsume(waypoint):
    length=0
    for i in range(len(waypoint)-1):
        length=length+math.sqrt(pow(waypoint[i][0]-waypoint[i+1][0],2)+pow(waypoint[i][1]-waypoint[i+1][1],2))
    turnnum=len(waypoint)-1
    return length/velocity+turnnum*turnpanalty

def rotategcoord(gcoord,angle):
    #rearrange coordinate(according to x)
    gcoord = gcoord[np.lexsort((gcoord[:,1],gcoord[:,0]))]
#convert to relative frame
    cgcoord=gcoord-gcoord[0]
    sin=math.sin(float(angle)/180*math.pi)
    cos=math.cos(float(angle)/180*math.pi)
    rcoord=np.zeros((0,2))
    for i in range(num):
        rpoint=rotate(cgcoord[i],sin,cos).reshape(1,2)
        rcoord=np.append(rcoord,rpoint,axis=0)
    coord=rcoord[np.lexsort((rcoord[:,1],rcoord[:,0]))]
    return coord
def coordtovertex(coord):   
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
            while(j<len(vertex)):
                if(vertex[j][0]<=coord[i][0]):
                    j=j+1
                else:
                    vertex=np.insert(vertex, j, values=coord[i], axis=0)
                    break
            if(j==len(vertex)):
                vertex=np.insert(vertex, j-1 , values=coord[i], axis=0)
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
    return vertex

def createwaypoint(gcoord,angle):
    #convert to rotated frame
    #angle=float(input("Enter the angle of the route(degree):"))
    coord=rotategcoord(gcoord,angle)    
    vertex=coordtovertex(coord)
    pivot=np.argwhere(vertex[:,0]==coord[num-1][0])[0][0]
    #define the width,generate cut position
#width=float(input("Enter width of the route:"))
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
        top=upsearchrange.max(0)[1]    
        buttom=downsearchrange.min(0)[1]
        
        newpoint=np.array([buttom,top])
        newpoint=newpoint.reshape(len(newpoint),1)
        if(i%2==1):
            newpoint=np.flipud(newpoint)
        newpoint=np.insert(newpoint,0,[pointpos[i]],axis=1)
        waypoint=np.append(waypoint,newpoint,axis=0)
    return waypoint



num=int(input("Enter the number of boundary points:"))
velocity=1
turnpanalty=30
gcoord=np.zeros((0,2))

for i in range(num):
    longitude_in=float(input("Enter boundary longitude:"))
    latitude_in=float(input("Enter boundary latitude:"))
    cin=np.array([[longitude_in,latitude_in]]).reshape(1,2)
    gcoord=np.append(gcoord,cin,axis=0)
    
width=float(input("Enter width of the route:"))
mintime=float('inf')
optimalangle=0
for angle in range(0,180,10):
    waypoint=createwaypoint(gcoord,angle)
    time=timeconsume(waypoint)
    if(time<mintime):
        mintime=time
        optimalangle=angle
for angle in range(optimalangle,optimalangle+10):
    waypoint=createwaypoint(gcoord,angle)
    time=timeconsume(waypoint)
    if(time<mintime):
        mintime=time
        optimalangle=angle
waypoint=createwaypoint(gcoord,optimalangle)
coord=gcoord[np.lexsort((gcoord[:,1],gcoord[:,0]))]
vertex=coordtovertex(coord)
#rotate back everything
rwaypoint=np.zeros(shape=[len(waypoint),2])
#rvertex=np.zeros(shape=[len(vertex),2])
sin=math.sin(float(-optimalangle)/180*math.pi)
cos=math.cos(float(-optimalangle)/180*math.pi)
for i in range(len(waypoint)):
    rwaypoint[i]=rotate(waypoint[i],sin,cos)
#for i in range(len(vertex)):
#    rvertex[i]=rotate(vertex[i],sin,cos)
rwaypoint=rwaypoint+coord[0]
#rvertex=rvertex+gcoord[0]

#plot the waypoints
plot(vertex,rwaypoint)


