#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:50:19 2019

@author: ysli
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import Point, Line, Polygon,Segment
import gurobipy as grb

def checkcrossline(polygon):
    sides=polygon.sides
    count=0
    for i in range(len(sides)-1):
        for j in range(i+1,len(sides)):
            if(len(sides[i].intersection(sides[j]))==1):
                count=count+1
    if(count==len(polygon.vertices)):
        return False
    else:
        return True
def formpolygon(coord):
    #rearrange coordinate(according to x)
    num=len(coord)
    coord=coord[np.lexsort((coord[:,1],coord[:,0]))]
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
def findpoint(points,goal):
    i=0
    index=[]
    for point in points:
        if(point.equals(goal)):
            index.append(i)
        else:
            i=i+1
    if(len(index)==1):
        return index[0]
    if(len(index)>1):
        return index
    if(len(index)==0):
        return -1
def splitpolygon(polygon,line):
    crosspoint=polygon.intersection(line)
    sides=polygon.sides
    vertices=polygon.vertices
    if(len(crosspoint)==2):
        if not(line.is_parallel(xaxis)):
            if(crosspoint[0].y>crosspoint[1].y):
                upcross=crosspoint[0]
                downcross=crosspoint[1]
            else:
                downcross=crosspoint[0]
                upcross=crosspoint[1]
            i=0
            for side in sides:
                if(side.contains(upcross)):
                    if(upcross.x<side.p1.x or(upcross.x==side.p1.x and upcross.x>side.p2.x)):                        
                        rightupindex=i
                        leftupindex=i+1
                    else:
                        leftupindex=i
                        rightupindex=i+1
                if(side.contains(downcross)):
                    if(downcross.x<side.p1.x or(downcross.x==side.p1.x and downcross.x>side.p2.x)):
                        rightdownindex=i
                        leftdownindex=i+1
                    else:
                        leftdownindex=i
                        rightdownindex=i+1
                i=i+1
            rightpt=list([upcross,downcross])
            if(rightdownindex<=rightupindex):
                rightpt.extend(vertices[rightdownindex:rightupindex+1])
            else:
                rightpt.extend(vertices[rightdownindex:len(vertices)])
                rightpt.extend(vertices[0:rightupindex+1])
            if(leftupindex<=leftdownindex):
                leftpt=list(vertices[leftupindex:leftdownindex+1])
            else:
                leftpt=list(vertices[leftupindex:len(vertices)])
                leftpt.extend(vertices[0:leftdownindex+1])
            leftpt.extend([downcross,upcross])
            leftpolygon=Polygon(*leftpt)
            rightpolygon=Polygon(*rightpt)
            return leftpolygon,rightpolygon
        else:
            if(crosspoint[0].x>crosspoint[1].x):
                rightcross=crosspoint[0]
                leftcross=crosspoint[1]
            else:
                rightcross=crosspoint[1]
                leftcross=crosspoint[0]
            for side in sides:
                if(side.contains(rightcross)):
                    if(rightcross.y>side.p1.y or(rightcross.y==side.p1.y and rightcross.y<side.p2.y)):
                        rightdownindex=findpoint(vertices,side.p1)
                        rightupindex=findpoint(vertices,side.p2)
                    else:
                        rightupindex=findpoint(vertices,side.p1)
                        rightdownindex=findpoint(vertices,side.p2)
                if(side.contains(leftcross)):
                    if(leftcross.y>side.p1.y or(leftcross.y==side.p1.y and leftcross.y<side.p2.y)):
                        leftdownindex=findpoint(vertices,side.p1)
                        leftupindex=findpoint(vertices,side.p2)
                    else:
                        leftupindex=findpoint(vertices,side.p1)
                        leftdownindex=findpoint(vertices,side.p2)
            uppt=list([leftcross,rightcross])
            if(rightupindex<=leftupindex):
                uppt.extend(vertices[rightupindex:leftupindex+1])
            else:
                uppt.extend(vertices[rightupindex:len(vertices)])
                uppt.extend(vertices[0:leftupindex+1])
            if(leftdownindex<=rightdownindex):
                downpt=list(vertices[leftdownindex:rightdownindex+1])
            else:
                downpt=list(vertices[leftdownindex:len(vertices)])
                downpt.extend(vertices[0:rightdownindex+1])
            downpt.extend([rightcross,leftcross])
            uppolygon=Polygon(*uppt)
            downpolygon=Polygon(*downpt)
            return uppolygon,downpolygon
    else:
        return polygon,None
def onpolygon(point,polygon):
    sides=polygon.sides
    for side in sides:
        if(side.contains(point)):
            return True
    return False
def containitem(item,xlist):
    for x in xlist:
        if(x==item):
            return True
    return False
def minuspolygon(polygon,obs,pt):
    verticesp=polygon.vertices
    sides=polygon.sides
    for side in sides:
        if(side.contains(Point(pt))):
            break       
    indexp1=findpoint(verticesp,side.p1)
    indexp2=findpoint(verticesp,side.p2)
    if(indexp1>indexp2):
        indexp=indexp1
    else:
        indexp=indexp2
    verticeso=obs.vertices
    verticeso.reverse()
    indexo=findpoint(verticeso,Point(pt))
    verticesp[indexp:indexp]=iter(verticeso[0:indexo+1])
    verticesp[indexp:indexp]=iter(verticeso[indexo:len(verticeso)])
    polygon=Polygon(*verticesp)
    return polygon
'''
area={
      1:[2,3],
      2:[4],
      3:[4]
      }
coord=[
       [(2,2),(5,2)],
       [(5,8),(9,11)],
       [(5,1),(10,1)],
       [(9,5),(13,3)]
       ]
'''
def addtomodel(cells,newcell,area,coord):
    i=1
    area[len(cells)+1]=[]
    for cell in cells:
        if(len(cell.intersection(newcell))):
            area[i].append(len(cells)+1)
        i=i+1
    a=newcell.sides[0].midpoint
    b=newcell.sides[-1].midpoint
    newcoord=[(a.x,a.y),(b.x,b.y)]
    coord.append(newcoord)
def distance(point1,point2):
    distance=math.sqrt(pow(point1[0]-point2[0],2)+pow(point1[1]-point2[1],2))
    return distance
def printsolution(m,pathvar):
    if m.status == grb.GRB.Status.OPTIMAL:
        print('totalCost: %g' % m.objVal)
        print('Path:')
        for f in pathvar:
            if (pathvar[f].X!=0):
                #print('%s' % (pathvar[f].VarName))
                if(f[2]==1):
                    print 'a%d to a%d' %(f[0],f[1])
                elif(f[2]==2):
                    print 'a%d to b%d' %(f[0],f[1])
                elif(f[2]==3):
                    print 'b%d to a%d' %(f[0],f[1])
                else:
                    print 'b%d to b%d' %(f[0],f[1])
    else:
        print('No solution')
def solveLP(area,coord):
    m=grb.Model('ILP')
    modes=[1,2,3,4]
    allpath=grb.tuplelist()
    cost={}
    for cell1 in area:
        for cell2 in area[cell1]:
            for mode in modes:
                allpath.append((cell1,cell2,mode))
                if(mode==1):
                    cost[(cell1,cell2,mode)]=distance(coord[cell1-1][0],coord[cell2-1][0])
                elif(mode==2):
                    cost[(cell1,cell2,mode)]=distance(coord[cell1-1][0],coord[cell2-1][1])
                elif(mode==3):
                    cost[(cell1,cell2,mode)]=distance(coord[cell1-1][1],coord[cell2-1][0])
                else:
                    cost[(cell1,cell2,mode)]=distance(coord[cell1-1][1],coord[cell2-1][1])
    pathvar=m.addVars(allpath,vtype=grb.GRB.BINARY, name="path") 
    m.setObjective(pathvar.prod(cost), grb.GRB.MINIMIZE)   
    for cell1 in area:
        m.addConstr(
            sum(pathvar.select(cell1,'*',1))+
            sum(pathvar.select(cell1,'*',2))+
            sum(pathvar.select('*',cell1,1))+
            sum(pathvar.select('*',cell1,3))==1,str(cell1)+'a')
        m.addConstr(
            sum(pathvar.select('*',cell1,2))+
            sum(pathvar.select(cell1,'*',3))+
            sum(pathvar.select('*',cell1,4))+
            sum(pathvar.select(cell1,'*',4))==1,str(cell1)+'b')
    m.optimize()
    return (m,pathvar)
'''
bnum=int(input("Enter the number of boundary points:"))
bptlist=[]
for i in range(bnum):
    longitude_in=float(input("Enter boundry longitude:"))
    latitude_in=float(input("Enter boundry latitude:"))
    bpt=(longitude_in, latitude_in)
    bptlist.append(bpt)
bpoints=map(Point,bptlist)
boundary=Polygon(*bpoints)
#check if inputs form a polygon
if(checkcrossline(boundary)):
    bptarray=np.array(bptlist)    
    verticearray=formpolygon(bptarray)
    verticelist=verticearray.tolist()
    boundary=Polygon(*verticelist)
    if not (boundary.is_convex()):
        print "please make sure input points are in correct sequence"
        '''
#make sure vertices are in ccw order
boundary=Polygon((0,0),(0,4),(4,4),(4,0))    
if(boundary.area<0):
    vertices=boundary.vertices
    vertices.reverse()
    boundary=Polygon(*vertices)

#obsnum=int(input("Enter the number of boundary points:"))
obsnum=1
obstacles=[]
allinnerpt=[]
obsleft=[]
obsright=[]
for i in range(obsnum):
    '''
    num=int(input( "Enter the number of No %d obstacle:" %(i+1)))
    left=float('inf')
    right=float('-inf')
    temp=[]
    for j in range(num):
        obsx=float(input("Enter boundry longitude:"))
        obsy=float(input("Enter boundry latitude:"))
        if(obsx<left):
            left=obsx
        if(obsx>right):
            right=obsx
        temp.append([obsx,obsy])
        allinnerpt.append([obsx,obsy])
    temppolygon=Polygon(*temp)
    '''
    temppolygon=Polygon((1,2),(2,3),(3,2),(2,1))
    allinnerpt=np.array([[1,2],[2,3],[3,2],[2,1]])
    left=1
    right=3
    if(temppolygon.area<0):
        vertices=temppolygon.vertices
        vertices.reverse()
        temppolygon=Polygon(*vertices)
    obstacles.append(temppolygon)
    obsleft.append(left)
    obsright.append(right)
allinnerpt=allinnerpt[np.lexsort((allinnerpt[:,1],allinnerpt[:,0]))]
xaxis=Line((0,0),(1,0))
yaxis=Line((0,0),(0,1))
remain=boundary
cells=[]
area={}
coord=[]
for pt in allinnerpt:
    isleft=False
    isright=False
    if(containitem(pt[0],obsleft)):
        isleft=True
    if(containitem(pt[0],obsright)):
        isright=True
    cut=xaxis.perpendicular_line(Point(pt))
    if(isright):
        crosspoint=remain.intersection(cut)
        seg1=Segment(crosspoint[0],crosspoint[1])
        seg2=Segment(crosspoint[1],crosspoint[2])
        cell,remain=splitpolygon(remain,seg1)
        addtomodel(cells,cell,area,coord)
        cells.append(cell)
        cell,remain=splitpolygon(remain,seg2)
        addtomodel(cells,cell,area,coord)
        cells.append(cell)
    if(isleft):
        cell,remain=splitpolygon(remain,cut)
        addtomodel(cells,cell,area,coord)
        cells.append(cell)
        for obs in obstacles:
            if(onpolygon(pt,obs)):
                break
        remain=minuspolygon(remain,obs,pt)
    if not (isleft or isright):
        cross=remain.intersection(cut)
        i=0
        for item in cross:
            if(isinstance(item,Point)):
                i=i+1
            else:
                break
        crosspoint=cross[0:i]
        crossline=cross[i:len(cross)]
        if not(len(crosspoint)==0):
            ptindex=findpoint(crosspoint,Point(pt))
            if((ptindex%2)==0):
                seg=Segment(Point(pt),crosspoint[ptindex+1])
            else:
                seg=Segment(Point(pt),crosspoint[ptindex-1])
            cell,remain=splitpolygon(remain,seg)
            addtomodel(cells,cell,area,coord)
            cells.append(cell)
        if not(len(crossline)==0):
            for line in crossline:
                if(line.contains(Point(pt))):
                    cell,remain=splitpolygon(remain,line)
                    addtomodel(cells,cell,area,coord)
                    cells.append(cell)
                    break
addtomodel(cells,remain,area,coord)
cells.append(remain)

m,pathvar=solveLP(area,coord)
printsolution(m,pathvar)

    