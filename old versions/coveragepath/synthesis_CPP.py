#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:50:19 2019

@author: ysli

features:
    implement TSP to optimize visiting sequence of cells(use lazy constrains for subloop to simplify solving procedure)
    consider centroid of each cell as nodes
    assume all nodes are connected
    consider Euclidean distance as cost of each edge
    for path inside each cell, 2 angles and 4 directions are considered
    the optimal path has minimal traveling time, taking length of the path itself, 
    distance from end point of last cell to first point of this cell, 
    distance from this cell to its centroid and turning cost into consideration
    
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
    if(len(crosspoint)==2 and isinstance(crosspoint[1],Point)):
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
cellcost=[0.25,4/3,4/3,4]
ane=[
     [1,0,0,0],
     [1,1,1,0],
     [0,1,1,1],
     [0,0,0,1]]
kn=[1,1,1,1]
ane,kn,criticalpoints=addtomodel(cells,cellcost,criticalpos)
'''
def addtomodel(cells,cellcost,criticalpos):
    i=0
    for cell in cells:
        bound=cell.bounds
        width=bound[2]-bound[0]
        cellcost.append(width*width/cell.area)
    i=i+1
    ane=np.zeros((len(criticalpos),len(cells)))
    kn=np.zeros((len(criticalpos),1))
    i=0
    for pos in criticalpos:
        j=0
        for cell in cells:       
            bound=cell.bounds
            if(bound[0]==pos or bound[2]==pos):
                ane[i][j]=1
            j=j+1
        i=i+1
    for i in range(len(criticalpos)):
        kn[i]=sum(ane[i])%2
    return ane,kn
def printsolution(m,var):
    if m.status == grb.GRB.Status.OPTIMAL:
        print('totalCost: %g' % m.objVal)
        print('cell copied:')
        for f in var:
            if (var[f].X!=0):
                #print('%s' % (pathvar[f].VarName))
                if(f[0]=='x'):
                    print 'copy cell %d for %d time' %(f[1],var[f].X)
    else:
        print('No solution')
def solveLP(cellcost,ane,kn):
    m=grb.Model('ILP_CPP')
    vars={}
    for i in range(len(cellcost)):
        vars['x',i]= m.addVar(obj=cellcost[i], vtype=grb.GRB.INTEGER,name='x'+str(i))
    for i in range(len(kn)):
        vars['w',i]=m.addVar( vtype=grb.GRB.INTEGER,name='w'+str(i))
    for i in range(len(kn)):
        m.addConstr((grb.quicksum(vars['x',j]*ane[i][j] for j in range(len(cellcost)))-2*vars['w',i])==kn[i])
    m.update()
    m._vars = vars
    m.optimize()
    return m,m._vars
def findbound(vertices):
    leftmost=float('inf')
    rightmost=float('-inf')
    topmost=float('-inf')
    buttommost=float('inf')
    for i in range(len(vertices)):
        vertice=vertices[i]
        if(vertice.x<leftmost):
            leftmost=vertice.x
            leftindex=i
        if(vertice.x>rightmost):
            rightmost=vertice.x
            rightindex=i
        if(vertice.y>topmost):
            topmost=vertice.y
            topindex=i
        if(vertice.y<buttommost):
            buttommost=vertice.y
            buttomindex=i
    return [leftindex,buttomindex,rightindex,topindex]
def formperpenicularlines(xpos):
    xpos=xpos.reshape(len(xpos),1)
    xpts=np.insert(xpos,1,[0],axis=1)
    xpts=map(Point,xpts)
    cutlines=[]
    for pt in xpts:
        cutline=xaxis.perpendicular_line(pt)
        cutlines.append(cutline)
    return cutlines

def createallwaypoint(rboundary,inputwidth):
    bounds=rboundary.bounds
    cutnum=float(math.ceil(round((round(bounds[2],6)-round(bounds[0],6))/inputwidth,6)))
    width=(round(bounds[2],6)-round(bounds[0],6))/cutnum 
    cutpos=np.arange(bounds[0]+width,bounds[2]+width,width)
    pointpos=np.arange(bounds[0]+width/2,bounds[2]+width/2,width)
    cutlines=formperpenicularlines(cutpos)
    #form sub-polygons cut by cutlines
    subpolygons=[]
    remain=rboundary
    for line in cutlines:
        subpolygon,remain=splitpolygon(remain,line)
        subpolygons.append(subpolygon)
        if not (remain):
            break
    #design waypoints in subpolygons
    i=0
    waypoints=[]
    for polygon in subpolygons:
        sides=polygon.sides
        vertices=polygon.vertices
        lines=[]
        lineindex=[]
        j=0
        for side in sides:
            if(side.is_parallel(yaxis)):
                lineindex.append(j)
                lines.append(side)
            j=j+1
        waypoint=[]
        if(len(lines)==2):#formed by 2 cuts
            #seprate vertices into upper part and lower part
            if(lines[0].p1.x<lines[1].p1.x):
                leftupindex=lineindex[0]
                leftdownindex=(lineindex[0]+1)%len(vertices)
                rightupindex=(lineindex[1]+1)%len(vertices)
                rightdownindex=lineindex[1]
            else:
                rightdownindex=lineindex[0]
                rightupindex=(lineindex[0]+1)%len(vertices)
                leftdownindex=(lineindex[1]+1)%len(vertices)
                leftupindex=lineindex[1]
            if(rightupindex<leftupindex):
                uppoints=list(vertices[rightupindex:leftupindex+1])
            else:
                uppoints=list(vertices[rightupindex:len(vertices)])
                uppoints.extend(vertices[0:leftupindex+1])
            if(leftdownindex<rightdownindex):
                downpoints=list(vertices[leftdownindex:rightdownindex+1])
            else:
                downpoints=list(vertices[leftdownindex:len(vertices)])
                downpoints.extend(vertices[0:leftupindex+1]) 
            upindex=findbound(uppoints)
            downindex=findbound(downpoints)
            upmin=uppoints[upindex[1]].y#ymin
            upmax=uppoints[upindex[3]].y #ymax
            downmin=downpoints[downindex[1]].y
            downmax=downpoints[downindex[3]].y
            if(upmax-downmin<=width):
                waypoint.append(Point(pointpos[i],(upmax+downmin)/2))
            else:
                if(upmax-upmin<=width/2):
                    waypoint.append(Point(pointpos[i],upmax-width/2))
                else:
                    cutpt=Point(pointpos[i],upmax-width/2)
                    cut=yaxis.perpendicular_line(cutpt)
                    crosspoint=polygon.intersection(cut)
                    if(crosspoint[0].x<crosspoint[1].x):
                        leftcross=crosspoint[0]
                        rightcross=crosspoint[1]
                    else:
                        leftcross=crosspoint[1]
                        rightcross=crosspoint[0]
                    if(leftcross.x<=pointpos[i] and rightcross.x>=pointpos[i]):
                        waypoint.append(Point(pointpos[i],upmax-width/2))
                    elif(leftcross.x>pointpos[i]):
                        waypoint.append(leftcross)
                        if(upmin-downmax<=width):
                            waypoint.append(Point(pointpos[i],(upmin+downmax)/2))
                        else:
                            waypoint.append(Point(pointpos[i],upmin-width/2))                                
                    else:
                        waypoint.append(rightcross)
                        if(upmin-downmax<=width):
                            waypoint.append(Point(pointpos[i],(upmin+downmax)/2))
                        else:
                            waypoint.append(Point(pointpos[i],upmin-width/2)) 
                if(downmax-downmin<=width/2):
                    waypoint.append(Point(pointpos[i],downmin+width/2))
                else:
                    cutpt=Point(pointpos[i],downmin+width/2)
                    cut=yaxis.perpendicular_line(cutpt)
                    crosspoint=polygon.intersection(cut)
                    if(crosspoint[0].x<crosspoint[1].x):
                        leftcross=crosspoint[0]
                        rightcross=crosspoint[1]
                    else:
                        leftcross=crosspoint[1]
                        rightcross=crosspoint[0]
                    if(leftcross.x<=pointpos[i] and rightcross.x>=pointpos[i]):
                        waypoint.append(Point(pointpos[i],downmin+width/2))
                    elif(leftcross.x>pointpos[i]):
                        if(upmin-downmax<=width):
                            waypoint.append(Point(pointpos[i],(upmin+downmax)/2))
                        else:
                            waypoint.append(Point(pointpos[i],downmax+width/2))
                        waypoint.append(leftcross)
                    else:
                        if(upmin-downmax<=width):
                            waypoint.append(Point(pointpos[i],(upmin+downmax)/2))
                        else:
                            waypoint.append(Point(pointpos[i],downmax+width/2)) 
                        waypoint.append(rightcross)
        elif(len(lines)==1):
            index=findbound(vertices)
            leftindex=index[0]
            line=lines[0]
            ymax=vertices[index[3]].y
            ymin=vertices[index[1]].y
            xmin=vertices[index[0]].x
            xmax=vertices[index[2]].x
            if(line.contains(vertices[leftindex])):#cut line at left
                if(ymax-ymin<=width):
                    waypoint.append(Point(xmax-width/2,(ymin+ymax)/2))
                else:
                    upcutpt=Point(pointpos[i],ymax-width/2)
                    upcut=yaxis.perpendicular_line(upcutpt)
                    downcutpt=Point(pointpos[i],ymin+width/2)
                    downcut=yaxis.perpendicular_line(downcutpt)
                    rightcutpt=Point(xmax-width/2,0)
                    rightcut=xaxis.perpendicular_line(rightcutpt)
                    upcrosspoint=rightcut.intersection(upcut)
                    if(polygon.encloses_point(upcrosspoint[0])):
                        waypoint.append(upcrosspoint[0])
                    else:
                        crosspoint=polygon.intersection(upcut)
                        if(crosspoint[0].x<crosspoint[1].x):
                            waypoint.append(crosspoint[0])
                        else:
                            waypoint.append(crosspoint[1])
                        crosspoint=polygon.intersection(rightcut)
                        if (crosspoint[0].y>crosspoint[1].y):
                            uppt=crosspoint[0]
                            downpt=crosspoint[1]
                        else:
                            uppt=crosspoint[1]
                            downpt=crosspoint[0]
                        if((uppt.y-downpt.y)<=width):
                            waypoint.append(Segment(uppt,downpt).midpoint)
                        else:
                            waypoint.append(Point(uppt.x,uppt.y-width/2))
                    downcrosspoint=rightcut.intersection(downcut)
                    if(polygon.encloses_point(downcrosspoint[0])):
                        waypoint.append(downcrosspoint[0])
                    else:                   
                        crosspoint=polygon.intersection(rightcut)
                        if (crosspoint[0].y>crosspoint[1].y):
                            uppt=crosspoint[0]
                            downpt=crosspoint[1]
                        else:
                            uppt=crosspoint[1]
                            downpt=crosspoint[0]
                        if((uppt.y-downpt.y)<=width):
                            waypoint.append(Segment(uppt,downpt).midpoint)
                        else:
                            waypoint.append(Point(downpt.x,downpt.y+width/2))
                        crosspoint=polygon.intersection(downcut)
                        if(crosspoint[0].x<crosspoint[1].x):
                            waypoint.append(crosspoint[0])
                        else:
                            waypoint.append(crosspoint[1])                                 
            else:#cut line at right
                if(ymax-ymin<=width):
                    waypoint=[Point(xmin+width/2,(ymin+ymax)/2)]
                else:
                    upcutpt=Point(pointpos[i],ymax-width/2)
                    upcut=yaxis.perpendicular_line(upcutpt)
                    downcutpt=Point(pointpos[i],ymin+width/2)
                    downcut=yaxis.perpendicular_line(downcutpt)
                    leftcutpt=Point(xmin+width/2,0)
                    leftcut=xaxis.perpendicular_line(leftcutpt)
                    upcrosspoint=leftcut.intersection(upcut)
                    if(polygon.encloses_point(upcrosspoint[0])):
                        waypoint.append(upcrosspoint[0])
                    else:
                        crosspoint=polygon.intersection(upcut)
                        if(crosspoint[0].x>crosspoint[1].x):
                            waypoint.append(crosspoint[0])
                        else:
                            waypoint.append(crosspoint[1])
                        crosspoint=polygon.intersection(leftcut)
                        if (crosspoint[0].y>crosspoint[1].y):
                            uppt=crosspoint[0]
                            downpt=crosspoint[1]
                        else:
                            uppt=crosspoint[1]
                            downpt=crosspoint[0]
                        if((uppt.y-downpt.y)<=width):
                            waypoint.append(Segment(uppt,downpt).midpoint)
                        else:
                            waypoint.append(Point(uppt.x,uppt.y-width/2))
                    downcrosspoint=leftcut.intersection(downcut)
                    if(polygon.encloses_point(downcrosspoint[0])):
                        waypoint.append(downcrosspoint[0])
                    else:                   
                        crosspoint=polygon.intersection(leftcut)
                        if (crosspoint[0].y>crosspoint[1].y):
                            uppt=crosspoint[0]
                            downpt=crosspoint[1]
                        else:
                            uppt=crosspoint[1]
                            downpt=crosspoint[0]
                        if((uppt.y-downpt.y)<=width):
                            waypoint.append(Segment(uppt,downpt).midpoint)
                        else:
                            waypoint.append(Point(downpt.x,downpt.y+width/2))
                        crosspoint=polygon.intersection(downcut)
                        if(crosspoint[0].x>crosspoint[1].x):
                            waypoint.append(crosspoint[0])
                        else:
                            waypoint.append(crosspoint[1])             
        else:#width of the area less than sensor width
            bound=polygon.bounds
            if(bound[3]-bound[1]<width):
                waypoint=[Point(pointpos[i],(bound[1]+bound[3])/2)]
            else:
                ypos=np.array(bound[3]-width/2,bound[1]+width/2,width)
                ypos=ypos.reshape(len(ypos),1)
                waypoint=ypos.insert(0,pointpos[i],axis=1)
        i=i+1
        if(i%2==1):
            waypoint.reverse()
        waypoints.extend(waypoint)
    return waypoints
def roundpolygon(polygon):
    newvertices=[]
    for vertice in polygon.vertices:
        newvertice=Point(round(vertice.x,6),round(vertice.y,6))
        newvertices.append(newvertice)
    newpolygon=Polygon(*newvertices)
    return newpolygon
def reflectpolygon(polygon):
    rfpolygon=polygon.reflect(xaxis)
    newvertices=rfpolygon.vertices
    newvertices.reverse()
    rfpolygon=Polygon(*newvertices)
    return rfpolygon
    
def singlearearoute(cell,inputwidth,Entrance):
    mintime=float('inf')
    optimalangle=0
    optimalpath=[]
    isoptrf=False
    for angle in range(0,180,180):
        angle=float(angle)
        rcell=cell.rotate(angle/180*math.pi)
        rentrance=Entrance.rotate(angle/180*math.pi)
        rcell=roundpolygon(rcell)
        waypoint=createallwaypoint(rcell,inputwidth)
        time=timeconsume(waypoint,rentrance)
        if(time<mintime):
            mintime=time
            optimalangle=angle
            optimalpath=waypoint
            isoptrf=False
        rvwaypoint=waypoint[:]
        rvwaypoint.reverse()
        time=timeconsume(rvwaypoint,rentrance)
        if(time<mintime):
            mintime=time
            optimalangle=angle
            optimalpath=rvwaypoint
            isoptrf=False   
        rfcell=reflectpolygon(rcell)
        rfentrance=rentrance.reflect(xaxis)
        waypoint=createallwaypoint(rfcell,inputwidth)
        time=timeconsume(waypoint,rfentrance)
        if(time<mintime):
            mintime=time
            optimalangle=angle
            optimalpath=waypoint 
            isoptrf=True
        rvwaypoint=waypoint[:]
        rvwaypoint.reverse()
        time=timeconsume(rvwaypoint,rfentrance)
        if(time<mintime):
            mintime=time
            optimalangle=angle
            optimalpath=rvwaypoint
            isoptrf=True
    for i in range(len(optimalpath)):
        if(isoptrf):
            optimalpath[i]=optimalpath[i].reflect(xaxis)
        optimalpath[i]=optimalpath[i].rotate(-optimalangle/180*math.pi)
    return optimalpath
def timeconsume(waypoint,rentrance):
    length=0
    for i in range(len(waypoint)-1):
        length=length+math.sqrt(pow(waypoint[i].x-waypoint[i+1].x,2)+pow(waypoint[i].y-waypoint[i+1].y,2))
    length=length+waypoint[0].distance(rentrance)
    turnnum=len(waypoint)-1
    return length/velocity+turnnum*turnpanalty
def isvisited(visited,index):
    for i in visited:
        if (index==i):
            return True
    return False
def findnext(newane,edge,current):
    temp=[x[int(edge)] for x in newane]
    nextnode=temp.index(1)
    if(nextnode==current):
        temp.reverse()
        nextnode=len(temp)-temp.index(1)-1
    return nextnode
def pointtolist(vertices):
    pointx=[]
    pointy=[]
    for vertice in vertices:
        pointx.append(vertice.x)
        pointy.append(vertice.y)
    return pointx,pointy                    
def plotpolygon(ax,boundary):    
    vertices=boundary.vertices
    vertices.append(vertices[0])
    verticex,verticey=pointtolist(vertices)
    ax.plot(verticex,verticey)
def plotline(ax,rwaypoints,width):
    pointx,pointy=pointtolist(rwaypoints) 
    ax.plot(pointx,pointy)
    for i in range(len(rwaypoints)-1):
        x=[]
        y=[]
        if(rwaypoints[i].x==rwaypoints[i+1].x and rwaypoints[i].y==rwaypoints[i+1].y):
            continue
        else:
            line=Line(rwaypoints[i],rwaypoints[i+1])
            diracp=line.direction.unit
            diracn=diracp.rotate(math.pi/2)
            rectpt1=rwaypoints[i].translate((diracn.x-diracp.x)*width/2,(diracn.y-diracp.y)*width/2)
            rectpt2=rwaypoints[i+1].translate((diracn.x+diracp.x)*width/2,(diracn.y+diracp.y)*width/2)
            rectpt3=rwaypoints[i+1].translate((diracp.x-diracn.x)*width/2,(diracp.y-diracn.y)*width/2)
            rectpt4=rwaypoints[i].translate((-diracp.x-diracn.x)*width/2,(-diracp.y-diracn.y)*width/2)
            x.extend([rectpt1.x,rectpt2.x,rectpt3.x,rectpt4.x])
            y.extend([rectpt1.y,rectpt2.y,rectpt3.y,rectpt4.y])
            ax.fill(x,y,'b',alpha=0.2)  
def sortvertice(vertices):
    coord=[]
    pts=[]
    for vertice in vertices:
        coord.append([vertice.x,vertice.y]) 
    coord=np.array(coord)
    coord=coord[coord[:,0].argsort()]      
    for pt in coord:
        pts.append(Point(pt[0],pt[1]))
    return pts
def splitCPPpolygon(cell,cutline):
    i=0
    uppoints=[]
    downpoints=[]
    for side in cell.sides:
        if(side.contains(cutline[0])):
            leftupindex=i
            leftdownindex=i+1
        if(side.contains(cutline[-1])):
            rightdownindex=i
            rightupindex=i+1
        i=i+1
    if(leftupindex>rightupindex):
        uppoints.extend(cell.vertices[rightupindex:leftupindex+1])
    else:
        uppoints.extend(cell.vertices[rightupindex:len(cell.vertices)])
        uppoints.extend(cell.vertices[0:leftupindex+1])
    uppoints.extend(cutline)
    if(leftdownindex<rightdownindex):
        downpoints.extend(cell.vertices[leftdownindex:rightdownindex+1])
    else:
        downpoints.extend(cell.vertices[leftdownindex:len(cell.vertices)])
        downpoints.extend(cell.vertices[0:rightdownindex+1])
    cutline.reverse()
    downpoints.extend(cutline)
    downpolygon=Polygon(*downpoints)
    uppolygon=Polygon(*uppoints)
    return downpolygon,uppolygon
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
'''
#sample1
boundary=Polygon((0,0),(4,0),(4,4),(0,4)) 
'''
'''
#sample2
boundary=Polygon((0,0),(9,0),(9,4),(0,4))
'''
'''
#sample3
boundary=Polygon((0,0),(7,0),(7,4),(0,4)) 
'''

#sample4
boundary=Polygon((0,0),(7,0),(7,4),(5,8),(0,4))  

'''
#sample5
boundary=Polygon((0,0),(7,0),(7,6),(0,6))
'''
'''
#sample6
boundary=Polygon((0,0),(9,0),(9,5),(0,5))
if(boundary.area<0):
    vertices=boundary.vertices
    vertices.reverse()
    boundary=Polygon(*vertices)
'''
#obsnum=int(input("Enter the number of boundary points:"))
obsnum=1
obstacles=[]
allinnerpt=[]
obsleft=[]
obsright=[]
'''
#sample1,2
temppolygon=Polygon((1,2),(2,1),(3,2),(2,3))
obstacles.append(temppolygon)
allinnerpt.extend([[1,2],[2,1],[3,2],[2,3]])
obsleft=[1]
obsright=[3]
iscell=[False]
'''
'''
#sample3
temppolygon=Polygon((0.8,2),(0.8,1),(2,1),(2,2))
obstacles.append(temppolygon)
allinnerpt.extend([[0.8,2],[0.8,1],[2,1],[2,2]])
obsleft=[0.8]
obsright=[2]
iscell=[False]
'''

#sample4
temppolygon=Polygon((0.8,2),(0.8,1),(2.3,0.5),(3,2),(2,3))
obstacles.append(temppolygon)
allinnerpt.extend([[0.8,2],[0.8,1],[2.3,0.5],[3,2],[2,3]])
obsleft=[0.8]
obsright=[3]
iscell=[False]

'''
#sample5
temppolygon=Polygon((1,3),(2,1),(4,2),(2,3))
obstacles.append(temppolygon)
allinnerpt.extend([[1,3],[2,1],[4,2],[2,3]])
obsleft=[1]
obsright=[4]
iscell=[True]
'''
'''
#sample6
temppolygon=Polygon((1,3),(2,1),(4,2),(2,3))
obstacles.append(temppolygon)
allinnerpt.extend([[1,3],[2,1],[4,2],[2,3]])
temppolygon=Polygon((5,1.5),(7,1.5),(7,4),(5,4))
obstacles.append(temppolygon)
allinnerpt.extend([[5,1.5],[7,1.5],[7,4],[5,4]])
obsleft=[1,5]
obsright=[4,7]
iscell=[False,True]
'''
'''
for i in range(obsnum):

    num=int(input( "Enter the number of No %d obstacle points:" %(i+1)))
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
    if(temppolygon.area<0):
        vertices=temppolygon.vertices
        vertices.reverse()
        temppolygon=Polygon(*vertices)
    obstacles.append(temppolygon)
    obsleft.append(left)
    obsright.append(right)
'''
allinnerpt=np.array(allinnerpt)
allinnerpt=allinnerpt[np.lexsort((allinnerpt[:,1],allinnerpt[:,0]))]
xaxis=Line((0,0),(1,0))
yaxis=Line((0,0),(0,1))
remain=boundary
cells=[]
cellcost=[]
criticalpos=[]
bounds=remain.bounds
leftbound=bounds[0]
rightbound=bounds[2]
criticalpos.append(leftbound)
for pt in allinnerpt:
    isleft=False
    isright=False
    if(containitem(pt[0],obsleft)):
        isleft=True
    if(containitem(pt[0],obsright)):
        isright=True
    cut=xaxis.perpendicular_line(Point(pt))
    if(isright):
        if not(pt[0]==criticalpos[-1]):
            criticalpos.append(pt[0])
        crosspoint=remain.intersection(cut)
        if(isinstance(crosspoint[-1],Point)):
            seg1=Segment(crosspoint[0],crosspoint[1])
            seg2=Segment(crosspoint[1],crosspoint[2])
            cell,remain=splitpolygon(remain,seg1)            
            cells.append(cell)
            cell,remain=splitpolygon(remain,seg2)
            cells.append(cell)
        else:
            if(len(crosspoint)==1):
                continue
            else:
                seg1=Segment(crosspoint[0],crosspoint[2].p2)
                seg2=Segment(crosspoint[2].p1,crosspoint[1])
                cell,remain=splitpolygon(remain,seg1)
                cells.append(cell)
                cell,remain=splitpolygon(remain,seg2)
                cells.append(cell)
    if(isleft):
        if not(pt[0]==criticalpos[-1]):
            criticalpos.append(pt[0])
        cell,remain=splitpolygon(remain,cut)
        if (remain):
            cells.append(cell)
            for obs in obstacles:
                if(onpolygon(pt,obs)):
                    break
            remain=minuspolygon(remain,obs,pt)
        else:
            remain=cell
cells.append(remain)
for i in range(obsnum):
    if(iscell[i]):
        cells.append(obstacles[i])
criticalpos.append(rightbound)
  
ane,kn=addtomodel(cells,cellcost,criticalpos)

m,var=solveLP(cellcost,ane,kn)
printsolution(m,var)

copiedcell=[]
for f in var:
    if (var[f].X!=0):
        if(f[0]=='x'):
            copiedcell.append(f[1])
 #update ane matrix   
copiedcell.sort() 
copiedcell.reverse()   
newane=ane  
for c in copiedcell:
    newane=np.insert(newane,c,ane[:,c],axis=1)  
copiedcell.reverse()   
criticalpoints=[]  
for pos in criticalpos:
    cutpt=[]
    cut=yaxis.parallel_line(Point(pos,0))
    cross=boundary.intersection(cut)
    if (len(cross)==1):
        if(isinstance(cross[0],Point)):
            criticalpoints.append(cross[0])
            continue
        else:
            criticalpoints.append(cross[0].midpoint)
            continue
    for obstacle in obstacles:
        cross=obstacle.intersection(cut)
        if not(len(cross)==0):
            if(isinstance(cross[0],Point)):
                cutpt.append(cross[0])
            else:
                cutpt.append(cross[0].midpoint)
    xsum=0
    ysum=0
    for pt in cutpt:
        xsum=xsum+pt.x
        ysum=ysum+pt.y
    criticalpoints.append(Point(xsum/len(cutpt),ysum/len(cutpt)))
j=0

for c in copiedcell:
    cutpt=[]   
    if(cells[c+j].is_convex()):
        cutline=[]
        for i in range(len(ane[:,c])):
            if not(ane[i,c]==0):            
                cutpt.append(criticalpoints[i])
        cutline.append(cutpt[0])
        for side in cells[c+j].sides:
            if(side.contains(cutpt[0])):
                if(side.p1.y>side.p2.y):
                    leftside=Segment(side.p1,side.p2)
                else:
                    leftside=Segment(side.p2,side.p1)
            if(side.contains(cutpt[1])):
                if(side.p1.y>side.p2.y):
                    rightside=Segment(side.p1,side.p2)
                else:
                    rightside=Segment(side.p2,side.p1)
        leftpartion=(leftside.p1.y-cutpt[0].y)/(leftside.p1.y-leftside.p2.y)
        rightpartion=(rightside.p1.y-cutpt[1].y)/(rightside.p1.y-rightside.p2.y)
        vertices=sortvertice(cells[c+j].vertices)
        for vertice in vertices:
            cut=yaxis.parallel_line(vertice)
            cross=cells[c+j].intersection(cut)
            if(isinstance(cross[-1],Point)):
                if(cross[0].y>cross[1].y):
                    upcross=cross[0]
                    downcross=cross[1]
                else:
                    upcross=cross[1]
                    downcross=cross[0]
                partion=(rightpartion-leftpartion)*(vertice.x-leftside.p1.x)/(rightside.p1.x-leftside.p1.x)+leftpartion
                cut=Point(vertice.x,upcross.y-(upcross.y-downcross.y)*partion)
                cutline.append(cut)
        cutline.append(cutpt[1])
        cell1,cell2=splitCPPpolygon(cells[c+j],cutline)
    else:
        cutpt=[]  
        bound=cells[c+j].bounds
        leftcut=yaxis.parallel_line(Point(bound[0],0))
        rightcut=yaxis.parallel_line(Point(bound[2],0))
        leftcross=cells[c+j].intersection(leftcut)
        rightcross=cells[c+j].intersection(rightcut)
        if(isinstance(leftcross,Point)):
            cutpt.append(leftcross[0])
        else:
            cutpt.append(leftcross[0].midpoint)
        if(isinstance(rightcross,Point)):
            cutpt.append(rightcross[0])
        else:
            cutpt.append(rightcross[0].midpoint)  
        cutseg=Segment(cutpt[0],cutpt[1])
        cell1,cell2=splitpolygon(cells[c+j],cutseg)
    cells.remove(cells[c+j])
    cells.insert(c+j,cell1)
    cells.insert(c+j,cell2)
    j=j+1
newane=newane.tolist()  
waypoints=[]  
velocity=1
turnpanalty=10
inputwidth=0.4
visitededge=[0]
visitednode=[0]
currentnode=0
cellwaypoint=singlearearoute(cells[0],inputwidth,criticalpoints[0])
waypoints.extend(cellwaypoint)
while(len(visitededge)<len(cells)):
    currentnode=findnext(newane,visitededge[-1],currentnode)
    if(len(visitednode)<len(criticalpoints)-1):      
        edgeoptions=[]
        for i in range(len(cells)):
            if(newane[currentnode][i]==1):
                edge=i
                if not(findnext(newane,edge,currentnode) in visitednode):
                    edgeoptions.append(edge)
        visitednode.append(currentnode)
    else:
        edgeoptions=[]
        nodedegree=[0]*len(criticalpoints)
        edgearray=[]
        for i in range(len(criticalpoints)):
            edgearray.append([])
        for i in range(len(cells)):
            edge=0
            if(newane[currentnode][i]==1):
                edge=i
            if not(edge in visitededge):
                nextnode=findnext(newane,edge,currentnode)
                nodedegree[nextnode]+=1
                edgearray[nextnode].append(edge)
        nextnode=nodedegree.index(max(nodedegree))    
        edgeoptions.extend(edgearray[nextnode])
    mindistance=[]
    for cell in edgeoptions:
        mindist=float('inf')
        for vertice in cells[int(cell)].vertices:           
            distance=Point(waypoints[-1]).distance(vertice)
            if(distance<mindist):
                mindist=distance
        mindistance.append(mindist)
    nextedge=edgeoptions[mindistance.index(min(mindistance))]
    if(nextedge>=len(cells)-2):
        width=0.4
    else:
        width=inputwidth
    cellwaypoint=singlearearoute(cells[nextedge],width,waypoints[-1])
    waypoints.extend(cellwaypoint)
    visitededge.append(nextedge)
    
    
    
fig=plt.figure()
ax = fig.add_subplot(1,1,1)
for cell in cells:
    plotpolygon(ax,cell)
plotline(ax,waypoints,inputwidth)

