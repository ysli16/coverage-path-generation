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
    bound=findbound(newcell.vertices)
    leftbound=newcell.vertices[bound[0]]
    rightbound=newcell.vertices[bound[2]]
    leftcut=yaxis.parallel_line(leftbound)
    rightcut=yaxis.parallel_line(rightbound)
    leftcross=newcell.intersection(leftcut)[0]
    rightcross=newcell.intersection(rightcut)[0]
    if(isinstance(leftcross,Point)):
        a=leftcross
    else:
        a=leftcross.midpoint
    if(isinstance(rightcross,Point)):
        b=rightcross
    else:
        b=rightcross.midpoint    
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

def createwaypoint(up1,down1,pos,width):
    upbound=up1.bounds
    upmin=upbound[1]#ymin
    upmax=upbound[3]#ymax        
    downbound=down1.bounds
    downmin=downbound[1]#ymin
    downmax=downbound[3]#ymax
    #to cover the rectangle part
    if(upmin-downmax<=width):
        waypoint=[Point(pos,(upmin+downmax)/2)]
    else:
        waypoint=[Point(pos,upmin-width/2),Point(pos,downmax+width/2)]
    #to cover upper irragular shape
    if(upmax-upmin<=width/2):
        del waypoint[0]
        waypoint.insert(0,Point(pos,upmax-width/2))
    else:
        cutpt=Point(pos,upmax-width/2)
        cut=yaxis.perpendicular_line(cutpt)
        crosspoint=up1.intersection(cut)
        if(crosspoint[0].x<crosspoint[1].x):
            leftcross=crosspoint[0]
            rightcross=crosspoint[1]
        else:
            leftcross=crosspoint[1]
            rightcross=crosspoint[0]
        if(leftcross.x<=pos and rightcross.x>=pos):
            waypoint.insert(0,Point(pos,upmax-width/2))
        elif(leftcross.x>pos):
            waypoint.insert(0,leftcross)
        else:
            waypoint.insert(0,rightcross)
    #to cover lower irragular shape
    if(downmax-downmin<=width/2):
        waypoint.append(Point(pos,downmin+width/2))
    else:
        cutpt=Point(pos,downmin+width/2)
        cut=yaxis.perpendicular_line(cutpt)
        crosspoint=down1.intersection(cut)
        if(crosspoint[0].x<crosspoint[1].x):
            leftcross=crosspoint[0]
            rightcross=crosspoint[1]
        else:
            leftcross=crosspoint[1]
            rightcross=crosspoint[0]
        if(leftcross.x<=pos and rightcross.x>=pos):
            waypoint.append(Point(pos,downmin+width/2))
        elif(leftcross.x>pos):
            waypoint.append(leftcross)
        else:
            waypoint.append(rightcross)
    return waypoint
def createallwaypoint(rboundary,inputwidth):    
    
    #define the width,generate cut position
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
        for side in sides:
            if(side.is_parallel(yaxis)):
                lines.append(side)
        if(len(lines)==2):#formed by 2 cuts
            if(lines[0].p1.x<lines[1].p1.x):
                leftline=lines[0]
                rightline=lines[1]
            else:
                leftline=lines[1]
                rightline=lines[0]
            if(leftline.p1.y<leftline.p2.y):
                leftline=Line(leftline.p2,leftline.p1)
            if(rightline.p1.y<rightline.p2.y):
                rightline=Line(rightline.p2,rightline.p1)
            #seprate vertices into upper part and lower part
            leftupindex=findpoint(vertices,leftline.p1)
            leftdownindex=findpoint(vertices,leftline.p2)
            rightupindex=findpoint(vertices,rightline.p1)
            rightdownindex=findpoint(vertices,rightline.p2)  
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
            index=findbound(uppoints)
            buttomindex=index[1]#upper min y
            if(index[1]==index[3]):#upper part is segment
                up1=None
            else:
                cut=xaxis.parallel_line(uppoints[buttomindex])
                up1,temp=splitpolygon(polygon,cut)
            index=findbound(downpoints)
            topindex=index[3]#lower max y
            if(index[1]==index[3]):#lower part is segment
                down1=None
            else:
                cut=xaxis.parallel_line(downpoints[topindex])
                temp,down1=splitpolygon(polygon,cut)
            waypoint=[]
            if (up1 and down1):
                waypoint=createwaypoint(up1,down1,pointpos[i],width)
            elif(up1):#lower part is segment
                upbound=up1.bounds
                upmin=upbound[1]#ymin
                upmax=upbound[3] #ymax
                buttom=downpoints[topindex].y#y of lower part
                if(upmax-upmin<=width/2):
                    waypoint=[Point(pointpos[i],upmax-width/2)]
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
                        waypoint=[Point(pointpos[i],upmax-width/2)]
                    elif(leftcross.x>pointpos[i]):
                        waypoint=[leftcross]
                    else:
                        waypoint=[rightcross]
                if(upmin-buttom>width):
                    waypoint.append(Point(pointpos[i],upmin-width/2))
                waypoint.append(Point(pointpos[i],buttom+width/2))
            elif(down1):
                downbound=down1.bounds
                downmin=downbound[1]#ymin
                downmax=downbound[3]#ymax
                top=uppoints[buttomindex].y#upper segment y
                if(downmax-downmin<=width/2):
                    waypoint=[Point(pointpos[i],downmin+width/2)]
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
                        waypoint=[Point(pointpos[i],downmin+width/2)]
                    elif(leftcross.x>pointpos[i]):
                        waypoint=[leftcross]
                    else:
                        waypoint=[rightcross]
                if(top-downmax>=width):
                    waypoint.append(Point(pointpos[i],downmax+width/2))
                waypoint.append(Point(pointpos[i],top-width/2))
                waypoint.reverse()
            else:
                buttom=downpoints[topindex].y
                top=uppoints[buttomindex].y
                if(top-buttom<=width):
                    waypoint=[Point(pointpos[i],(top+buttom)/2)]
                else:
                    waypoint=[Point(pointpos[i],top-width/2),Point(pointpos[i],buttom+width/2)]
        elif(len(lines)==1):
            index=findbound(vertices)
            leftindex=index[0]
            rightindex=index[2]
            line=lines[0]
            if(line.p1.y<line.p2.y):
                line=Line(line.p2,line.p1)
            bound=polygon.bounds
            ymax=bound[3]
            ymin=bound[1]
            xmin=bound[0]
            xmax=bound[2]
            if(line.contains(vertices[leftindex])):#cut line at left
                if(ymax-ymin<=width):
                    waypoint=[Point(xmax-width/2,(ymin+ymax)/2)]
                else:
                    cutpt=Point(pointpos[i],ymax-width/2)
                    cut=yaxis.perpendicular_line(cutpt)
                    crosspoint=polygon.intersection(cut)
                    if(crosspoint[0].x<crosspoint[1].x):
                        crosspoint=crosspoint[1]
                    else:
                        crosspoint=crosspoint[0]
                    if(crosspoint.x>pointpos[i]):
                        waypoint=[Point(pointpos[i],ymax-width/2)]
                    else:
                        waypoint=[crosspoint]
                    cutpt=Point(pointpos[i],ymin+width/2)
                    cut=yaxis.perpendicular_line(cutpt)
                    crosspoint=polygon.intersection(cut)
                    if(crosspoint[0].x<crosspoint[1].x):
                        crosspoint=crosspoint[1]
                    else:
                        crosspoint=crosspoint[0]
                    if(crosspoint.x>pointpos[i]):
                        waypoint.append(Point(pointpos[i],ymin+width/2))
                    else:
                        waypoint.append(crosspoint)
                    #see if segment consist of 2 waypoints can cover rightmost part
                    way=Segment(*waypoint)
                    wp=Point(vertices[rightindex].x-width/2,vertices[rightindex].y)
                    limit=xaxis.perpendicular_line(wp)
                    if not(way.intersection(limit)):
                        waypoint.insert(1,wp)
            else:#cut line at right
                if(ymax-ymin<=width):
                    waypoint=[Point(xmin+width/2,(ymax+ymin)/2)]
                else:
                    cutpt=Point(pointpos[i],ymax-width/2)
                    cut=yaxis.perpendicular_line(cutpt)
                    crosspoint=polygon.intersection(cut)
                    if(crosspoint[0].x<crosspoint[1].x):
                        crosspoint=crosspoint[0]
                    else:
                        crosspoint=crosspoint[1]
                    if(crosspoint.x<pointpos[i]):
                        waypoint=[Point(pointpos[i],ymax-width/2)]
                    else:
                        waypoint=[crosspoint]
                    cutpt=Point(pointpos[i],ymin+width/2)
                    cut=yaxis.perpendicular_line(cutpt)
                    crosspoint=polygon.intersection(cut)
                    if(crosspoint[0].x<crosspoint[1].x):
                        crosspoint=crosspoint[0]
                    else:
                        crosspoint=crosspoint[1]
                    if(crosspoint.x<pointpos[i]):
                        waypoint.append(Point(pointpos[i],ymin+width/2))
                    else:
                        waypoint.append(crosspoint)
                    way=Segment(*waypoint)
                    wp=Point(vertices[leftindex].x+width/2,vertices[leftindex].y)
                    limit=xaxis.perpendicular_line(wp)
                    if not(way.intersection(limit)):
                        waypoint.insert(1,wp)           
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

def singlearearoute(cell,inputwidth):
    mintime=float('inf')
    optimalangle=0
    for angle in range(0,180,90):
        angle=float(angle)
        rcell=cell.rotate(angle/180*math.pi)
        waypoint=createallwaypoint(rcell,inputwidth)
        Entrance=cell.sides[0]
        rentrance=Entrance.rotate(angle/180*math.pi)
        Exit=cell.sides[-1]
        rexit=Exit.rotate(angle/180*math.pi)
        time=timeconsume(waypoint,rentrance,rexit)
        if(time<mintime):
            mintime=time
            optimalangle=angle
    rcell=cell.rotate(optimalangle/180*math.pi) 
    waypoints=createallwaypoint(rcell,inputwidth)   
    rwaypoints=[]
    for i in range(len(waypoints)):
        rwaypoints.append(waypoints[i].rotate(-optimalangle/180*math.pi))
    return rwaypoints
def timeconsume(waypoint,rentrance,rexit):
    length=0
    for i in range(len(waypoint)-1):
        length=length+math.sqrt(pow(waypoint[i].x-waypoint[i+1].x,2)+pow(waypoint[i].y-waypoint[i+1].y,2))
    length=length+waypoint[0].distance(rentrance)
    length=length+waypoint[len(waypoint)-1].distance(rexit)
    turnnum=len(waypoint)-1
    return length/velocity+turnnum*turnpanalty
def isvisited(visited,index):
    for i in visited:
        if (index==i):
            return True
    return False
def findnext(pathvar,visited,current):
    for p in pathvar:
        if not(pathvar[p].X==0):
            if(p[0]==current):
                if(len(visited)==1):
                    if (p[2]==3 or p[2]==4):
                        return p
                else:
                    if not(isvisited(visited,p[1])):
                        return p
            elif(p[1]==current):
                if not(isvisited(visited,p[0])):
                    return p
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
boundary=Polygon((0,0),(0,4),(9,4),(9,0))    
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
    
    temppolygon=Polygon((1,2),(2,3),(3,2),(2,1))
    allinnerpt=np.array([[1,2],[2,3],[3,2],[2,1]])
    '''
    temppolygon=Polygon((1,1),(2,1),(2,2),(1,2))
    allinnerpt=np.array([[1,1],[2,1],[2,2],[1,2]])
    left=1
    right=2
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
        if(isinstance(crosspoint[-1],Point)):
            seg1=Segment(crosspoint[0],crosspoint[1])
            seg2=Segment(crosspoint[1],crosspoint[2])
            cell,remain=splitpolygon(remain,seg1)
            addtomodel(cells,cell,area,coord)
            cells.append(cell)
            cell,remain=splitpolygon(remain,seg2)
            addtomodel(cells,cell,area,coord)
            cells.append(cell)
        else:
            if(len(crosspoint)==1):
                continue
            else:
                seg1=Segment(crosspoint[0],crosspoint[2].p2)
                seg2=Segment(crosspoint[2].p1,crosspoint[1])
                cell,remain=splitpolygon(remain,seg1)
                addtomodel(cells,cell,area,coord)
                cells.append(cell)
                cell,remain=splitpolygon(remain,seg2)
                addtomodel(cells,cell,area,coord)
                cells.append(cell)
    if(isleft):
        cell,remain=splitpolygon(remain,cut)
        if (remain):
            addtomodel(cells,cell,area,coord)
            cells.append(cell)
            for obs in obstacles:
                if(onpolygon(pt,obs)):
                    break
            remain=minuspolygon(remain,obs,pt)
        else:
            remain=cell
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

velocity=1
turnpanalty=10
inputwidth=0.4
fig=plt.figure()
ax = fig.add_subplot(1,1,1)
waypoints=[]
current=1
visited=[current]
cellwaypoint=singlearearoute(cells[current-1],inputwidth)
plotpolygon(ax,cells[current-1])
waypoints.extend(cellwaypoint)
while(len(visited)<len(cells)):      
    nextvar=findnext(pathvar,visited,current)
    if(nextvar[0]==current):
        current=nextvar[1]
        cellwaypoint=singlearearoute(cells[current-1],inputwidth)
        if(nextvar[2]==2 or nextvar[2]==4):
            cellwaypoint.reverse()
    else:
        current=nextvar[0]
        cellwaypoint=singlearearoute(cells[current-1],inputwidth)
        if(nextvar[2]==3 or nextvar[2]==4):
            cellwaypoint.reverse()
    waypoints.extend(cellwaypoint)
    visited.append(current)
    plotpolygon(ax,cells[current-1])
plotline(ax,waypoints,inputwidth)
