#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:50:19 2019

@author: ysli

features:
    implement TSP to optimize visiting sequence of cells(use lazy constrains for subloop to simplify solving procedure)
    consider entrance and exit(defined as midpoint) of each cell as nodes
    assume all nodes are connected
    consider Euclidean distance as cost of each edge
    constrain each edge inside cell(i.e. edge connecting entrance and exit) at least once
    for path inside each cell, 2 angles and 2 directions are considered,also consider reverse the path when horizental
    the optimal path has minimal traveling time, taking length of the path itself, 
    distance from end point of last cell to first point of this cell, 
    distance from this cell to its exit and turning cost
    
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
points=[(0,2),(1,2),(1,1),(2,0.5),(1,3),(2,3.5),(2,0.5),(3,1),(2,3.5),(3,3),(3,2),(4,2)]
num=12
'''
def addtomodel(newcell,points):
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
    points.extend([(a.x,a.y),(b.x,b.y)])
def distance(point1,point2):
    distance=math.sqrt(pow(point1[0]-point2[0],2)+pow(point1[1]-point2[1],2))
    return distance
# Callback - use lazy constraints to eliminate sub-tours
def subtourelim(model, where):
  if where == grb.GRB.callback.MIPSOL:
    selected = []
    # make a list of edges selected in the solution
    for i in range(num):
      sol = model.cbGetSolution([model._vars[i,j] for j in range(num)])
      selected += [(i,j) for j in range(num) if sol[j] > 0.5]
    # find the shortest cycle in the selected edge list
    tour = subtour(selected)
    if len(tour) < num:
      # add a subtour elimination constraint
      expr = 0
      for i in range(len(tour)):
        for j in range(i+1, len(tour)):
          expr += model._vars[tour[i], tour[j]]
      model.cbLazy(expr <= len(tour)-1)
# Given a list of edges, finds the shortest subtour     
def subtour(edges):
  visited = [False]*num
  cycles = []
  lengths = []
  selected = [[] for i in range(num)]
  for x,y in edges:
    selected[x].append(y)
  while True:
    current = visited.index(False)
    thiscycle = [current]
    while True:
      visited[current] = True
      neighbors = [x for x in selected[current] if not visited[x]]
      if len(neighbors) == 0:
        break
      current = neighbors[0]
      thiscycle.append(current)
    cycles.append(thiscycle)
    lengths.append(len(thiscycle))
    if sum(lengths) == num:
      break
  return cycles[lengths.index(min(lengths))]
def printsolution(m,pathvar):
    if m.status == grb.GRB.Status.OPTIMAL:
        print('totalCost: %g' % m.objVal)
        print('Path:')
        for f in pathvar:
            if(round(pathvar[f].X,5)==1):
                if(f[0]<f[1]):
                    no1=f[0]/2
                    no2=f[1]/2
                    if(no1==no2):
                        print('explore cell %d' %(no1))
                    else:
                        if(f[0]%2==0 and f[1]%2==0):
                            print('from cell %d a to cell %d a' %((no1),(no2)))
                        elif(f[0]%2==0 and f[1]%2==1):
                            print('from cell %d a to cell %d b' %((no1),(no2)))
                        elif(f[0]%2==1 and f[1]%2==0):
                            print('from cell %d b to cell %d a' %((no1),(no2)))
                        else:
                            print('from cell %d b to cell %d b' %((no1),(no2)))
    else:
        print('No feasible solution')
def solveLP(points):
    # Create variables
    m = grb.Model()
    vars = {}
    for i in range(num):
       for j in range(i+1):
           vars[i,j] = m.addVar(obj=distance(points[i],points[j]), vtype=grb.GRB.BINARY,
                              name='e'+str(i)+'_'+str(j))
           vars[j,i] = vars[i,j]
       m.update()
    # Add degree-2 constraint, and forbid loops   
    for i in range(num):
      m.addConstr(grb.quicksum(vars[i,j] for j in range(num)) == 2)
      vars[i,i].ub = 0
    #all egdes inside each cell must be visited at least once
    for i in range(0,num,2):
        vars[i,i+1].lb=1
    m.update()
    m._vars = vars
    m.params.LazyConstraints = 1
    m.optimize(subtourelim)
    return m
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
                waypoint=[Point(pointpos[i],bound[3]-width/2),Point(pointpos[i],bound[1]+width/2)]
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
    newvertices=list(rfpolygon.vertices)
    newvertices.reverse()
    rfpolygon=Polygon(*newvertices)
    return rfpolygon
    
def singlearearoute(cell,inputwidth,Entrance,Exit):
    mintime=float('inf')
    optimalangle=0
    optimalpath=[]
    isoptrf=False
    Entrance=Point(Entrance)
    Exit=Point(Exit)
    for angle in range(0,180,90):
        angle=float(angle)
        rcell=cell.rotate(angle/180*math.pi)
        rentrance=Entrance.rotate(angle/180*math.pi)
        rexit=Exit.rotate(angle/180*math.pi)
        rcell=roundpolygon(rcell)
        waypoint=createallwaypoint(rcell,inputwidth)
        time=timeconsume(waypoint,rentrance,rexit)
        if(time<mintime):
            mintime=time
            optimalangle=angle
            optimalpath=waypoint
            isoptrf=False
        rvwaypoint=waypoint[:]
        rvwaypoint.reverse()
        time=timeconsume(rvwaypoint,rentrance,rexit)
        if(time<mintime):
            mintime=time
            optimalangle=angle
            optimalpath=rvwaypoint
            isoptrf=False   
        rfcell=reflectpolygon(rcell)
        rfentrance=rentrance.reflect(xaxis)
        rfexit=rexit.reflect(xaxis)
        waypoint=createallwaypoint(rfcell,inputwidth)
        time=timeconsume(waypoint,rfentrance,rfexit)
        if(time<mintime):
            mintime=time
            optimalangle=angle
            optimalpath=waypoint 
            isoptrf=True
        rvwaypoint=waypoint[:]
        rvwaypoint.reverse()
        time=timeconsume(rvwaypoint,rfentrance,rfexit)
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
def timeconsume(waypoint,rentrance,rexit):
    length=0
    for i in range(len(waypoint)-1):
        length=length+math.sqrt(pow(waypoint[i].x-waypoint[i+1].x,2)+pow(waypoint[i].y-waypoint[i+1].y,2))
    length=length+waypoint[0].distance(rentrance)
    length=length+waypoint[-1].distance(rexit)
    turnnum=len(waypoint)-1
    return length/velocity+turnnum*turnpanalty
def isvisited(visited,index):
    for i in visited:
        if (index==i):
            return True
    return False
def findnext(pathvar,visited,current):
    for p in pathvar:
        if (round(pathvar[p].X,5)==1):
            if(p[0]<p[1]):
                no1=p[0]/2
                no2=p[1]/2
                if(no1==current):
                    if(len(visited)==1):
                        if (p[0]%2==1):
                            return p
                    else:
                        if not(isvisited(visited,no2)):
                            return p
                elif(no2==current):
                    if not(isvisited(visited,no1)):
                        return p
def pointtoarray(vertices):
    point=[]
    for vertice in vertices:
        point.append([vertice.x,vertice.y])
    return np.array(point)                 
def plotpolygon(ax,boundary):    
    vertices=list(boundary.vertices)
    vertices.append(vertices[0])
    vertice=pointtoarray(vertices)
    ax.plot(vertice[:,0],vertice[:,1])
def plotline(ax,rwaypoints,width):
    point=pointtoarray(rwaypoints) 
    ax.plot(point[:,0],point[:,1])
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
boundary=Polygon((0,0),(0,4),(7,4),(7,0)) 
'''
'''
#sample4
boundary=Polygon((0,0),(7,0),(7,4),(5,8),(0,4))  
'''
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
ratio=10
#lake Mascoma
bpt=[(0,0.4),(0.55,0.1),(2,0.1),(2.93,0.65),(3.75,0.6),
        (4.05,0.3),(4.55,0.4),(5.25,0.1),(7.1,0.5),(9.3,0.2),
        (10.2,0.75),(11.1,0.75),(11.5,0.4),(13.7,0),(14.5,0.5),
        (15.4,0.1),(16,0.8),(16,1.9),(14.3,1.7),(13.65,2.2),
        (13.3,2.1),(11.4,2.4),(10.7,2),(10.3,2.8),(8.2,3.2),
        (7.8,2.95),(6.5,2.95),(6.75,1.3),(5.3,0.8),(4,1.95),
        (3.2,2.1),(1,0.8),(0.3,0.7)]
bpt=np.array(bpt)
bpt=bpt*ratio
boundary=map(Point,bpt)
boundary=Polygon(*boundary)
#obsnum=int(input("Enter the number of boundary points:"))
obsnum=4
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
'''
#sample4
temppolygon=Polygon((0.8,2),(0.8,1),(2.3,0.5),(3,2),(2,3))
obstacles.append(temppolygon)
allinnerpt.extend([[0.8,2],[0.8,1],[2.3,0.5],[3,2],[2,3]])
obsleft=[0.8]
obsright=[3]
iscell=[False]
'''
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
#lake Mascoma
temppt=[(8.9,2.5),(9.1,2.4),(9.25,2.5),(9,2.7)]
temppt=np.array(temppt)
temppt=temppt*ratio
temppolygon=map(Point,temppt)
temppolygon=Polygon(*temppolygon)
obstacles.append(temppolygon)
inner=np.array([[8.9,2.5],[9.1,2.4],[9.25,2.5],[9,2.7]])
inner=inner*ratio
allinnerpt.extend(inner)

temppt=[(9.3,2.55),(9.4,2.4),(9.6,2.4),(9.6,2.5),(9.5,2.6)]
temppt=np.array(temppt)
temppt=temppt*ratio
temppolygon=map(Point,temppt)
temppolygon=Polygon(*temppolygon)
obstacles.append(temppolygon)
inner=np.array([[9.3,2.55],[9.4,2.4],[9.6,2.4],[9.6,2.5],[9.5,2.6]])
inner=inner*ratio
allinnerpt.extend(inner)

temppt=[(7,2.4),(7.15,0.9),(8.45,0.5),(9.2,0.5),(10.15,0.9),(10.15,1.7),(8.7,1.9),(7.25,2.6)]
temppt=np.array(temppt)
temppt=temppt*ratio
temppolygon=map(Point,temppt)
temppolygon=Polygon(*temppolygon)
obstacles.append(temppolygon)
inner=np.array([[7,2.4],[7.15,0.9],[8.45,0.5],[9.2,0.5],[10.15,0.9],[10.15,1.7],[8.7,1.9],[7.25,2.6]])
inner=inner*ratio
allinnerpt.extend(inner)

temppt=[(11.75,1.5),(12.4,1),(12.75,1.25),(12.45,2),(11.9,1.85)]
temppt=np.array(temppt)
temppt=temppt*ratio
temppolygon=map(Point,temppt)
temppolygon=Polygon(*temppolygon)
obstacles.append(temppolygon)
inner=np.array([[11.75,1.5],[12.4,1],[12.75,1.25],[12.45,2],[11.9,1.85]])
inner=inner*ratio
allinnerpt.extend(inner)
obsleft=np.array([8.9,9.3,7,11.75])*ratio
obsright=np.array([9.25,9.6,10.15,12.75])*ratio
iscell=[False,False,True,True]
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
points=[]
fig=plt.figure()
ax = fig.add_subplot(1,1,1)
for pt in allinnerpt:
    isleft=False
    isright=False
    if(containitem(pt[0],obsleft)):
        isleft=True
    if(containitem(pt[0],obsright)):
        isright=True
    cut=xaxis.perpendicular_line(Point(pt))
    cross=remain.intersection(cut)
    crosspt=[]
    for item in cross:
        if(isinstance(item,Point)):
            crosspt.append(item)
        else:
            i=0
            for cpt in crosspt:
                if cpt.y<item.p1.y:
                    i=i+1
                    continue
                else:                    
                    crosspt.insert(i,item.p1)
                    crosspt.insert(i,item.p2)
                    break
    if(isright):
        if not (Point(pt) in crosspt):
            continue
        elif not(Point(pt) in cross):
            index=crosspt.index(Point(pt))
            seg1=Segment(crosspt[index-1],crosspt[index])
            seg2=Segment(crosspt[index+1],crosspt[index+2])
            cell,remain=splitpolygon(remain,seg1)
            cells.append(cell)
            plotpolygon(ax,cell)
            cell,remain=splitpolygon(remain,seg2)
            cells.append(cell)
            plotpolygon(ax,cell)
        else:
            index=crosspt.index(Point(pt))
            seg1=Segment(crosspt[index-1],crosspt[index])
            seg2=Segment(crosspt[index],crosspt[index+1])
            cell,remain=splitpolygon(remain,seg1)
            cells.append(cell)
            plotpolygon(ax,cell)
            cell,remain=splitpolygon(remain,seg2)
            cells.append(cell)  
            plotpolygon(ax,cell)
    if(isleft):
        if (Point(pt) in crosspt):
            continue
        else:
            index=0
            for cpt in crosspt:
                if(cpt.y>pt[1]):
                   break 
                index=index+1
            seg=Segment(crosspt[index-1],crosspt[index])
            cell,remain=splitpolygon(remain,seg)
            if (remain):
                cells.append(cell)
                plotpolygon(ax,cell)
            else:
                remain=cell
            for obs in obstacles:
                if(onpolygon(pt,obs)):
                    break
            remain=minuspolygon(remain,obs,pt)
    if not (isleft or isright):
        for obs in obstacles:
            if(onpolygon(pt,obs)):
                break
        if (obs.angles[Point(pt)]<math.pi):
            index=findpoint(crosspt,Point(pt))
            if(index%2==0):
                seg=Segment(crosspt[index],crosspt[index+1])
            else:
                seg=Segment(crosspt[index-1],crosspt[index])
            cell,remain=splitpolygon(remain,seg)
            cells.append(cell)   
            plotpolygon(ax,cell)

cells.append(remain)
plotpolygon(ax,remain)
width=np.ones(len(cells))*0.5
for i in range(obsnum):
    if(iscell[i]):
        cells.append(obstacles[i])
        plotpolygon(ax,obstacles[i])
        width=np.append(width,2.4)

cindex=0
count=0 
convexcells=list(cells) 
for cell in cells:
    cutpt=[]
    if not(cell.is_convex()):
        for vertice in cell.vertices:
            if(cell.angles[vertice]>math.pi):
                cutpt.append(vertice)
        pts=pointtoarray(cutpt)
        pts = pts[pts[:,0].argsort()]
        for pt in pts:
            cut=yaxis.parallel_line(Point(pt))
            cross=convexcells[cindex+count].intersection(cut)
            if(len(cross)==2):
                cell1,cell2=splitpolygon(convexcells[cindex+count],cut)    
                convexcells.remove(convexcells[cindex+count])
                convexcells.insert(cindex+count,cell2)
                convexcells.insert(cindex+count,cell1) 
                width=np.insert(width,cindex+count,width[cindex+count])
                plotpolygon(ax,cell1)
                count=count+1
            elif(len(cross)==3):
                seg1=Segment(cross[0],cross[1])
                cell1,cell2=splitpolygon(convexcells[cindex+count],seg1)    
                convexcells.remove(convexcells[cindex+count])
                convexcells.insert(cindex+count,cell2)
                convexcells.insert(cindex+count,cell1) 
                width=np.insert(width,cindex+count,width[cindex+count])
                plotpolygon(ax,cell1)
                count=count+1
                seg2=Segment(cross[1],cross[2])
                cell1,cell2=splitpolygon(convexcells[cindex+count],seg2)    
                convexcells.remove(convexcells[cindex+count])
                convexcells.insert(cindex+count,cell2)
                convexcells.insert(cindex+count,cell1)  
                width=np.insert(width,cindex+count,width[cindex+count])
                plotpolygon(ax,cell1)
                count=count+1
    cindex=cindex+1
for cell in convexcells:
    addtomodel(cell,points)


num=len(points)
m=solveLP(points)
printsolution(m,m._vars)

velocity=1
turnpanalty=10

waypoints=[]
current=0
visited=[current]
cellwaypoint=singlearearoute(convexcells[current],width[current],points[(current)*2],points[(current)*2+1])
waypoints.extend(cellwaypoint)
while(len(visited)<len(convexcells)):      
    nextvar=findnext(m._vars,visited,current)
    no1=nextvar[0]/2
    no2=nextvar[1]/2
    if(no1==current):
        current=no2     
        if(nextvar[1]%2==1):
            #cellwaypoint=singlearearoute(cells[current-1],width,points[(current-1)*2+1],points[(current-1)*2])
            cellwaypoint=singlearearoute(convexcells[current],width[current],waypoints[-1],points[(current)*2])
            #cellwaypoint.reverse()
        else:
            #cellwaypoint=singlearearoute(cells[current-1],width,points[(current-1)*2],points[(current-1)*2+1])
            cellwaypoint=singlearearoute(convexcells[current],width[current],waypoints[-1],points[(current)*2+1])
    else:
        current=no1
        if(nextvar[0]%2==1):
            #cellwaypoint=singlearearoute(cells[current-1],width,points[(current-1)*2+1],points[(current-1)*2])
            cellwaypoint=singlearearoute(convexcells[current],width[current],waypoints[-1],points[(current)*2])
            #cellwaypoint.reverse()
        else:
            #cellwaypoint=singlearearoute(cells[current-1],width,points[(current-1)*2],points[(current-1)*2+1])
            cellwaypoint=singlearearoute(convexcells[current],width[current],waypoints[-1],points[(current)*2+1])
    waypoints.extend(cellwaypoint)
    visited.append(current)
plotline(ax,waypoints,width[0])
waypoints.append(waypoints[0])
totalcost=timeconsume(waypoints,waypoints[0],waypoints[-1])
print ("total cost: %f" %totalcost)