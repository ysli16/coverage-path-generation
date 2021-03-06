#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 11:50:19 2019

@author: ysli

reference:
    codes for modelling TSP problem are from Gurobi Official Document
features:
    implement a variant of TSP to optimize visiting sequence of cells(use lazy constrains for subloop to simplify solving procedure)
    consider all vertices of each cell that are adjacent to other cells as nodes
    assume all nodes are connected
    consider Euclidean distance as cost of edges connected to different cells
    consider total time cost of the coverage path entering from vertice i, leaving from vertice j as the cost of the edge connected to vertice i and vertice j of this cell
    
    for path inside each cell, consider optimal angle as the one that minimize width of the cell,
    compare 4 possibilities(start from leftup,leftdown,rightup,rightdown)
    
    the optimal path has minimal traveling time, taking length of the path itself, 
    distance from end point of last cell to first point of this cell, 
    distance from this cell to its centroid and turning cost into consideration 
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import Point, Line, Polygon,Segment
import gurobipy as grb
import time


#find all Points that equals to the target
#return -1 for none, single point or a list of points
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
    
#split input polygon into two parts by input line
#return leftpolygon,rightpolygon or uppolygon,downpolygon or original polygon,none(if no intersection)
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
    
#remove the area of obstacle from the area of interest
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
#merge two adjacent cells
def mergecell(cell0,cell1,index0,index1):
    vertice0=cell0.vertices
    vertice1=cell1.vertices
    if(index0==0):
        vertice=vertice0[0:len(vertice0)-1]
    else:
        vertice=vertice0[index0:len(cell0.vertices)]+vertice0[0:index0-1]
    if(index1==0):
        vertice+=vertice1[0:len(vertice1)-1]
    else:
        vertice+=vertice1[index1:len(cell1.vertices)]+vertice1[0:index1-1]
    return Polygon(*vertice)

def addtoTSPmodel(newcell,nodes):
    node=newcell.centroid
    nodes.append([node.x,node.y])
    
def solveTSPLP(points):
    # Create variables
    m = grb.Model()
    vars = {}
    for i in range(num):
       for j in range(i+1):
           vars[i,j] = m.addVar(obj=TSPdistance(points[i],points[j]), vtype=grb.GRB.BINARY,
                              name='e'+str(i)+'_'+str(j))
           vars[j,i] = vars[i,j]
       m.update()
    # Add degree-2 constraint, and forbid loops   
    for i in range(num):
      m.addConstr(grb.quicksum(vars[i,j] for j in range(num)) == 2)
      vars[i,i].ub = 0
    m.update()
    m._vars = vars
    m.params.LazyConstraints = 1
    m.optimize(TSPsubtourelim)
    return m

def TSPdistance(point1,point2):
    distance=math.sqrt(pow(point1[0]-point2[0],2)+pow(point1[1]-point2[1],2))
    return distance

def TSPsubtourelim(model, where):
  if where == grb.GRB.callback.MIPSOL:
    selected = []
    # make a list of edges selected in the solution
    for i in range(num):
      sol = model.cbGetSolution([model._vars[i,j] for j in range(num)])
      selected += [(i,j) for j in range(num) if sol[j] > 0.5]
    # find the shortest cycle in the selected edge list
    tour = TSPsubtour(selected)
    if len(tour) < num:
      # add a subtour elimination constraint
      expr = 0
      for i in range(len(tour)):
        for j in range(i+1, len(tour)):
          expr += model._vars[tour[i], tour[j]]
      model.cbLazy(expr <= len(tour)-1)

def TSPsubtour(edges):
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
'''
LP model input:
vpoints=[
        [(1,0),(1,4)],
        [(1,0),(1,2),(2,0),(2,1)],
        [(1,2),(1,4),(2,3),(2,4)],
        [(2,0),(2,1),(3,0),(3,2)],
        [(2,3),(2,4),(3,2),(3,4)],
        [(3,0),(3,4)]
        ]
innercost=[
        [[0,20],[20,0]],
        [[0,10,10,10],[10,0,10,10],[10,10,0,10],[10,10,10,0]],
        [[0,10,10,10],[10,0,10,10],[10,10,0,10],[10,10,10,0]],
        [[0,10,10,10],[10,0,10,10],[10,10,0,10],[10,10,10,0]],
        [[0,10,10,10],[10,0,10,10],[10,10,0,10],[10,10,10,0]],
        [[0,10],[10,0]]
        ]
totalnum=0
num=[]
for i in vpoints:
    totalnum=totalnum+len(i)
    num.append(len(i))

'''

#compute data required by the LP model
#vpoints:each row records adjacent vertices of one cell
#index:the rank of current cell in cells
#pathindex:record the rank of the path in list possiblepaths responding to the minimal cost of current combination of vertices
def addtomodel(newcell,possiblepaths,cells,vpoints,innercost,pathindex):
    vpoints.append([])
    innercost.append([])
    pathindex.append([])
    index=len(vpoints)-1
    '''
    for newvertice in newcell.vertices:
        i=0
        for cell in cells:   
            if not(i==index):
                cross=cell.intersection(newvertice)
                if(cross):
                    vpoints[-1].append(newvertice)
                    break
            i=i+1
    '''
    for path in possiblepaths:
        vpoints[-1].append(Point(path[0]))
    cost=np.zeros((len(vpoints[-1]),len(vpoints[-1])))
    pindex=np.zeros((len(vpoints[-1]),len(vpoints[-1])))
    for i in range(len(vpoints[-1])):
        for j in range(i):
            mincost=float('inf')
            k=0
            for path in possiblepaths:
                pathcost=timeconsume(path,vpoints[-1][i],vpoints[-1][j])
                if(pathcost<mincost):
                    mincost=pathcost
                    index=k
                k=k+1           
            cost[i][j]=mincost
            cost[j][i]=cost[i][j]
            pindex[i][j]=int(index)
            if(index==0 or index==1):
                pindex[j][i]=int(1-index)
            else:
                pindex[j][i]=int(5-index)
    innercost[-1].extend(cost.tolist())
    pathindex[-1].extend(pindex.tolist())


#compute distance between two nodes. 
#consider Euclidean distance of two nodes belong to different cells
#consider path time cost if two nodes belong to the same cell
def distance(indexi,indexj):
    cellno1=indexi/4
    cellno2=indexj/4
    verticeno1=indexi%4
    verticeno2=indexj%4
    if not(cellno1==cellno2):
        distance=math.sqrt(pow(vpoints[cellno1][verticeno1][0]-vpoints[cellno2][verticeno2][0],2)+pow(vpoints[cellno1][verticeno1][1]-vpoints[cellno2][verticeno2][1],2))
    else:
        distance=innercost[cellno1][verticeno1][verticeno2]
    return distance

# Callback - use lazy constraints to eliminate sub-tours
def subtourelim(model, where):
  if where == grb.GRB.callback.MIPSOL:
    selected = []
    # make a list of edges selected in the solution
    for i in range(totalnum):
      sol = model.cbGetSolution([model._vars[i,j] for j in range(totalnum)])
      selected += [(i,j) for j in range(totalnum) if sol[j] > 0.5]
    # find the shortest cycle in the selected edge list
    tour = subtour(selected)
    if len(tour) < 2*len(vnum):
      # add a subtour elimination constraint
      expr = 0
      for i in range(len(tour)):
        for j in range(i+1, len(tour)):
          expr += model._vars[tour[i], tour[j]]
      model.cbLazy(expr <= len(tour)-1)
      
# Given a list of edges, finds the shortest subtour     
def subtour(edges):
  visited = [False]*totalnum
  cycles = []
  lengths = []
  selected = [[] for i in range(totalnum)]
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
    if not (len(thiscycle)==1):
        cycles.append(thiscycle)
        lengths.append(len(thiscycle))
    if sum(lengths) == len(edges)/2:
      break
  return cycles[lengths.index(min(lengths))]

def printsolution(m,pathvar):
    if m.status == grb.GRB.Status.OPTIMAL:
        print('totalCost: %g' % m.objVal)
        print('Path:')
        for j in range(totalnum):
            for i in range(j):
                if (round(pathvar[(i,j)].X,5)==1):
                    cellno1=i/4
                    cellno2=j/4
                    verticeno1=i%4
                    verticeno2=j%4
                    if(cellno1==cellno2):
                        print('explore cell %d, from vertice %d to %d' %(cellno1,verticeno1,verticeno2))
                    else:
                        print('from cell %d vertice %d to cell %d vertice %d' %((cellno1),(verticeno1),(cellno2),(verticeno2)))
    else:
        print('No feasible solution')
        
#formulate and solve LP model
def solveLP(points,startvar):
    # Create variables
    m = grb.Model()
    vars = {}
    for i in range(totalnum):
       for j in range(i+1):           
           vars[i,j] = m.addVar(obj=distance(i,j), vtype=grb.GRB.BINARY,
                              name='e'+str(i)+'_'+str(j))
           if(startvar):
               vars[i,j].start=0
           vars[j,i] = vars[i,j]
       m.update()
    for i in range(totalnum):     
      vars[i,i].ub = 0
    
    # Add degree-2 constraint 
    for y in range(totalnum):
        vars[y]=m.addVar(vtype=grb.GRB.BINARY,name='y')
        if(startvar):
            vars[y].start=0
    m.update()
    for i in range(totalnum):
        m.addConstr((vars[i]==0)>>(grb.quicksum(vars[i,j] 
                    for j in range(totalnum)
                    )==0),name='unselected '+str(i))
        
        m.addConstr((vars[i]==1)>>(grb.quicksum(vars[i,j] 
                    for j in range(totalnum)
                    )==2),name='selected '+str(i))
    for k in range(len(vnum)): 
        indexrange=[4*k,4*k+3]
        m.addConstr(grb.quicksum(vars[i,j] 
                    for i in range(indexrange[0],indexrange[1]+1) 
                    for j in range(indexrange[0],indexrange[1]+1))==2,name='inside '+str(k))
    m.update()
    
    if(startvar):
        visited=[]
        currentc=0
        for j in range(num):
            vars[4*j,4*j+1].start=1
            vars[4*j+1,4*j].start=1
            vars[4*j].start=1
            vars[4*j+1].start=1
            nextc=findnextTSP(startvar,visited,currentc)
            vars[4*currentc+1,4*nextc].start=1
            vars[4*nextc,4*currentc+1].start=1
            currentc=nextc
            visited.append(currentc)

        m.update() 
    m._vars = vars           
    m.params.LazyConstraints = 1
    m.params.MIPFocus=3
    m.params.ImproveStartTime=300
    m.params.MIPgap=0.005
    m.optimize(subtourelim)
    return m
def findnextTSP(pathvar,visited,current):
    for p in pathvar:
        if (round(pathvar[p].X,5)==1):
            if(p[0]==current):
                if not(p[1] in visited):
                    return p[1]
            elif(p[1]==current):
                if not(p[0] in visited):
                    return p[0]
#given a sequence of points, return the indexs of left,bottom,right,top extremities in a list
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

#form a list of vertical lines given a list of points
def formperpenicularlines(xpos):
    xpos=xpos.reshape(len(xpos),1)
    xpts=np.insert(xpos,1,[0],axis=1)
    xpts=map(Point,xpts)
    cutlines=[]
    for pt in xpts:
        cutline=xaxis.perpendicular_line(pt)
        cutlines.append(cutline)
    return cutlines

#find coverage path inside a convex polygon, return two paths, start from leftup and leftdown
def createallwaypoint(rboundary,inputwidth,height):
    bounds=rboundary.bounds
    cutnum=float(math.ceil(round((round(bounds[2],6)-round(bounds[0],6))/inputwidth,6)))
    width=(round(bounds[2],6)-round(bounds[0],6))/cutnum 
    cutpos=np.arange(bounds[0]+width,bounds[2]+width,width).astype(np.double)
    cutpos=np.around(cutpos,6)
    pointpos=np.arange(bounds[0]+width/2,bounds[2]+width/2,width).astype(np.double)
    pointpos=np.around(pointpos,6)
    cutlines=formperpenicularlines(cutpos)
    #form sub-polygons by cutlines
    subpolygons=[]
    remain=rboundary
    for line in cutlines:
        subpolygon,remain=splitpolygon(remain,line)
        subpolygons.append(subpolygon)
        if not (remain):
            break
    #generate waypoints in subpolygons
    i=0
    waypoints1=[]
    waypoints2=[]
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
            if(upmax-downmin<=height):
                waypoint.append(Point(pointpos[i],(upmax+downmin)/2))
            else:
                if(upmax-upmin<=height/2):
                    waypoint.append(Point(pointpos[i],upmax-height/2))
                else:
                    cutpt=Point(pointpos[i],upmax-height/2)
                    cut=yaxis.perpendicular_line(cutpt)
                    crosspoint=polygon.intersection(cut)
                    if(crosspoint[0].x<crosspoint[1].x):
                        leftcross=crosspoint[0]
                        rightcross=crosspoint[1]
                    else:
                        leftcross=crosspoint[1]
                        rightcross=crosspoint[0]
                    if(leftcross.x<=pointpos[i] and rightcross.x>=pointpos[i]):
                        waypoint.append(Point(pointpos[i],upmax-height/2))
                    elif(leftcross.x>pointpos[i]):
                        waypoint.append(leftcross)
                        if(upmin-downmax<=height):
                            waypoint.append(Point(pointpos[i],(upmin+downmax)/2))
                        else:
                            waypoint.append(Point(pointpos[i],upmin-height/2))                                
                    else:
                        waypoint.append(rightcross)
                        if(upmin-downmax<=height):
                            waypoint.append(Point(pointpos[i],(upmin+downmax)/2))
                        else:
                            waypoint.append(Point(pointpos[i],upmin-height/2)) 
                if(downmax-downmin<=height/2):
                    waypoint.append(Point(pointpos[i],downmin+height/2))
                else:
                    cutpt=Point(pointpos[i],downmin+height/2)
                    cut=yaxis.perpendicular_line(cutpt)
                    crosspoint=polygon.intersection(cut)
                    if(crosspoint[0].x<crosspoint[1].x):
                        leftcross=crosspoint[0]
                        rightcross=crosspoint[1]
                    else:
                        leftcross=crosspoint[1]
                        rightcross=crosspoint[0]
                    if(leftcross.x<=pointpos[i] and rightcross.x>=pointpos[i]):
                        waypoint.append(Point(pointpos[i],downmin+height/2))
                    elif(leftcross.x>pointpos[i]):
                        if(upmin-downmax<=height):
                            waypoint.append(Point(pointpos[i],(upmin+downmax)/2))
                        else:
                            waypoint.append(Point(pointpos[i],downmax+height/2))
                        waypoint.append(leftcross)
                    else:
                        if(upmin-downmax<=height):
                            waypoint.append(Point(pointpos[i],(upmin+downmax)/2))
                        else:
                            waypoint.append(Point(pointpos[i],downmax+height/2)) 
                        waypoint.append(rightcross)
        elif(len(lines)==1):#formed by one cut
            index=findbound(vertices)
            leftindex=index[0]
            line=lines[0]
            ymax=vertices[index[3]].y
            ymin=vertices[index[1]].y
            xmin=vertices[index[0]].x
            xmax=vertices[index[2]].x
            if(line.contains(vertices[leftindex])):#cutline at left
                if(ymax-ymin<=height):
                    waypoint.append(Point(xmax-width/2,(ymin+ymax)/2))
                else:
                    upcutpt=Point(pointpos[i],ymax-height/2)
                    upcut=yaxis.perpendicular_line(upcutpt)
                    downcutpt=Point(pointpos[i],ymin+height/2)
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
                        if((uppt.y-downpt.y)<=height):
                            waypoint.append(Segment(uppt,downpt).midpoint)
                        else:
                            waypoint.append(Point(uppt.x,uppt.y-height/2))
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
                        if((uppt.y-downpt.y)<=height):
                            waypoint.append(Segment(uppt,downpt).midpoint)
                        else:
                            waypoint.append(Point(downpt.x,downpt.y+height/2))
                        crosspoint=polygon.intersection(downcut)
                        if(crosspoint[0].x<crosspoint[1].x):
                            waypoint.append(crosspoint[0])
                        else:
                            waypoint.append(crosspoint[1])                                 
            else:#cutline at right
                if(ymax-ymin<=height):
                    waypoint=[Point(xmin+width/2,(ymin+ymax)/2)]
                else:
                    upcutpt=Point(pointpos[i],ymax-height/2)
                    upcut=yaxis.perpendicular_line(upcutpt)
                    downcutpt=Point(pointpos[i],ymin+height/2)
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
                        if((uppt.y-downpt.y)<=height):
                            waypoint.append(Segment(uppt,downpt).midpoint)
                        else:
                            waypoint.append(Point(uppt.x,uppt.y-height/2))
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
                        if((uppt.y-downpt.y)<=height):
                            waypoint.append(Segment(uppt,downpt).midpoint)
                        else:
                            waypoint.append(Point(downpt.x,downpt.y+height/2))
                        crosspoint=polygon.intersection(downcut)
                        if(crosspoint[0].x>crosspoint[1].x):
                            waypoint.append(crosspoint[0])
                        else:
                            waypoint.append(crosspoint[1])             
        else:#width of the area less than sensor coverage range
            bound=polygon.bounds
            if(bound[3]-bound[1]<height):
                waypoint=[Point(pointpos[i],(bound[1]+bound[3])/2)]
            else:
                waypoint=[Point(pointpos[i],bound[3]-height/2),Point(pointpos[i],bound[1]+height/2)]
        i=i+1
        if(i%2==1):
            waypoint.reverse()
        waypoints1.extend(waypoint)
        waypoint.reverse()
        waypoints2.extend(waypoint)
    return waypoints1,waypoints2

#eliminate the error caused by rotation
def roundpolygon(polygon):
    newvertices=[]
    for vertice in polygon.vertices:
        newvertice=Point(round(vertice.x,6),round(vertice.y,6))
        newvertices.append(newvertice)
    newpolygon=Polygon(*newvertices)
    return newpolygon

#find the rotation angle that minimize the width of the cell  
def findoptimalangle(cell):
    minwidth=float('inf')
    if(cell.is_convex()):
        for side in cell.sides:
            angle=math.atan(side.slope)
            rcell=cell.rotate(math.pi/2-angle)
            bound=rcell.bounds
            width=(bound[2]-bound[0])
            if(width<minwidth):
                minwidth=width
                optimalangle=math.pi/2-angle
        return optimalangle
    else:
        return 0

#store all possible paths
def possiblepath(cell,inputwidth,height):
    angle=float(findoptimalangle(cell))
    possiblepaths=[]
    rcell=cell.rotate(angle)
    rcell=roundpolygon(rcell)
    waypoint1,waypoint2=createallwaypoint(rcell,inputwidth,height)
    for i in range(len(waypoint1)):
        waypoint1[i]=waypoint1[i].rotate(-angle)
    possiblepaths.append(waypoint1)   
    
    rvwaypoint=waypoint1[:]
    rvwaypoint.reverse()
    possiblepaths.append(rvwaypoint)
        
    for i in range(len(waypoint2)):
        waypoint2[i]=waypoint2[i].rotate(-angle)
    possiblepaths.append(waypoint2)
    
    rvwaypoint=waypoint2[:]
    rvwaypoint.reverse()
    possiblepaths.append(rvwaypoint)
    return possiblepaths

#estimate time cost of the path (take number of turns into consideration)
def timeconsume(waypoint,rentrance,rexit):
    length=0
    for i in range(len(waypoint)-1):
        length=length+math.sqrt(pow(waypoint[i].x-waypoint[i+1].x,2)+pow(waypoint[i].y-waypoint[i+1].y,2))
    length=length+waypoint[0].distance(rentrance)
    length=length+waypoint[-1].distance(rexit)
    turnnum=len(waypoint)-1
    return length/velocity+turnnum*turnpanalty

#get sequence from the model solution
def findnext(pathvar,visited,current):
    for i in range(totalnum):
        if (round(pathvar[(i,current)].X,5)==1):
            if not(i in visited):
                return i
            
#convert Point to array
def pointtoarray(vertices):
    point=[]
    for vertice in vertices:
        point.append([vertice.x,vertice.y])
    return np.array(point)    

#plot cells                     
def plotpolygon(ax,boundary):    
    vertices=list(boundary.vertices)
    vertices.append(vertices[0])
    vertice=pointtoarray(vertices)
    ax.plot(vertice[:,0],vertice[:,1],linewidth=2,color='black',linestyle=':')
    
#plot boundarys
def plotboundary(ax,boundary):    
    vertices=list(boundary.vertices)
    vertices.append(vertices[0])
    vertice=pointtoarray(vertices)
    ax.plot(vertice[:,0],vertice[:,1],linewidth=2,color='black',linestyle='-')    

#plot waypoints
def plotline(ax,rwaypoints):
    point=pointtoarray(rwaypoints) 
    ax.plot(point[:,0],point[:,1],linewidth=1,color='g')

#the following are some test samples, consisting of boundary part and obstacle part. 
#Input coordinates are in counterclockwise
#boundary part
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

ratio=1000/2.4
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

'''
#simplified Lake Mascoma
boundary=Polygon((0,1),(4,0),(6.5,0.5),(8,1.8),(4.6,3.5),(3,3),(2.8,1),(1,1.5))
'''
#obstacle part
obstacles=[]#store all obstacles(deeper areas) as polygons
allinnerpt=[]#store all coordinates([x,y]) of obstacles and deeper areas
obsleft=[]#store all leftmost position of obstacles and deeper areas
obsright=[]#store all rightmost position of obstacles and deeper areas

'''
#sample1,2
obsnum=1
temppolygon=Polygon((1,2),(2,1),(3,2),(2,3))
obstacles.append(temppolygon)
allinnerpt.extend([[1,2],[2,1],[3,2],[2,3]])
obsleft=[1]
obsright=[3]
iscell=[False]
'''
'''
#sample3
obsnum=1
temppolygon=Polygon((0.8,2),(0.8,1),(2,1),(2,2))
obstacles.append(temppolygon)
allinnerpt.extend([[0.8,2],[0.8,1],[2,1],[2,2]])
obsleft=[0.8]
obsright=[2]
iscell=[False]
'''
'''
#sample4
obsnum=1
temppolygon=Polygon((0.8,2),(0.8,1),(2.3,0.5),(3,2),(2,3))
obstacles.append(temppolygon)
allinnerpt.extend([[0.8,2],[0.8,1],[2.3,0.5],[3,2],[2,3]])
obsleft=[0.8]
obsright=[3]
iscell=[False]
'''
'''
#sample5
obsnum=1
temppolygon=Polygon((1,3),(2,1),(4,2),(2,3))
obstacles.append(temppolygon)
allinnerpt.extend([[1,3],[2,1],[4,2],[2,3]])
obsleft=[1]
obsright=[4]
iscell=[True]
'''
'''
#sample6
obsnum=2
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
obsnum=4
temppt=[(8.9,2.5),(9.1,2.4),(9.25,2.5),(9,2.7)]
temppt=np.array(temppt)
temppt=temppt*ratio
temppolygon=map(Point,temppt)
temppolygon=Polygon(*temppolygon)
obstacles.append(temppolygon)
allinnerpt.extend(temppt)

temppt=[(9.3,2.55),(9.4,2.4),(9.6,2.4),(9.6,2.5),(9.5,2.6)]
temppt=np.array(temppt)
temppt=temppt*ratio
temppolygon=map(Point,temppt)
temppolygon=Polygon(*temppolygon)
obstacles.append(temppolygon)
allinnerpt.extend(temppt)

temppt=[(7,2.4),(7.15,0.9),(8.45,0.5),(9.2,0.5),(10.15,0.9),(10.15,1.7),(8.7,1.9),(7.25,2.6)]
temppt=np.array(temppt)
temppt=temppt*ratio
temppolygon=map(Point,temppt)
temppolygon=Polygon(*temppolygon)
obstacles.append(temppolygon)
allinnerpt.extend(temppt)

temppt=[(11.75,1.5),(12.4,1),(12.75,1.25),(12.45,2),(11.9,1.85)]
temppt=np.array(temppt)
temppt=temppt*ratio
temppolygon=map(Point,temppt)
temppolygon=Polygon(*temppolygon)
obstacles.append(temppolygon)
allinnerpt.extend(temppt)
obsleft=np.array([8.9,9.3,7,11.75])*ratio
obsright=np.array([9.25,9.6,10.15,12.75])*ratio
iscell=[False,False,True,True]

'''
obsleft=np.array([8.9,9.3])*ratio
obsright=np.array([9.25,9.6,])*ratio
iscell=[False,False]
'''
#simplified Lake Mascoma
'''
obsnum=2
temppolygon=Polygon((3.5,1.5),(4,0.8),(5,1.5),(5,2.2))
obstacles.append(temppolygon)
allinnerpt.extend([[3.5,1.5],[4,0.8],[5,1.5],[5,2.2]])
temppolygon=Polygon((3.8,3.1),(3.8,2.7),(4.6,2.7),(4.6,3))
obstacles.append(temppolygon)
allinnerpt.extend([[3.8,3.1],[3.8,2.7],[4.6,2.7],[4.6,3]])
obsleft=[3.5,3.8]
obsright=[5,4.6]
iscell=[True,False]
'''

#user settings
postmerge=True
smartstart=True
width1=16.8
height1=12.5
width2=84
height2=62.5
velocity=1
turnpanalty=10

#algorithm begins
starttime=time.time()
allinnerpt=np.array(allinnerpt)
allinnerpt=allinnerpt[np.lexsort((allinnerpt[:,1],allinnerpt[:,0]))]
xaxis=Line((0,0),(1,0))
yaxis=Line((0,0),(0,1))
remain=boundary
cells=[]
vpoints=[]
innercost=[]
pathindex=[]
#area decomposition
for pt in allinnerpt:
    isleft=False
    isright=False
    if(pt[0] in obsleft):
        isleft=True
    if(pt[0] in obsright):
        isright=True
    cut=xaxis.perpendicular_line(Point(pt))
    cross=remain.intersection(cut)
    crosspt=[]
    for item in cross:
        if(isinstance(item,Point)):
            crosspt.append(item)
        else:          
            if(len(crosspt)==0):
                crosspt.append(item.p2)
                crosspt.append(item.p1)
            else:
                i=0
                for cpt in crosspt:
                    if cpt.y<item.p1.y:
                        i=i+1
                    else:
                        break
                crosspt.insert(i,item.p1)
                crosspt.insert(i,item.p2)
    if(isright):
        if not (Point(pt) in crosspt):
            continue
        elif not(Point(pt) in cross):
            index=crosspt.index(Point(pt))
            seg1=Segment(crosspt[index-1],crosspt[index])
            seg2=Segment(crosspt[index+1],crosspt[index+2])
            cell,remain=splitpolygon(remain,seg1)
            cells.append(cell)
            cell,remain=splitpolygon(remain,seg2)
            cells.append(cell)
        else:
            index=crosspt.index(Point(pt))
            seg1=Segment(crosspt[index-1],crosspt[index])
            seg2=Segment(crosspt[index],crosspt[index+1])
            cell,remain=splitpolygon(remain,seg1)
            cells.append(cell)
            cell,remain=splitpolygon(remain,seg2)
            cells.append(cell)  
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
            else:
                remain=cell
            for obs in obstacles:
                if(obs.contains(Point(pt))):
                    break
            remain=minuspolygon(remain,obs,pt)
    if not (postmerge):        
        if not (isleft or isright):
            for obs in obstacles:
                if(obs.contains(Point(pt))):
                    break
            if (obs.angles[Point(pt)]<math.pi):
                index=findpoint(crosspt,Point(pt))
                if(index%2==0):
                    seg=Segment(crosspt[index],crosspt[index+1])
                else:
                    seg=Segment(crosspt[index-1],crosspt[index])
                cell,remain=splitpolygon(remain,seg)
                cells.append(cell)  
            

cells.append(remain)
width=np.ones(len(cells))*width1
height=np.ones(len(cells))*height1
#add deeper areas
for i in range(obsnum):
    if(iscell[i]):
        cells.append(obstacles[i])
        width=np.append(width,width2)
        height=np.append(height,height2)

#split concave cell into convex cells 
cindex=0
count=0 
convexcells=list(cells) 
for cell in cells:
    cutpt1=[]
    cutpt2=[]
    if not(cell.is_convex()):
        for vertice in cell.vertices:
            if(cell.angles[vertice]>math.pi):
                cutpt1.append(vertice)
            else:
                cutpt2.append(vertice)
        if(len(cutpt1)>len(cutpt2)):
            pts=pointtoarray(cutpt2)
        else:
            pts=pointtoarray(cutpt1)
        pts = pts[pts[:,0].argsort()]
        if(postmerge):
            i=0
            for i in range(len(pts)+1):       
                if(i==len(pts)):
                    if(i>1):
                        cell0=convexcells[cindex+count-1]
                        cell1=convexcells[cindex+count]
                        if((findoptimalangle(cell0)==0) and (findoptimalangle(cell1)==0)):
                            seg=cell1.intersection(cell0)
                            index1=findpoint(cell1.vertices,seg[0].p1)
                            index0=findpoint(cell0.vertices,seg[0].p2)
                            if not (index1==-1 or index0==-1):
                                mergedcell=mergecell(cell0,cell1,index0,index1)
                                convexcells.remove(convexcells[cindex+count-1])
                                convexcells.insert(cindex+count-1,mergedcell)
                                convexcells.remove(convexcells[cindex+count])
                                width=np.delete(width,cindex+count)
                                height=np.delete(height,cindex+count)
                                count=count-1
                    else: 
                        continue
                else:
                    pt=pts[i]
                    cut=yaxis.parallel_line(Point(pt))
                    cross=convexcells[cindex+count].intersection(cut)
                    if(len(cross)==2):
                        cell1,cell2=splitpolygon(convexcells[cindex+count],cut)
                        if(len(pts)==1):
                            if not((findoptimalangle(cell1)==0) and (findoptimalangle(cell2)==0)):
                                convexcells.remove(convexcells[cindex+count])
                                convexcells.insert(cindex+count,cell2)
                                convexcells.insert(cindex+count,cell1) 
                                width=np.insert(width,cindex+count,width[cindex+count])
                                height=np.insert(height,cindex+count,height[cindex+count]) 
                                count=count+1
                        else:
                            if(i>0):
                                cell0=convexcells[cindex+count-1]
                                seg=cell1.intersection(cell0)
                                if(len(seg)):
                                    if(isinstance(seg[0],Segment)):
                                        if ((findoptimalangle(cell1)==0) and (findoptimalangle(cell0)==0)):                        
                                            index1=findpoint(cell1.vertices,seg[0].p1)
                                            index0=findpoint(cell0.vertices,seg[0].p2)
                                            if not (index1==-1 or index0==-1):
                                                mergedcell=mergecell(cell0,cell1,index0,index1)
                                                convexcells.remove(convexcells[cindex+count-1])
                                                convexcells.insert(cindex+count-1,mergedcell)
                                                convexcells.remove(convexcells[cindex+count])
                                                convexcells.insert(cindex+count,cell2)                                           
                                                continue
                            convexcells.remove(convexcells[cindex+count])
                            convexcells.insert(cindex+count,cell2)
                            convexcells.insert(cindex+count,cell1) 
                            width=np.insert(width,cindex+count,width[cindex+count])
                            height=np.insert(height,cindex+count,height[cindex+count]) 
                            count=count+1
                    elif(len(cross)==3):
                        seg1=Segment(cross[0],cross[1])
                        cell1,cell2=splitpolygon(convexcells[cindex+count],seg1)    
                        convexcells.remove(convexcells[cindex+count])
                        convexcells.insert(cindex+count,cell2)
                        convexcells.insert(cindex+count,cell1) 
                        width=np.insert(width,cindex+count,width[cindex+count])
                        height=np.insert(height,cindex+count,height[cindex+count])
                        count=count+1
                        seg2=Segment(cross[1],cross[2])
                        cell1,cell2=splitpolygon(convexcells[cindex+count],seg2)    
                        convexcells.remove(convexcells[cindex+count])
                        convexcells.insert(cindex+count,cell2)
                        convexcells.insert(cindex+count,cell1)  
                        width=np.insert(width,cindex+count,width[cindex+count])
                        height=np.insert(height,cindex+count,height[cindex+count])
                        count=count+1
        else:
            for pt in pts:
                cut=yaxis.parallel_line(Point(pt))
                cross=convexcells[cindex+count].intersection(cut)
                if(len(cross)==2):
                    cell1,cell2=splitpolygon(convexcells[cindex+count],cut)    
                    convexcells.remove(convexcells[cindex+count])
                    convexcells.insert(cindex+count,cell2)
                    convexcells.insert(cindex+count,cell1) 
                    width=np.insert(width,cindex+count,width[cindex+count])
                    height=np.insert(height,cindex+count,height[cindex+count]) 
                    count=count+1                         
                elif(len(cross)==3):
                    seg1=Segment(cross[0],cross[1])
                    cell1,cell2=splitpolygon(convexcells[cindex+count],seg1)    
                    convexcells.remove(convexcells[cindex+count])
                    convexcells.insert(cindex+count,cell2)
                    convexcells.insert(cindex+count,cell1) 
                    width=np.insert(width,cindex+count,width[cindex+count])
                    height=np.insert(height,cindex+count,height[cindex+count]) 
                    count=count+1
                    seg2=Segment(cross[1],cross[2])
                    cell1,cell2=splitpolygon(convexcells[cindex+count],seg2)    
                    convexcells.remove(convexcells[cindex+count])
                    convexcells.insert(cindex+count,cell2)
                    convexcells.insert(cindex+count,cell1)  
                    width=np.insert(width,cindex+count,width[cindex+count])
                    height=np.insert(height,cindex+count,height[cindex+count]) 
                    count=count+1
    cindex=cindex+1
#get TSP solution as MIP start
if(smartstart): 
    nodes=[]
    for cell in convexcells:
        addtoTSPmodel(cell,nodes)   
    num=len(nodes)
    tspmodel=solveTSPLP(nodes)  
    
#get all possible coverage path inside each cell
#compute data required by the LP model
allpossiblepath=[]
i=0
for cell in convexcells:
    path=possiblepath(cell,width[i],height[i])
    allpossiblepath.append(path)
    addtomodel(cell,path,convexcells,vpoints,innercost,pathindex)
    i=i+1
totalnum=0
vnum=[]
for i in vpoints:
    totalnum=totalnum+len(i)
    vnum.append(len(i))

#solve model
if(smartstart):
    startvar=tspmodel._vars
    m=solveLP(vpoints,startvar)
else:
    m=solveLP(vpoints,None)

printsolution(m,m._vars)

#generate whole path
waypoints=[]
for i in range(4):
    for j in range(i):
        if(round(m._vars[(i,j)].X,5)==1):
            current=j
visited=[current]
next1=findnext(m._vars,visited,current)
if(current/4==next1/4):
    cellno=current/4
    verticeno1=current%4
    verticeno2=next1%4
    cellwaypoint=allpossiblepath[cellno][int(pathindex[cellno][verticeno1][verticeno2])]
    waypoints.extend(cellwaypoint)
current=next1
visited.append(current)
next1=findnext(m._vars,visited,current)
while(next1):    
    if(current/4==next1/4):
        cellno=current/4
        verticeno1=current%4
        verticeno2=next1%4
        cellwaypoint=allpossiblepath[cellno][int(pathindex[cellno][verticeno1][verticeno2])]
        waypoints.extend(cellwaypoint)
    current=next1
    visited.append(current)
    next1=findnext(m._vars,visited,current)
    
computetime=time.time()-starttime

#plot result
fig=plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_aspect('equal', 'box')
for cell in convexcells:
    plotpolygon(ax,cell)
plotboundary(ax,boundary)
for i in range(obsnum):
    if not(iscell[i]):
        plotboundary(ax,obstacles[i])
        
plotline(ax,waypoints)
waypoints.append(waypoints[0])
totaltime=timeconsume(waypoints,waypoints[0],waypoints[-1])
turnnum=len(waypoints)-1

print "distance: %f" %(totaltime-turnnum*turnpanalty)
print "number of turns:%f" %turnnum
print "algorithm runtime: %f s" %computetime