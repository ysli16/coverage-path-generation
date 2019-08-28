#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 16:57:50 2019
ref:http://examples.gurobi.com/traveling-salesman-problem/
@author: ysli
"""

import gurobipy as grb
import math
def findcno(index):
    count=0
    for i in range(len(num)):
        count=count+num[i]
        if(count>index):
            return i
def findvno(cellno,index):
    count=0
    for i in range(cellno):
        count=count+num[i]
    return index-count
def findcindex(cellno):
    count=0
    for i in range(cellno):
        count=count+num[i]
    return [count,count+num[cellno]-1]
def distance(indexi,indexj):
    cellno1=findcno(indexi)
    cellno2=findcno(indexj)
    verticeno1=findvno(cellno1,indexi)
    verticeno2=findvno(cellno2,indexj)
    if not(cellno1==cellno2):
        distance=math.sqrt(pow(points[cellno1][verticeno1][0]-points[cellno2][verticeno2][0],2)+pow(points[cellno1][verticeno1][1]-points[cellno2][verticeno2][1],2))
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
    if len(tour) < len(num):
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
    cycles.append(thiscycle)
    lengths.append(len(thiscycle))
    if sum(lengths) >= len(num):
      break
  return cycles[lengths.index(min(lengths))]

def printsolution(m,pathvar):
    if m.status == grb.GRB.Status.OPTIMAL:
        print('totalCost: %g' % m.objVal)
        print('Path:')
        for j in range(totalnum):
            for i in range(j):
                if (pathvar[(i,j)].X!=0):
                    cellno1=findcno(i)
                    cellno2=findcno(j)
                    verticeno1=findvno(cellno1,i)
                    verticeno2=findvno(cellno2,j)
                    if(cellno1==cellno2):
                        print('explore cell %d, from vertice %d to %d' %(cellno1,verticeno1,verticeno2))
                    else:
                        print('from cell %d vertice %d to cell %d vertice %d' %((cellno1),(verticeno1),(cellno2),(verticeno2)))
    else:
        print('No feasible solution')
def solveLP(points):
    # Create variables
    m = grb.Model()
    vars = {}
    for i in range(totalnum):
       for j in range(i+1):           
           vars[i,j] = m.addVar(obj=distance(i,j), vtype=grb.GRB.BINARY,
                              name='e'+str(i)+'_'+str(j))
           vars[j,i] = vars[i,j]
       m.update()
    for i in range(totalnum):     
      vars[i,i].ub = 0
    
    # Add degree-2 constraint 
    for y in range(totalnum):
        vars[y]=m.addVar(vtype=grb.GRB.BINARY,name='y')
    m.update()
    for i in range(totalnum):
        m.addConstr((vars[i]==1)>>(grb.quicksum(vars[i,j] 
                    for j in range(totalnum)
                    )==0))
        
        m.addConstr((vars[i]==0)>>(grb.quicksum(vars[i,j] 
                    for j in range(totalnum)
                    )==2))
    #all egdes inside each cell must be visited at least once
    for cellno in range(len(num)): 
        indexrange=findcindex(cellno)
        m.addConstr(grb.quicksum(vars[i,j] 
                    for i in range(indexrange[0],indexrange[1]+1) 
                    for j in range(indexrange[0],indexrange[1]+1))>=1)
    m.update()
    m._vars = vars
    m.params.LazyConstraints = 1
    m.optimize(subtourelim)
    return m
points=[
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
for i in points:
    totalnum=totalnum+len(i)
    num.append(len(i))
    
m=solveLP(points)
printsolution(m,m._vars)