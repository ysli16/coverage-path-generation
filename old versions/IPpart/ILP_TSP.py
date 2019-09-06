#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 16:57:50 2019
ref:http://examples.gurobi.com/traveling-salesman-problem/
@author: ysli
"""

import gurobipy as grb
import math
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
            if (pathvar[f].X!=0):
                if(f[0]<f[1]):
                    no1=f[0]/2
                    no2=f[1]/2
                    if(no1==no2):
                        print('explore cell %d' %(no1+1))
                    else:
                        if(f[0]%2==0 and f[1]%2==0):
                            print('from cell %d a to cell %d a' %((no1+1),(no2+1)))
                        elif(f[0]%2==0 and f[1]%2==1):
                            print('from cell %d a to cell %d b' %((no1+1),(no2+1)))
                        elif(f[0]%2==1 and f[1]%2==0):
                            print('from cell %d b to cell %d a' %((no1+1),(no2+1)))
                        else:
                            print('from cell %d b to cell %d b' %((no1+1),(no2+1)))
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
    
points=[(0,2),(1,2),(1,1),(2,0.5),(1,3),(2,3.5),(2,0.5),(3,1),(2,3.5),(3,3),(3,2),(4,2)]
num=len(points)
m=solveLP(points)
printsolution(m,m._vars)