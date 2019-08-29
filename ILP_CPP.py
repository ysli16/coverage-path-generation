#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 06:43:04 2019

@author: ysli
"""

import gurobipy as grb
import math

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
    return (m,m._vars)
cellcost=[0.25,4/3,4/3,4]
ane=[
     [1,0,0,0],
     [1,1,1,0],
     [0,1,1,1],
     [0,0,0,1]]
kn=[1,1,1,1]
m,var=solveLP(cellcost,ane,kn)
printsolution(m,var)