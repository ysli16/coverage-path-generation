#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 06:43:04 2019

@author: ysli
"""

import gurobipy as grb
import math
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
    m=grb.Model('ILP_CPP')
    modes=[1,2,3,4]
    allpath=grb.tuplelist()
    relax=grb.tuplelist()
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
        relax.append((cell1,'a'))
        relax.append((cell1,'b'))
    pathvar=m.addVars(allpath,vtype=grb.GRB.BINARY, name="path") 
    relaxvar=m.addVars(relax,vtype=grb.GRB.INTEGER,name="r")
    m.setObjective(pathvar.prod(cost), grb.GRB.MINIMIZE)   
    for cell1 in area:
        m.addConstr(
            sum(pathvar.select(cell1,'*',1))+
            sum(pathvar.select(cell1,'*',2))+
            sum(pathvar.select('*',cell1,1))+
            sum(pathvar.select('*',cell1,3))==1+2*relaxvar[cell1,'a'],str(cell1)+'a')
        m.addConstr(
            sum(pathvar.select('*',cell1,2))+
            sum(pathvar.select(cell1,'*',3))+
            sum(pathvar.select('*',cell1,4))+
            sum(pathvar.select(cell1,'*',4))==1+2*relaxvar[cell1,'b'],str(cell1)+'b')       
    m.optimize()
    return (m,pathvar)

m,pathvar=solveLP(area,coord)
printsolution(m,pathvar)