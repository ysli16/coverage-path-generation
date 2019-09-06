# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import Point, Line, Polygon,Segment
global velocity
global turnpanalty
global num

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
#estimate time comsumption   
def timeconsume(waypoint,rentrance,rexit):
    length=0
    for i in range(len(waypoint)-1):
        length=length+math.sqrt(pow(waypoint[i].x-waypoint[i+1].x,2)+pow(waypoint[i].y-waypoint[i+1].y,2))
    length=length+waypoint[0].distance(rentrance)
    length=length+waypoint[len(waypoint)-1].distance(rexit)
    turnnum=len(waypoint)-1
    return length/velocity+turnnum*turnpanalty
#rearrange the sequence of boundry points so that they form a polygon
def formpolygon(coord):
    #rearrange coordinate(according to x)
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
def formperpenicularlines(xpos):
    xpos=xpos.reshape(len(xpos),1)
    xpts=np.insert(xpos,1,[0],axis=1)
    xpts=map(Point,xpts)
    cutlines=[]
    for pt in xpts:
        cutline=xaxis.perpendicular_line(pt)
        cutlines.append(cutline)
    return cutlines
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
def findpoint(points,goal):
    i=0
    for point in points:
        if(point.equals(goal)):
            return i
        else:
            i=i+1
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
            for side in sides:
                if(side.contains(upcross)):
                    if(upcross.x<side.p1.x or(upcross.x==side.p1.x and upcross.x>side.p2.x)):
                        rightupindex=findpoint(vertices,side.p1)
                        leftupindex=findpoint(vertices,side.p2)
                    else:
                        leftupindex=findpoint(vertices,side.p1)
                        rightupindex=findpoint(vertices,side.p2)
                if(side.contains(downcross)):
                    if(downcross.x<side.p1.x or(downcross.x==side.p1.x and downcross.x>side.p2.x)):
                        rightdownindex=findpoint(vertices,side.p1)
                        leftdownindex=findpoint(vertices,side.p2)
                    else:
                        leftdownindex=findpoint(vertices,side.p1)
                        rightdownindex=findpoint(vertices,side.p2)
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
    
def pointtolist(vertices):
    pointx=[]
    pointy=[]
    for vertice in vertices:
        pointx.append(vertice.x)
        pointy.append(vertice.y)
    return pointx,pointy
def createallwaypoint(rboundary):
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
            leftupindex=lineindex[0]
            leftdownindex=(lineindex[0]+1)%len(vertices)
            rightupindex=(lineindex[1]+1)%len(vertices)
            rightdownindex=lineindex[1]  
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

velocity=1
turnpanalty=10
'''
#input coordinates of boundary points
num=int(input("Enter the number of boundary points:"))
bptlist=[]
for i in range(num):
    longitude_in=float(input("Enter boundry longitude:"))
    latitude_in=float(input("Enter boundry latitude:"))
    bpt=(longitude_in, latitude_in)
    bptlist.append(bpt)
'''
num=4
bptlist=[(0,4),(3,0),(5,7),(8,4)]

p1x=float(input("Enter entrance point1x"))
p1y=float(input("Enter entrance point1y"))
p2x=float(input("Enter entrance point2x"))
p2y=float(input("Enter entrance point2y"))
Entrance=Line((p1x,p1y),(p2x,p2y))
p1x=float(input("Enter exit point1x"))
p1y=float(input("Enter exit point1y"))
p2x=float(input("Enter exit point2x"))
p2y=float(input("Enter exit point2y"))
Exit=Line((p1x,p1y),(p2x,p2y))
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
#make sure vertices are in ccw order
if(boundary.area<0):
    vertices=boundary.vertices
    vertices.reverse()
    boundary=Polygon(*vertices)

inputwidth=float(input("Enter width of the route:"))  
xaxis=Line((0,0),(1,0))
yaxis=Line((0,0),(0,1))
mintime=float('inf')
optimalangle=0
for angle in range(0,180,10):
    angle=float(angle)
    rboundary=boundary.rotate(angle/180*math.pi)
    bounds=rboundary.bounds
    cutnum=float(math.ceil(round((round(bounds[2],6)-round(bounds[0],6))/inputwidth,6)))
    width=(round(bounds[2],6)-round(bounds[0],6))/cutnum  
    waypoint=createallwaypoint(rboundary)
    rentrance=Entrance.rotate(angle/180*math.pi)
    rexit=Exit.rotate(angle/180*math.pi)
    time=timeconsume(waypoint,rentrance,rexit)
    if(time<mintime):
        mintime=time
        optimalangle=angle
rboundary=boundary.rotate(optimalangle/180*math.pi)
bounds=rboundary.bounds
cutnum=float(math.ceil(round((round(bounds[2],6)-round(bounds[0],6))/inputwidth,6)))
width=(round(bounds[2],6)-round(bounds[0],6))/cutnum  
waypoints=createallwaypoint(rboundary)

rwaypoints=[]
for i in range(len(waypoints)):
    rwaypoints.append(waypoints[i].rotate(-optimalangle/180*math.pi))
    
fig=plt.figure()
ax = fig.add_subplot(1,1,1)
pointx,pointy=pointtolist(rwaypoints) 
ax.plot(pointx,pointy)
vertices=boundary.vertices
vertices.append(vertices[0])
verticex,verticey=pointtolist(vertices)
ax.plot(verticex,verticey)

for i in range(len(rwaypoints)-1):
    x=[]
    y=[]
    line=Line(rwaypoints[i],rwaypoints[i+1])
    diracp=line.direction.unit
    diracn=diracp.rotate(math.pi/2)
    rectpt1=rwaypoints[i].translate((diracn.x-diracp.x)*inputwidth/2,(diracn.y-diracp.y)*inputwidth/2)
    rectpt2=rwaypoints[i+1].translate((diracn.x+diracp.x)*inputwidth/2,(diracn.y+diracp.y)*inputwidth/2)
    rectpt3=rwaypoints[i+1].translate((diracp.x-diracn.x)*inputwidth/2,(diracp.y-diracn.y)*inputwidth/2)
    rectpt4=rwaypoints[i].translate((-diracp.x-diracn.x)*inputwidth/2,(-diracp.y-diracn.y)*inputwidth/2)
    x.extend([rectpt1.x,rectpt2.x,rectpt3.x,rectpt4.x])
    y.extend([rectpt1.y,rectpt2.y,rectpt3.y,rectpt4.y])
    ax.fill(x,y,'b',alpha=0.2)
