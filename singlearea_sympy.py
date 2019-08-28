# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import Point, Line, Polygon,Segment

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
    xpos=xpos.reshape(len(cutpos),1)
    xpts=np.insert(xpos,1,[0],axis=1)
    xpts=map(Point,xpts)
    xaxis=Line((0,0),(1,0))
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
            i=0
            for side in sides:
                if(side.contains(rightcross)):
                    if(rightcross.y>side.p1.y or(rightcross.y==side.p1.y and rightcross.y<side.p2.y)):
                        rightdownindex=i
                        rightupindex=i+1
                    else:
                        rightupindex=i
                        rightdownindex=i+1
                if(side.contains(leftcross)):
                    if(leftcross.y>side.p1.y or(leftcross.y==side.p1.y and leftcross.y<side.p2.y)):
                        leftdownindex=i
                        leftupindex=i+1
                    else:
                        leftupindex=i
                        leftdownindex=i+1
                i=i+1
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
    
def createwaypoint(up1,down1,pos):
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
def pointtolist(vertices):
    pointx=[]
    pointy=[]
    for vertice in vertices:
        pointx.append(vertice.x)
        pointy.append(vertice.y)
    return pointx,pointy
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
num=8
bptlist=[(0,0),(2,0),(2,2),(0,2),(0,1.5),(1,1.5),(1,0.5),(0,0.5)]

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
#convert to rotated frame
angle=float(input("Enter the angle of the route(degree):"))
angle=angle/180*math.pi
rboundary=boundary.rotate(angle)

#define the width,generate cut position
width=float(input("Enter width of the route:"))
bounds=rboundary.bounds
cutnum=float(math.ceil((bounds[2]-bounds[0])/width))
width=(bounds[2]-bounds[0])/cutnum
xaxis=Line((0,0),(1,0))
yaxis=Line((0,0),(0,1))
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
            waypoint=createwaypoint(up1,down1,pointpos[i])
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
rwaypoints=[]
for i in range(len(waypoints)):
    rwaypoints.append(waypoints[i].rotate(-angle))
    
fig=plt.figure()
ax = fig.add_subplot(1,1,1)
pointx,pointy=pointtolist(rwaypoints) 
ax.plot(pointx,pointy)
vertices=boundary.vertices
vertices.append(vertices[0])
verticex,verticey=pointtolist(vertices)
ax.plot(verticex,verticey)
x=[]
y=[]
for i in range(len(rwaypoints)-1):
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
pivot=np.argwhere(vertex[:,0]==coord[num-1][0])[0][0]2-

#find cross point with boundry
allvertex=vertex
upaddnum=0
for i in range(len(cutpos)-1):
    j=0
    k=pivot
    while(cutpos[i]>vertex[j][0]):
        j=j+1
    lefttop=vertex[j-1]
    righttop=vertex[j]
    while(k<len(vertex)):
        if(cutpos[i]<vertex[k][0]):
            k=k+1
        else:
            rightbuttom=vertex[k-1]
            leftbuttom=vertex[k]
            break
    if(k==len(vertex)):
        rightbuttom=vertex[k-1]
        leftbuttom=vertex[0]
    upperpos=[(cutpos[i],(cutpos[i]-lefttop[0])/(righttop[0]-lefttop[0])*(righttop[1]-lefttop[1])+lefttop[1])]
    lowerpos=[(cutpos[i],(cutpos[i]-leftbuttom[0])/(rightbuttom[0]-leftbuttom[0])*(rightbuttom[1]-leftbuttom[1])+leftbuttom[1])]
    allvertex=np.insert(allvertex, upaddnum+j, values=upperpos, axis=0)
    allvertex=np.insert(allvertex, upaddnum+k+1, values=lowerpos, axis=0)
    upaddnum=upaddnum+1
#find range of waypoints ineach cut
pointrange=np.zeros(shape=[0,2])
waypoint=np.zeros(shape=[0,2])
allvertex=np.append(allvertex,allvertex[0].reshape(1,2),axis=0)
allpivot=np.argwhere(allvertex[:,0]==coord[num-1][0])[0][0]
for i in range(len(cutpos)):
    if(i==0):
        leftupindex=0
        leftdownindex=len(allvertex)-1
    else:
        leftupindex=np.argwhere(allvertex[:,0]>=cutpos[i-1])[0][0]
        leftdownindex=np.argwhere(allvertex[allpivot:len(allvertex),0]<cutpos[i-1])[0][0]+allpivot-1
    if(i==len(cutpos)-1):
        rightupindex=allpivot
        rightdownindex=allpivot
    else:
        rightupindex=np.argwhere(allvertex[:,0]>cutpos[i])[0][0]-1
        rightdownindex=np.argwhere(allvertex[allpivot:len(allvertex),0]<cutpos[i])[0][0]+allpivot-1       
    upsearchrange=allvertex[leftupindex:rightupindex+1,:]
    downsearchrange=allvertex[rightdownindex:leftdownindex+1,:]   
    topmaxindex=upsearchrange.argmax(0)[1]
    topmin=upsearchrange.min(0)
    buttommax=downsearchrange.max(0)
    buttomminindex=downsearchrange.argmin(0)[1]
    if(topmin[1]-width/2>buttommax[1]+width/2):
        newpoint=np.array([topmin[1]-width/2,buttommax[1]+width/2])
    else:
        newpoint=np.array([(topmin[1]+buttommax[1])/2])
    newpoint=newpoint.reshape(len(newpoint),1)    
    newpoint=np.insert(newpoint,0,[pointpos[i]],axis=1)
    if(topmaxindex==0):
        topx1=upsearchrange[topmaxindex][0]
        if(len(upsearchrange==2)):
            topx2=cutpos[i]
        else:
            topx2=findcross(upsearchrange[topmaxindex],upsearchrange[topmaxindex+1],upsearchrange[topmaxindex][1]-width/2)
    elif(topmaxindex==len(upsearchrange)-1):
        if(len(upsearchrange==2)):
            topx1=cutpos[i-1]
        else:
            topx1=findcross(upsearchrange[topmaxindex-1],upsearchrange[topmaxindex],upsearchrange[topmaxindex][1]-width/2)
        topx2=upsearchrange[topmaxindex][0]
    else:
        topx1=findcross(upsearchrange[topmaxindex-1],upsearchrange[topmaxindex],upsearchrange[topmaxindex][1]-width/2)
        topx2=findcross(upsearchrange[topmaxindex],upsearchrange[topmaxindex+1],upsearchrange[topmaxindex][1]-width/2)
    if(topx1<=pointpos[i] and topx2>=pointpos[i]):
        newpoint=np.insert(newpoint,0,[[pointpos[i],upsearchrange[topmaxindex][1]-width/2]],axis=0)
    elif(topx1>pointpos[i]):
        newpoint=np.insert(newpoint,0,[[topx1,upsearchrange[topmaxindex][1]-width/2]],axis=0)
    else:
        newpoint=np.insert(newpoint,0,[[topx2,upsearchrange[topmaxindex][1]-width/2]],axis=0)
        
    if(buttomminindex==len(upsearchrange)-1):
        buttomx1=downsearchrange[buttomminindex][0]
        buttomx2=findcross(downsearchrange[buttomminindex],downsearchrange[buttomminindex-1],downsearchrange[buttomminindex][1]-width/2)
    elif(buttomminindex==0):
        buttomx1=findcross(downsearchrange[buttomminindex+1],downsearchrange[buttomminindex],downsearchrange[buttomminindex][1]-width/2)
        buttomx2=downsearchrange[buttomminindex][0]
    else:
        buttomx1=findcross(downsearchrange[buttomminindex+1],downsearchrange[buttomminindex],downsearchrange[buttomminindex][1]-width/2)
        buttomx2=findcross(downsearchrange[buttomminindex],downsearchrange[buttomminindex-1],downsearchrange[buttomminindex][1]-width/2)
    if(buttomx1<=pointpos[i] and buttomx2>=pointpos[i]):
        newpoint=np.append(newpoint,[[pointpos[i],downsearchrange[buttomminindex][1]-width/2]],axis=0)
    elif(topx1>pointpos[i]):
        newpoint=np.append(newpoint,[[buttomx1,downsearchrange[buttomminindex][1]-width/2]],axis=0)
    else:
        newpoint=np.append(newpoint,[[buttomx2,downsearchrange[buttomminindex][1]-width/2]],axis=0)
    if(i%2==1):
        newpoint=np.flipud(newpoint)
    
    waypoint=np.append(waypoint,newpoint,axis=0)
#rotate back everything
rwaypoint=np.zeros(shape=[len(waypoint),2])
rvertex=np.zeros(shape=[len(vertex),2])
sin=math.sin(-angle/180*math.pi)
cos=math.cos(-angle/180*math.pi)
for i in range(len(waypoint)):
    rwaypoint[i]=rotate(waypoint[i])
for i in range(len(vertex)):
    rvertex[i]=rotate(vertex[i])
rwaypoint=rwaypoint+gcoord[0]
rvertex=rvertex+gcoord[0]
#plot the waypoints
plt.plot(rwaypoint[:,0],rwaypoint[:,1])
poly = geometry.Polygon(rvertex)
x,y = poly.exterior.xy
plt.plot(x,y)
'''