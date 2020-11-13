#This module contains useful functions to convert between coordinates
# and how to get angular sizes and diameters

from math import degrees,atan,atan2,acos,sin,cos,inf
from numpy import sqrt,dot
import sys

def to_pol(x,y,z):
      r = sqrt(x**2.+y**2.+z**2.)
      theta = atan2(y,x)  #atan(y/x)
      phi = acos(z/r)
      return r,theta,phi

def to_car(r,theta,phi):
      x = r*sin(phi)*cos(theta)
      y = r*sin(phi)*sin(theta)
      z = r*cos(phi)
      return x,y,z
            
#get_size(rp[p],pplanet,ISPt)
def get_size(r,p1,p2):
      s = 0.
      v = p1 - p2
      d = sqrt(v.dot(v))
      if d>0.:
            s = 2.*degrees(atan(r/(1.496E11*d)))
      else:
            s = inf
      return s

#get_angle(rp[p],ISPt,psun,pplanet)
# dos = 1 - 2; dbo = 3 - 1; dbs = 3 - 2
#note that this angle is p3-p1-p2 (vertex is located at p1)
def get_angle(r,p1,p2,p3):
      a = 0.
      rau = r/1.496E11 #convert r from km to au
      v12 = p1 - p2
      d12 = sqrt(v12.dot(v12))
      v23 = p2 - p3
      d23 = sqrt(v23.dot(v23))
      v31 = p3 - p1
      d31 = sqrt(v31.dot(v31)) #max(rau,sqrt(v31.dot(v31)))
      anglearg = (d12**2. + d31**2. - d23**2.)/(2.*d12*d31)
      if -1.<anglearg<1.:
            a = degrees(acos(anglearg))
      else:
            a = inf
      return a

def get_magnitude(planet,r,p1,p2,p3):
      angle = get_angle(r,p2,p1,p3) #flip arguments to get alpha w/ vertex at planet
      
      return mag

