from math import atan,atan2,acos,log10
from numpy import sqrt,pi,logspace
import sys
from subroutines.functs import to_car,to_pol

GMsun = 3.986E14 #<Earth #1.327E20 #<Sun #m3/s2
m3s2toau3day2 = (24.*3600.)**2.*(1./1.496E11)**3. #convert from m3/s2 to au3/day2
GMsun = GMsun * m3s2toau3day2

#Bisecting root finder for C3', since C3 is a function of angle as well
def find_C3(c1,c2,ddays,dayskip):
      C3bot = -1.E-5; C3top = 3000.
      arrival_right = False; tries = 0
      vrad = [None for i in range(ddays+1)]
      vrad2 = [None for i in range(ddays+1,dayskip)]

      x1,y1,z1 = c1
      x2,y2,z2 = c2
      r1,theta1,phi1 = to_pol(x1,y1,z1)
      r2,theta2,phi2 = to_pol(x2,y2,z2)
      while not arrival_right:
            rguess = r1
            C3guess = (C3bot+C3top)/2.
            vrad[0] = sqrt(C3guess + 2*GMsun/rguess)
            for j in range(0,ddays):
                  rguess = rguess+vrad[j]*1. #vrad is in au/day, rguess in au
                  vrad[j+1] = sqrt(C3guess + 2*GMsun/rguess)
            #print(r2,rguess)
            if tries> 100:
                  print(C3guess,rguess,r2)
                  sys.exit('No root found. Consider modifying find_C3 bounds.')
            if abs(rguess-r2)<1E-6:
                  arrival_right = True
            elif rguess<r2:
                  tries+=1
                  C3bot = C3guess
            else:
                  tries+=1
                  C3top = C3guess
            vrad2 = vrad[::dayskip]
      return C3guess,vrad2


def plot_probe(c1,c2,v,ddays,dayskip):
      #0=uniform dtheta, dphi sweep with r [big looping path]; 1=uniform dx,dy,dz sweep [smoother, hooks in before J]; 
      #2=straight line w/ r+dr constraint [skips over path interior to E]; 3=parabolic [broken];
      #4=min. rad?
      method = -1 #2 is optimal at this point
      #if ddays != len(v):
      #      sys.exit('Days in path do not divide evenly into velocity vector.') 
      x = [None for i in range(len(v))]
      y = [None for i in range(len(v))]
      z = [None for i in range(len(v))]
      dd1 = ddays 
      x[0],y[0],z[0] = c1
      r,theta,phi = to_pol(x[0],y[0],z[0])
      x2,y2,z2 = c2
      r2,theta2,phi2 = to_pol(x2,y2,z2)
      if theta2 - theta < -1*pi:
            dtheta = (2.*pi + theta2 - theta)/dd1
      elif theta2 - theta > pi:
            dtheta = (theta2 - theta - 2.*pi)/dd1
      else:
            dtheta = (theta2 - theta)/dd1
      dphi = (phi2 - phi)/dd1
      signdtheta = dtheta/abs(dtheta)
      signt = theta/abs(theta); signt2 = theta2/abs(theta2)
      dt2 = logspace(log10(theta/signt),log10(theta2/signt2),num=len(v))
      dx = x2 - x[0]; dy = y2 - y[0]; dz = z2 - z[0]
      dnx = dx/dd1
      dny = dy/dd1
      dnz = dz/dd1
      if method == 2:
            a = (x2 - x[0])**2. + (y2 - y[0])**2. + (z2 - z[0])**2.
            b = 2*(dx*x[0] + dy*y[0] + dz*z[0])
      elif method == 3:
            a = 1.
            b = dx**2./dz + dy**2./dz
      #b2a2 = (y2**2. - y[0]**2.)/(x2**2. - x[0]**2.)
      #print(b2a2,sqrt(1+b2a2))
      r0 = r
      #print(c2,[x[0]+dnx*ddays,y[0]+dny*ddays,z[0]+dnz*ddays])
      for i in range(0,len(v)-1):
            #thetanew = atan2(y[i]+dny*dayskip,x[i]+dnx*dayskip)
            #dtheta2 = (thetanew - theta)
            dr = v[i]*1.*dayskip
            rold = r
            r += dr
            if method == -1:
                  x[i+1] = x[i] + dnx*dayskip
                  y[i+1] = y[i] + dny*dayskip
                  z[i+1] = z[i] + dnz*dayskip
            if method == 0:
                  phi += dphi*dayskip
                  theta += dtheta*dayskip
                  x[i+1],y[i+1],z[i+1] = to_car(r,theta,phi)
                  continue
            if method == 4:
                  phi += dphi+dayskip
                  theta += signdtheta*(dt2[-i-1] - dt2[-i-2])
                  x[i+1],y[i+1],z[i+1] = to_car(r,theta,phi)
                  continue
            if method == 1:
                  phi = acos((z[i]+dnz*dayskip)/r)
                  theta = atan2(y[0]+(i+1)*dny*dayskip,x[0]+(i+1)*dnx*dayskip)
                  x[i+1],y[i+1],z[i+1] = to_car(r,theta,phi)
                  continue
            if method in [2,3]:
                  if method == 2:
                        c = r0**2. - r**2.
                  elif method == 3:
                        c = -1.*b*z[0] - r**2.
                  rad = b**2. - 4*a*c
                  if rad < 0:
                        sys.exit('No root found. Revise plot_probe!')
                  elif rad == 0:
                        t = -b/(2*a)
                  else:
                        t1 = (-b+sqrt(rad))/(2*a)
                        t2 = (-b-sqrt(rad))/(2*a)
                        t = max(t1,t2)
                  if method == 2:
                        x[i+1],y[i+1],z[i+1] = [x[0]+dx*t,y[0]+dy*t,z[0]+dz*t]
                  elif method == 3:
                        x[i+1],y[i+1],z[i+1] = [x[0]+dx*sqrt(t-z[0])/sqrt(dz),y[0]+dy*sqrt(t-z[0])/sqrt(dz),t]
      r3,theta3,phi3 = to_pol(x[-1],y[-1],z[-1])
      if abs(r2-r3) > 1.E-5:
            print(r2,r3,'...',theta2,theta3,'...',phi2,phi3)
            sys.exit('New radius in plot_probe does not match destination. '+str(r3)+' '+str(r2))
      return x,y,z

