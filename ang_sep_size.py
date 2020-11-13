#options include skyfield(using SPKs from NAIF), spiceyPy, raw HORIZONS data, CALCEPH/INPOP
#migrate to poliastro? It also has access to SPK from NAIF, and can compute orbital elements

from subroutines.set_course import read_course,Course
from subroutines.objects import get_planet_details
from subroutines.functs import to_pol,to_car
from subroutines.functs import get_angle,get_size
from subroutines.probe import find_C3,plot_probe,GMsun
from subroutines.make_gif import fig2data,fig2img
from subroutines.camera import project

import calendar
import time
import numpy as np
from math import floor,degrees,acos,atan,asin
from shlex import split
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from skyfield.api import load,Loader
from skyfield.positionlib import ICRF
#for using the FuncAnimation example at file end
from matplotlib.animation import FuncAnimation, PillowWriter
from PIL import Image

#from scipy.interpolate import interp1d

plotstill = True
plot = False
plot_gif = False
animate = False
animated = False
qdate = []; qxyz = []
dayskip = 1

#C3 math: solar escape velocity = 42.12 km/s
#         earth orbital velocity = 29.78 km/s
#         velocity needed to get free of Earth (11.2 km/s) = 12.7 km/s
#         velocity needed to get free of Sun   (42.12-29.78-12.7) = 6.72 km/s 
# Total delta-v sunk into launch: 12.7+6.72 = 19.42 km/s

load = Loader('./skyfield_data')
ts = load.timescale()
t = ts.utc(2028, 1, range(0,18000,dayskip))
dates = t.utc_datetime()

plist,pcol,plen,psym,rp,mp = get_planet_details(dayskip)

#also need Quaoar: getting from
#https://naif.jpl.nasa.gov/pub/naif/pds/data/nh-j_p_ss-spice-6-v1.0/nhsp_1000/data/spk/kbo_centaur_20170422.bsp
#kbosdata = load('kbo_centaur_20170422.bsp')
#print(kbosdata)
#Approximate location of quaoar at intercept for trajectory1
with open('Data/quaoar.txt','r') as qfile:
      lines = qfile.readlines()[64:1202]
dates = lines[::2]
coord = lines[1::2]
qx = [0 for i in range(len(coord))]
qy = [0 for i in range(len(coord))]
qz = [0 for i in range(len(coord))]
for i in range(len(dates)):
      d1 = dates[i].split()[3].split('-')
      qdate.append(ts.utc(int(d1[0]),list(calendar.month_abbr).index(d1[1]),int(d1[2])))
      qxyz.append(coord[i].split())
      qx[i],qy[i],qz[i] = [float(j) for j in qxyz[-1]]
#quaoar = ICRF(qxyz,t=qdate)

solarsystem = {}; planetsonly = {}
planetdata = load('de438t.bsp')
for p in plist:
      solarsystem[p] = planetdata[plist[p]]
      if p != 'S':
            planetsonly[p] = planetdata[plist[p]]

#This works!!!
crs = read_course(1,planetsonly,ts)
if 'q' in crs.places:
      qind = [0 for i in range(len(t))]
      for j in range(len(t)):
            qind[j] = min(range(len(qdate)), key=lambda i: abs(t[j]-qdate[i]))
#print(crs.stops)
#print(crs.places[-1],crs.coords[-1],crs.at(crs.times[-1]))
cname = 'Course '+'-'.join(crs.places)+crs.name

#Match launch date to a spot in the orrery's time range
ind = {}
for k in crs.places:
      j = crs.places.index(k)
      tj = crs.times[j]
      ind[k] = min(range(len(t)), key=lambda i: abs(t[i]-tj))
sind = ind[crs.places[0]]-floor(365/dayskip)

#build probe trajectory:
di = [[(crs.coords[j+1][i] - crs.coords[j][i])/float(ind[crs.places[j+1]] - ind[crs.places[j]]) for i in range(len(crs.coords[0]))] for j in range(crs.len-1)]
pos = []; pos2 = []; vel = []; tp = []; dprobe = []; px = []; py = []; pz = []
tt = [ct.utc_datetime() for ct in crs.times]
vsep = [crs.coords[j+1] - crs.coords[j] for j in range(crs.len-1)]
dsep = [np.sqrt(x.dot(x)) for x in vsep]
C3 = [None for i in range(crs.len)]
for i in range(crs.len-1):
      delta = (crs.times[i+1].utc_datetime() - crs.times[i].utc_datetime()).days
      c1 = crs.coords[i]; c2 = crs.coords[i+1]
      C3temp,vtemp = find_C3(c1,c2,delta,dayskip)
      xtemp,ytemp,ztemp = plot_probe(c1,c2,vtemp,delta,dayskip)
      #print(crs.places[i]+' - '+crs.places[i+1],C3temp)
      #print(ind)
      for j in range(len(vtemp)):
            #vx,vy,vz = to_car(vtemp[j],0,0)
            pos2.append([xtemp[j],ytemp[j],ztemp[j]])
      pos2.append([c2[0],c2[1],c2[2]])
      #print(pos2[-3:])
#      fig5 = plt.figure()
#      plt.plot(xtemp,ytemp)
#      plt.show()
      
frad= 0. ; ftheta = 0. ; fphi = 0.
new_to_town = True
for i in range(len(t)):
      tutc = t[i].utc_datetime()
      tp.append(t[i])
      if tutc<tt[0]:
            pos.append(solarsystem[crs.places[0]].at(t[i]).position.au)
            vel.append(solarsystem[crs.places[0]].at(t[i]).velocity)
      elif tt[0]<=tutc and tutc<tt[-1]:
            pos.append(np.array(pos2[i-ind[crs.places[0]]]))
            vel.append(vel[-1])
      #for j in range(crs.len-1):
      #      if tt[j]<=tutc and tutc<tt[j+1]:
      #            pos.append(crs.coords[j]+[x*float(i-ind[crs.places[j]]) for x in di[j]])
      #            #vel = ?
      elif tt[-1]<=tutc:
            if new_to_town:
                  new_to_town = False
                  dx,dy,dz = pos[-1] - pos[-2]
            frad,ftheta,fphi = to_pol(pos[-1][0],pos[-1][1],pos[-1][2])
            vrad = np.sqrt(C3temp + 2.*GMsun/frad)
            frad += vrad*dayskip
            #ftheta += dftheta
            #fphi += dfphi
            x,y,z = to_car(frad,ftheta,fphi)
            pos.append(np.array([x,y,z]))#pos[-1]+np.array([dx,dy,dz]))#pos[-1]+[x*float(i-ind[crs.places[-1]]) for x in di[-1]])
            vel.append(vel[-1])
            #vel = ?
      px.append(pos[-1][0])
      py.append(pos[-1][1])
      pz.append(pos[-1][2])
      dprobe.append(np.sqrt(pos[-1].dot(pos[-1])))
##This also now works!
ISP = ICRF(pos,vel,t=tp,center=0)
#print(ISP.t[ind[m]])
#ICRF positions don't support slicing
#print(ISP.position.au[ind[m]])

#Need to think about viewing - angular size and separation from the Sun
sizes = {p:[None]*len(t) for p in planetsonly} #angular size in arcseconds 
angles = {p:[None]*len(t) for p in planetsonly} #angle in degrees between sun and planet p
mag = {p:[None]*len(t) for p in planetsonly} #apparent V-band magnitude of planet p
camera = {p:[x[:] for x in [[None]*3]*len(t)] for p in planetsonly} #angle in degrees between sun and planet p
#out_of_frame = {p:[None]*len(t) for p in planetsonly} #records size when planet p is more than 90 degrees away from Sun
#in_frame = {p:[None]*len(t) for p in planetsonly} #records size when planet p is less than 90 degrees away from Sun
#alpha = {p:[None]*len(t) for p in planetsonly} #angle in degrees between sun and observer from planet p
#include Quaoar
if 'q' in plist:
      sizes['q'] = [None]*len(t)
      angles['q'] = [None]*len(t)
for j in range(len(t[ind[crs.places[0]]+1:])):
      i = j + ind[crs.places[0]]+1
      ISPt = ISP.position.au[i]
      psun = solarsystem['S'].at(t[i]).distance().au
      for p in planetsonly:
            pplanet = planetsonly[p].at(t[i]).position.au
            angles[p][i] = 3600*get_angle(rp[p],ISPt,psun,pplanet)
            #angles[p][i] = 3600*get_angle(rp[p],ISPt,psun,pplanet)
            sizes[p][i] = 3600*get_size(rp[p],pplanet,ISPt)
            camera[p][i] = project(ISPt,ISPt,pplanet)
            #mag[p][i] = mp['v'] #V-band mag data read in 
#fig4 = plt.figure()
#k = 500
#x1=ISP.position.au[ind[crs.places[0]]+k][0]
#y1=ISP.position.au[ind[crs.places[0]]+k][1]
#sx = solarsystem['S'].at(t[ind[crs.places[0]]+k]).position.au[0]
#sy = solarsystem['S'].at(t[ind[crs.places[0]]+k]).position.au[1]
#plt.plot(sx,sy,c='y',marker='o')
#for p in planetsonly:
#      x=planetsonly[p].at(t[ind[crs.places[0]]+k]).position.au[0]
#      y=planetsonly[p].at(t[ind[crs.places[0]]+k]).position.au[1]
#      plt.plot(x,y,c=pcol[p],marker='o')
#      plt.plot(x1,y1,marker='o')
#      plt.plot([x, x1],[y, y1],c=pcol[p])
#      plt.text(x,y,"%4.1f" %(angles[p][ind[crs.places[0]]+k]/3600.))
#plt.show()
#exit()
tele_d = 1. #telescope D in cm
fig2,ax2 = plt.subplots(2,1,sharex=True)
#plt.suptitle(cname)
for p in planetsonly:
      #ax2[0].loglog(dprobe,in_frame[p][:],color=pcol[p])
      ax2[0].loglog(dprobe[ind[crs.places[0]]+1:],sizes[p][ind[crs.places[0]]+1:],color=pcol[p])
      ax2[0].loglog(dprobe[ind[crs.places[0]]+1],sizes[p][ind[crs.places[0]]+1],color=pcol[p],marker='$'+p+'$')
      ax2[0].set_ylabel('Ang. size (arcsec)')
      ax2[1].loglog(dprobe[ind[crs.places[0]]+1:],angles[p][ind[crs.places[0]]+1:],color=pcol[p])
      ax2[1].loglog(dprobe[ind[crs.places[0]]+1],angles[p][ind[crs.places[0]]+1],color=pcol[p],marker='$'+p+'$')
      ax2[1].loglog([min(dprobe[ind[crs.places[0]]+1:]),max(dprobe[ind[crs.places[0]]+1:])],[4*1.22*206265*2.5E-4/tele_d, 4*1.22*206265*2.5E-4/tele_d],'k:')
      ax2[1].text(min(dprobe[ind[crs.places[0]]+1:]),1.22*206265*2.5E-4/1,'4$\lambda$/D @ 2.5 um, D='+str(tele_d)+' cm')
      ax2[1].set_xlabel('Barycentric distance of ISP (au)')
      ax2[1].set_ylabel('Ang. sep. Sun-planet (arcsec)')
      ax2[1].set_ylim(206265.*2.5E-4/tele_d,3600.*90.)
      #if plot_date:
      #      ax2[1].set_xlabel('Date (UTC)')
fig2.savefig('distance_'+cname+'2.png',bbox_inches='tight',dpi=300)
plt.show()
#remember, diffraction limit is ~206265*n*lambda/D (arcseconds); n is usually on the order of a few

#I want to make a movie of what the RVC would see - map positions of planets onto 
# the camera plane.
max_fov = 90. # in degrees
max_fov *= 3600. #To arcseconds

if plotstill:
      fig = plt.figure(figsize=(7,7),dpi=100)
      if crs.name in ['1', '2']:
            ax = plt.axes(xlim=(-10, 20), ylim=(-50, 10))
      elif crs.name == '3':
            ax = plt.axes(xlim=(-20, 100), ylim=(-10, 50))
      ax.set_aspect('equal','box')
      #ax.set_title(cname)
      start = ind[crs.places[0]]; end = ind[crs.places[-1]]
      ax.plot(px[start:end],py[start:end],'k',lw=2)
      for j in range(crs.len):
            ax.plot(crs.coords[j][0],crs.coords[j][1],c='gray',marker='x')
      for p in plist:
            sweep = max(0,end-plen[p])
            x,y,z = solarsystem[p].at(t[sweep:end]).position.au
            ax.plot(x,y,pcol[p],lw=1)
            x,y,z = solarsystem[p].at(t[end]).position.au
            ax.plot(x,y,pcol[p],marker='$'+p+'$')
      if 'q' in crs.places:
            #ax.plot(qx[:qind[end]],qy[:qind[end]],'pink',lw=1)
            ax.plot(qx[qind[start]:qind[end]],qy[qind[start]:qind[end]],'pink',lw=1)
            ax.plot(qx[qind[end]],qy[qind[end]],'pink',marker='$q$')
      if 'e' in crs.places:
            ax.plot(crs.coords[-1][0],crs.coords[-1][1],'gray',marker='$e$')
      fig.savefig('still_'+cname+'.png',bbox_inches='tight',dpi=300)
      plt.show()

if plot:
      images = []
      unplot = [True for i in range(crs.len)]
      launch = {}
      #Make moving figure of orbits, trajectory
      line = {}; dot ={}
      if crs.name in ['1', '2']:
            fig3 = plt.figure(figsize=(5,5),dpi=150)
            ax = plt.axes(xlim=(-50, 50), ylim=(-50, 50))
      elif crs.name == '3':
            fig3 = plt.figure(figsize=(12.5,7.5),dpi=150)
            ax = plt.axes(xlim=(-25, 100), ylim=(-50, 50))
            #ax = plt.axes(xlim=(-30, 30), ylim=(-30, 30))
      ax.set_aspect('equal')
      ax.set_title(cname)#'Course '+str([c+'-' for c in crs.places]))
      end = ind[crs.places[-1]]+10 #len(t)#-sind
      for j in range(crs.len):
            launch[j], = ax.plot(crs.coords[j][0],crs.coords[j][1],c='gray',marker='o')
      for p in plist:
            dot[p], = ax.plot([],[],pcol[p],marker='$'+p+'$') #psym[p])
            line[p], = ax.plot([],[],pcol[p],lw=1)
      if 'q' in crs.places:
            qua, = ax.plot([],[],'pink',lw=1)
            quadot, = ax.plot([],[],'pink',marker='$q$')
      timetext = plt.text(0.12, 0.95,t[sind-1].utc_strftime('%m-%d-%Y'), ha='center', va='center', transform=ax.transAxes)
      for i in range(sind,end,7):
            timetext.set_text(t[i].utc_strftime('%m-%d-%Y'))
            if i>ind[crs.places[0]]+1:
                  ax.plot(px[ind[crs.places[0]]:i],py[ind[crs.places[0]]:i],'k',lw=2)
                  #maxISP = max(5,1.5*dprobe[i])
                  #ax.set_xlim=(-1*maxISP, maxISP)
                  #ax.set_ylim=(-1*maxISP, maxISP)
                  #plt.axes(xlim=(-1*maxISP, maxISP),ylim=(-1*maxISP, maxISP))
            for p in plist:
                  sweep = max(i-plen[p],0); swept = i
                  x,y,z = solarsystem[p].at(t[sweep:swept]).position.au
                  line[p].set_ydata(y)
                  line[p].set_xdata(x)
                  x,y,z = solarsystem[p].at(t[i]).position.au
                  dot[p].set_ydata(y)
                  dot[p].set_xdata(x)
            if 'q' in crs.places:
                  qua.set_ydata(qy[:qind[i]])
                  qua.set_xdata(qx[:qind[i]])
                  quadot.set_ydata(qy[qind[i]])
                  quadot.set_xdata(qx[qind[i]])
            for c in crs.times:
                  j = crs.times.index(c)
                  if c.utc_datetime()<=t[i].utc_datetime() and unplot[j]:
                        launch[j].set_color('b')
                        unplot[j] = False
            plt.pause(0.00005)
            if plot_gif:
                  im = fig2img(fig3)
                  images.append(im)
      if plot_gif:
            images[0].save(cname+'.gif',
                  save_all=True, append_images=images[1:], optimize=False, duration=10, loop=0)
if animate:
      #gif example
      i=0
      #line = []
      fig = plt.figure()
      ax = plt.axes(xlim=(-40, 40), ylim=(-40, 40))
      line, = ax.plot([],[],pcol['S'],lw=2)
      line2, = ax.plot([],[],pcol['E'],lw=2)
      line3, = ax.plot([],[],pcol['J'],lw=2)
      line4, = ax.plot([],[],pcol['N'],lw=2)
      #for p in plist:
      #      test,  = ax.plot([],[],pcol[p],lw=2)
      #      line.append(test)

      def init():
            #j=0
            #for p in plist:
            #      line[j].set_data([],[])
            #      j+=1
            line.set_data([], [])#,[],[],[],[],[],[])
            line2.set_data([], [])#,[],[],[],[],[],[])
            line3.set_data([], [])#,[],[],[],[],[],[])
            line4.set_data([], [])#,[],[],[],[],[],[])
            return line,line2,line3,line4,
      def animate(i):
            #j = 0
            #for p in plist:
            #      x,y,z = solarsystem[p].at(t[sind+i-plen[p]:sind+i]).position.au
            #      line[j].set_data(x,y)
            #      j+=1
            x,y,z = solarsystem['S'].at(t[sind-5:sind+i]).position.au
            x2,y2,z2 = planetsonly['E'].at(t[sind-plen['E']:sind+i]).position.au
            x3,y3,z3 = planetsonly['J'].at(t[sind-plen['J']:sind+i]).position.au
            x4,y4,z4 = planetsonly['N'].at(t[sind-plen['N']:sind+i]).position.au
            line.set_data(x, y)
            line2.set_data(x2,y2)
            line3.set_data(x3,y3)
            line4.set_data(x4,y4)
            return line,line2,line3,line4,

      anim = FuncAnimation(fig, animate, init_func=init,
                               frames=500, interval=20, blit=True)

      anim.save('orbits.gif', writer=PillowWriter())
      plt.show()

if animated:
      unplot = [True for i in range(crs.len)]
      launch = {}
      #Make moving figure of orbits, trajectory
      line = {}; dot ={}
      fig = plt.figure(figsize=(7,7),dpi=100)
      ax = plt.axes(xlim=(-5, 5), ylim=(-5, 5))
      ax.set_aspect('equal','box')
      ax.set_title(cname)#'Course '+str([c+'-' for c in crs.places]))
      end = len(t)#-sind
      def init():
            for j in range(crs.len):
                  launch[j], = ax.plot(crs.coords[j][0],crs.coords[j][1],c='gray',marker='o')
            for p in plist:
                  dot[p], = ax.plot([],[],pcol[p],marker='$'+p+'$') #psym[p])
                  line[p], = ax.plot([],[],pcol[p],lw=1)
                  qua, = ax.plot([],[],'pink',lw=1)
            timetext = plt.text(0.1, 0.95,t[sind-1].utc_strftime('%m-%d-%Y'), ha='center', va='center', transform=ax.transAxes)
      def animate(i):
      #for i in range(sind,end,14):
            timetext.set_text(t[i].utc_strftime('%m-%d-%Y'))
            #ax.plot(1,1,'+')
            if i>ind[crs.places[0]]+1:
                  ax.plot(px[ind[crs.places[0]]:i],py[ind[crs.places[0]]:i],'k',lw=2)
                  #maxISP = max(5,1.5*dprobe[i])
                  #ax.set_xlim=(-1*maxISP, maxISP)
                  #ax.set_ylim=(-1*maxISP, maxISP)
                  #plt.axes(xlim=(-1*maxISP, maxISP),ylim=(-1*maxISP, maxISP))
            for p in plist:
                  sweep = max(i-plen[p],0); swept = i
                  x,y,z = solarsystem[p].at(t[sweep:swept]).position.au
                  line[p].set_ydata(y)
                  line[p].set_xdata(x)
                  x,y,z = solarsystem[p].at(t[i]).position.au
                  dot[p].set_ydata(y)
                  dot[p].set_xdata(x)
            qua.set_ydata(qy[:qind[i]])
            qua.set_xdata(qx[:qind[i]])
            for c in crs.times:
                  j = crs.times.index(c)
                  if c.utc_datetime()<=t[i].utc_datetime() and unplot[j]:
                        launch[j].set_color('b')
                        unplot[j] = False
            #plt.pause(0.00005)
            return line,dot,qua,
      anim = FuncAnimation(fig, animate, init_func=init,
                               frames=len(t), interval=10, blit=True)

      anim.save(cname+'.gif', writer=PillowWriter())
      
