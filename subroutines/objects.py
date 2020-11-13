#this module contains utilities to read the user-specified planet data file,
#and will allow the user to create custom objects like probes and planetoids, etc.
from math import floor

class probe:
      kind = 'spacecraft'
      def __init__(self,launch,C3,host):
            self.launch = launch #date of launch
            position = host.at(launch).position.au #x,y,z of host at launch date
            self.C3 = C3
            self.x = position[0]
            self.y = position[1]
            self.z = position[2]
            self.vel = [None,None,None]

def get_planet_details(dayskip): 
      #Read in the planet data 
      #Set up the planet dictionaries keyed to the letter of each planet
      plist = {}; pcol = {}; plen = {}; psym = {}; rp = {}; mp = {}
      with open('Data/planet_data.txt','r') as pfile:
            for line in pfile:
                  if line[0]=='#':
                        continue
                  n,d,r,a,u,b,v,rr,i,rc,ic,spk1,spk2,c,s,l = line.split()
                  if spk2 != '#':
                        spk = spk1+' '+spk2
                  else:
                        spk = spk1
                  plist[l] = spk
                  pcol[l] = c
                  plen[l] = max(1,floor(120./dayskip*float(d)**1.5))
                  psym[l] = '$'+str(s)+'$'
                  rp[l] = float(r)
                  mp[l] = {'u':u,'v':v,'b':b,'r':rr,'i':i,'rc':rc,'ic':ic}
      return plist,pcol,plen,psym,rp,mp

#pfile = open('Data/planet_data.txt','r')
