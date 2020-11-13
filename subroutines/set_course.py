#This file contains the functions necessary to 'set course',
# translating between calendar dates and solar system locations
import os
import sys
from subroutines.functs import to_pol,to_car

class Course:
      kind = 'stops'
      def __init__(self,number):
            self.name = str(number)
            self.len = 0
            self.times = []
            self.places = []
            self.coords = []
            self.stops = {}
            #self.att = {}
      def add_stop(self,dest,time,place):
            if(time not in self.times):
                  self.stops[dest] = [time, place] 
                  self.times.append(time)
                  self.places.append(dest)
                  self.coords.append(place)
                  #self.att[time] = place
                  self.len += 1
            else:
                  sys.exit("Time value already exists in Course"+self.name)
      def at(self,time):
            t = self.times.index(time)
            return self.coords[t]

def read_course(course_number,planets,time):
      course = Course(course_number)
      cc = [None,None,None]
      
      f = 'Courses/'+str(course_number)+'.txt'
      if not os.path.isfile(f):
            sys.exit("Course file not found.")
      with open(f,'r') as cfile:
            for line in cfile:
                  custom_coord = False
                  if line[0]=='#':
                        continue
                  if len(line.split())==4:
                        p,x1,x2,x3 = line.split()
                  elif len(line.split())==7:
                        p,x1,x2,x3,c1,c2,c3 = line.split()
                        cc = [float(i) for i in [c1,c2,c3]]
                        custom_coord = True
                  else:
                        sys.exit("Course data in"+f+"not compatible.")
                  y,m,d = [int(i) for i in [x1,x2,x3]]
                  tp = time.utc(y,m,d)
                  if p not in planets:
                        if not custom_coord:
                              sys.exit("Planet not in data and custom coord. not provided.")
                        course.add_stop(p, tp, cc)
                  else:
                        course.add_stop(p, tp, planets[p].at(tp).position.au)
      return course
