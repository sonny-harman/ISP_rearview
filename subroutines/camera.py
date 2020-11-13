from numpy import array,sqrt,dot

def project(loc,sun,target):
      #project takes the current position location, normal to the image plane
      # (dictated either by the normal determined by subtracting the current loc from 
      # the previous loc) or by the ray connecting the probe to the Sun. Target coords
      # are then translated onto the camera plane. We could then add a pixel transform
      # if it's useful.
      b = [0.,0.,0.]
      #z axis is normal
      normal = loc - sun
      unit_normal = [[n / sqrt(normal.dot(normal))] for n in normal]
      #assume that the view window's x axis is parallel to barycentric x axis
      unit_w = [1.,0,0]
      unit_x = [1.,0,0]
      #compute y axis as perpendicular to unit x and unit normal
      unit_y = [0,1,0]
      return b

