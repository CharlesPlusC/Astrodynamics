import numpy as np
import pandas as pd
import matplotlib as plt
import math


 # Converted Cartesian
 #          x         -7134.401598
 #          y         -1344.205350
 #          z          2616.199171
 #          u             2.737023
 #          v            -2.641275
 #          w             6.099438
          
          
class ECI_ECEF: 
    
    #Defining general variables for all classes
    GM = 398600.4415 #Km^3 s^-2
    Re = 6378.137 #Km
    w_cross = (2*(np.pi) / 86400) #Rad/Sec
    
    def __init__(self, 
                     x_i,y_i,z_i,u_i,v_i,w_i,
                     x_f,y_f,z_f,u_f,v_f,w_f,
                 ):
        
        self.x_i = x_i
        self.y_i = y_i
        self.z_i = z_i
        self.u_i = u_i
        self.v_i = v_i
        self.w_i = w_i
        
        self.x_f = x_f
        self.y_f = y_f
        self.z_f = z_f
        self.u_f = u_f
        self.v_f = v_f
        self.w_f = w_f
        
        
    def ECI2ECEF(self,d):
        
        OGAST = 280.4606+360.9856473662*d
        #Making the coordinates into the relevant vectors

        X_i = np.array([x_i,y_i,z_i])
        v_i = np.array([u_i,v_i,w_i])
        
        #Constructing a rotation matrix for position vectors
        R_z = ([[np.cos(OGAST),  np.sin(OGAST),  0,
        -np.sin(OGAST), np.cos(OGAST),  0,
        0,              0,              1]])
        
        #Constructing a rotation matrix for velocity vectors
        mat_a = ([[np.sin(OGAST),  -np.cos(OGAST),  0,
            np.cos(OGAST), np.sin(OGAST),   0,
            0,              0,              1]]) 

        R_dot_z = w_cross.dot(mat_a)
        
        #Converted position coordinates from eci to ecef = multiply coordinates by rotation matrix
        X_f = R_z.dot(X_i)
        
        #The differenetial of Xf (i.e. v_f) = (R_dot_z *X_i) + (R_z *V_i) --Product rule
        v_f = R_dot_z.dot(X_i) + R_z.dot(V_i)
        
        ECI_2_ECEF = pd.DataFrame ({"Coordinate":["x","y","z","u","v","w"],
                            "ECI": [x.i, y.i, z.i, u.i, v.i, w.i],
                            "ECEF":[x.f, y.f, z.f, u.f, v.f, w.f]})

        print ('Table of results:')
        ECI_2_ECEF
        
    # def ECEF2ECI(self):
        