import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
# A set of Keplerian elements define the satellite's state at an instant in time
# They can also be used to predict the subsequent trajectory assuming two-body problem dynamics

###Declaring all required variables###


# n - Mean motion
# Mo - Mean anomaly at start epoch
# Mi - Mean anomaly at required epoch
# a - Semi-major axis
# e - Eccentricity
# i - Inclination
# w - Right ascension of the ascending node
# W - Argument of perigee
# Vo - True anomaly at the start epoch
# r - Initial starting range
# Eo - Eccentric anomaly at the initial epoch
# Ei - Eccentric anomaly at the required epoch
# sinEo, cosEo - variables to hold the cos and sin of eccentric anomaly at the initial epoch
# x, y - Magnitude of the gaussian vectors at (t+dt)
# xdot, ydot - derivative with respect to time of x and y at t+dt
# P,Q - Gaussian vectors
# X_bar - state vector at time t+dt

### Stating the method ###

# 1) Convert inertial state to Keplerian elements
# 2) Compute mean anomaly at start point
# 3) Compute mean anomaly at t+dt
# 4) use M(t+dt) to determine E(t+dt) - Solving Kepler's equation
# 5) Calculate the magnitude of the x,y vector components at (t+dt)
# 6) Project these onto the P and Q vectors as in Kep2Cart: position(t+dt)
# 7) Calculate xdot,ydot, project onto P and Q: velocity (t+dt)


#semi-major axis:7719.637186
#eccentricity:0.000493
#inclination:1.152689
#argument of perigee:1.161349
#RAAN:3.167019
#True Anomaly:5.501897
# Input Kep2Car function here to get keplerian elements

#Example input:
#orbit_1 = OrbitCoordinates(7719.637186,0.000493,1.152689,1.161349,3.167019,5.501897)

class OrbitCoordinates:
    
    def __init__(self,a, e, i, w, W, V):
        self.a = a
        self.e = e
        self.i = i
        self.w = w
        self.W = W
        self.V = V
        
#Example input:
#orbit_1.ReturnInitialOrbit()
        
    def ReturnInitialOrbit(self):
        Specified_Orbit_Geometry = pd.DataFrame({"Orbital Parameter": ["semi-major axis - a", "eccentricity - e",
                                                                       "inclination - i", "argument of perigee - w",
                                                                       "RAAN - W", "TRAN - V"],
                                                 "Keplerian": [self.a, self.e, self.i, self.w, self.W, self.V]})
        print('\nSpecified Orbit Geometry:')
        print(Specified_Orbit_Geometry)  # Printing input data as a check
        
#Example input:
#orbit_1.Kep2Cart()     
        
    def Kep2Cart(self):  # in km and km/s
        GM = 398600.4415
    
        # Compute the in-orbital plane Gaussian Vectors
        # This gives P and Q in ECI components
        P = np.matrix([[np.cos(self.W) * np.cos(self.w) - np.sin(self.W) * np.cos(self.i) * np.sin(self.w)],
                       [np.sin(self.W) * np.cos(self.w) + np.cos(self.W) * np.cos(self.i) * np.sin(self.w)],
                       [np.sin(self.i) * np.sin(self.w)]])
    
        Q = np.matrix([[-np.cos(self.W) * np.sin(self.w) - np.sin(self.W) * np.cos(self.i) * np.cos(self.w)],
                       [-np.sin(self.W) * np.sin(self.w) + np.cos(self.W) * np.cos(self.i) * np.cos(self.w)],
                       [np.sin(self.i) * np.cos(self.w)]])
    
        # Compute the semi-latus rectum and the radial distance
        global p
        p = self.a * (1 - (self.e ** 2))
        global r
        r = p / (1 + self.e * np.cos(self.V))
    
        # x and y are the coordinates of the satellite in the orbital basis
        x = r * np.cos(self.V)
        y = r * np.sin(self.V)
    
        # We know the inertial vector components along the P and Q vectors.
        # Thus we can project the satellite position onto the ECI basis.
        global cart_pos_x
        cart_pos_x = (x * P.item(0)) + (y * Q.item(0))
        global cart_pos_y
        cart_pos_y = (x * P.item(1)) + (y * Q.item(1))
        global cart_pos_z
        cart_pos_z = (x * P.item(2)) + (y * Q.item(2))
    
        # For the velocity components of the state vector we require cosE,sinE,f and g as defined below
    
        cos_E = ((x / self.a) + self.e)
        sin_E = (y / self.a * np.sqrt(1 - self.e ** 2))
    
        f = (np.sqrt(self.a * GM)) / r
        g = np.sqrt(1 - self.e ** 2)
    
        global cart_vel_x
        cart_vel_x = (-f * sin_E * P.item(0)) + (f * g * cos_E * Q.item(0))  # x component of velocity - aka 'u'
        global cart_vel_y
        cart_vel_y = (-f * sin_E * P.item(1)) + (f * g * cos_E * Q.item(1))  # y component of velocity - aka 'v'
        global cart_vel_z
        cart_vel_z = (-f * sin_E * P.item(2)) + (f * g * cos_E * Q.item(2))  # z component of velocity - aka 'w'
    
        Kep2Car = pd.DataFrame({"Orbital Parameter": ["semi-major axis - a", "eccentricity - e", "inclination - i",
                                                      "argument of perigee - w", "RAAN - W", "TRAN - V"],
                                "Keplerian": [self.a, self.e, self.i, self.w, self.W, self.V], "Coordinate": ['x', 'y', 'z', 'u', 'v', 'w'],
                                "Converted Cartesian": [cart_pos_x, cart_pos_y, cart_pos_z, cart_vel_x, cart_vel_y,
                                                        cart_vel_z],
                                })
        print(Kep2Car)  
        
#Example Input: orbit_1.OrbitPropagation(3600)
#Set desired time interval dt to get orbit parameters (in parentheses)
    
    def OrbitPropagation(self,dt):
        GM = 398600.4415
        #Compute the mean motion
        n = np.sqrt(GM/(self.a**3))
        
        # Compute the eccentric anomaly at t=t0
        cos_Eo = ((r * np.cos(self.V)) / self.a) + self.e
        sin_Eo = (r * np.sin(self.V)) / (self.a * np.sqrt(1 - self.e ** 2))
    
        # adding 2 pi for for very small values ensures that Eo stays in the range 0< Eo <2Pi
        Eo = math.atan2(sin_Eo, cos_Eo)
        if Eo < 0.0:
            Eo = Eo + 2 * np.pi
        else:
            Eo = Eo
    
        # Compute mean anomaly at start point
        Mo = Eo - self.e * np.sin(Eo)  # From Kepler's equation
    
        # Compute the mean anomaly at t+dt
        Mi = Mo + n * dt
    
        # Solve Kepler's equation to compute the eccentric anomaly at t+dt
        M = Mi

        #Minimal value approaching 0 (level of accuracy)
        min_val = 1E-7
        
        #Initial Guess at Eccentric Anomaly (taken these conditions from Fundamentals of Astrodynamics by Roger.E.Bate)
        if M < np.pi:
            E = M + (self.e / 2)
        if M > np.pi:
            E = M - (self.e / 2)
            
        #Initial Conditions
        f = E - self.e*np.sin(E) - M
        f_prime = 1 - self.e*np.cos(E)
        ratio = f / f_prime
        
        #Numerical iteration for ratio compared to level of accuracy wanted
        
        iteration_array = []
        while abs(ratio) > min_val:
            f = E - self.e*np.sin(E) - M
            f_prime = 1 - self.e*np.cos(E)
            ratio = f / f_prime
            
            iteration_array.append(1)
         
            if abs(ratio) > min_val:
                 E = E - ratio
            if abs(ratio) < min_val:
                break

        Ei = E
    
        # Compute the gaussian vector component x,y
        x_new = self.a * (np.cos(Ei) - self.e)
        y_new = self.a * ((np.sqrt(1 - self.e ** 2)) * (np.sin(Ei)))
    
        # Compute the in-orbital plane Gaussian Vectors
        # This gives P and Q in ECI components
    
        P = np.matrix([[np.cos(self.W) * np.cos(self.w) - np.sin(self.W) * np.cos(self.i) * np.sin(self.w)],
                        [np.sin(self.W) * np.cos(self.w) + np.cos(self.W) * np.cos(self.i) * np.sin(self.w)],
                        [np.sin(self.i) * np.sin(self.w)]])
    
        Q = np.matrix([[-np.cos(self.W) * np.sin(self.w) - np.sin(self.W) * np.cos(self.i) * np.cos(self.w)],
                        [-np.sin(self.W) * np.sin(self.w) + np.cos(self.W) * np.cos(self.i) * np.cos(self.w)],
                        [np.sin(self.i) * np.cos(self.w)]])
    
        # Compute the position vector at t+dt
    
        # We know the inertial vector components along the P and Q vectors.
        # Thus we can project the satellite position onto the ECI basis.
    
        cart_pos_x_new = (x_new * P.item(0)) + (y_new * Q.item(0))
        cart_pos_y_new = (x_new * P.item(1)) + (y_new * Q.item(1))
        cart_pos_z_new = (x_new * P.item(2)) + (y_new * Q.item(2))
    
        print('x-coordinate:', cart_pos_x_new)
        print('y-coordinate:', cart_pos_y_new)
        print('z-coordinate:', cart_pos_z_new)
    
        # Compute the range at t+dt
        r_new = self.a * (1 - self.e * (np.cos(Ei)))
    
        # Compute the gaussian velocity components
        cos_Ei = ((x_new / self.a) + self.e)
        sin_Ei = (y_new / self.a * np.sqrt(1 - self.e ** 2))
    
        f_new = (np.sqrt(self.a * GM)) / r_new
        g_new = np.sqrt(1 - self.e ** 2)
    
        cart_vel_x_new = (-f_new * sin_Ei * P.item(0)) + (
                f_new * g_new * cos_Ei * Q.item(0))  # x component of velocity - a.k.a 'u'
        cart_vel_y_new = (-f_new * sin_Ei * P.item(1)) + (
                f_new * g_new * cos_Ei * Q.item(1))  # y component of velocity - a.k.a 'v'
        cart_vel_z_new = (-f_new * sin_Ei * P.item(2)) + (
                f_new * g_new * cos_Ei * Q.item(2))  # z component of velocity - a.k.a 'w'
    
        print('u-coordinate:', cart_vel_x_new)
        print('v-coordinate:', cart_vel_y_new)
        print('w-coordinate:', cart_vel_z_new)
    
        Ephermeris = pd.DataFrame(
            {"Orbital Parameter": ["semi-major axis", "eccentricity", "inclination", "argument of perigee", "RAAN", "TRAN"],
              "Keplerian": [self.a, self.e, self.i, self.w, self.W, self.V],
              "Initial Cartesian Coordinates(t=t0)": [cart_pos_x, cart_pos_y, cart_pos_z, cart_vel_x, cart_vel_y,
                                                      cart_vel_z],
              "New Cartesian Coordinates(t+dt)": [cart_pos_x_new, cart_pos_y_new, cart_pos_z_new, cart_vel_x_new,
                                                  cart_vel_y_new, cart_vel_z_new],
              })
        Ephermeris
    
    
