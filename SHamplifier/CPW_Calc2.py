# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 14:12:56 2019

@author: iwh
"""

### See R.N. Simons, "Coplanar Waveguide Circuits, Components, and Systems" Wilwy and Sons (2001); as referecne and for derivations

import phidl as pd
import numpy as np
from scipy.special import ellipk
from scipy.optimize import minimize


thickness1 = 350 #substrate thickness (sapphire)
oxide = 0.3

#eps1 = 12
eps2 = 4
light = 2.997e8


### Finite Thickness substrate with single dielectric constant
def epseff(g, w, h, eps):
    return 1 + (eps-1)*0.5*(ellipk(np.sinh(np.pi*w*0.25/h)/np.sinh(np.pi*(w + 2*g)*0.25/h))/ellipk(np.sqrt(1-(np.sinh(np.pi*w*0.25/h)/np.sinh(np.pi*(w+2*g)*0.25/h))**2)))*(ellipk(np.sqrt(1-(w/(w+2*g))**2))/ellipk(w/(w+2*g)))

def epseffft(g, w, h, eps,t=0.1):
    delta = (1.25*t/np.pi)*(1 + np.log(4*np.pi*w/t)) 
    w = w+ delta
    g = g-delta    
    k1 = w/(w+2*g)
    k2 = np.sinh(np.pi*w*0.25/h)/np.sinh(np.pi*(w + 2*g)*0.25/h)
    epseff= 1 + (eps-1)*0.5*(ellipk(k2)/ellipk(np.sqrt(1-k2**2)))*(ellipk(np.sqrt(1-k1**2))/ellipk(k1))
    epsefft = epseff- (0.7*(1-epseff)*t/g)/(ellipk(k1)/ellipk(np.sqrt(1-k1**2))+0.7*t/g)
    return epsefft

### Finite width ground planes, with Finite Thickness substrate with single dielectric constant. Eq 4.20 (!) in Simons
def epseff_g(g, w, h, eps,gndw=1):
    # hard code gnd as function of w until I get better at param lists
    gnd = w*gndw #assume gnds = w *gndw
    a = w/2
    b = (2*g + 2*a)/2
    c = (2*gnd + 2*b)/2
    k = (c/b)*np.sqrt((b**2 - a**2)/(c**2 - a**2))
    kp = np.sqrt(1-k**2)
    
    k5 =(np.sinh(np.pi*c*0.5/h)/np.sinh(np.pi*b*0.5/h))*np.sqrt((np.sinh(np.pi*b*0.5/h)**2 - np.sinh(np.pi*a*0.5/h)**2)/(np.sinh(np.pi*c*0.5/h)**2 - np.sinh(np.pi*a*0.5/h)**2)) ### ian has 0.25 where 0.5 is, as in Eq 2.3?...check limit
    kp5 = np.sqrt(1-k5**2)
    return 1 + (-1+eps)*0.5*ellipk(k)*ellipk(kp5)/(ellipk(kp)*ellipk(k5))

def epseff_ft(g, w, h, eps,t=0.1, gndw=1):
    # hard code gnd as function of w until I get better at param lists
    # calculate epseff based on modified W and gap due to thickness of conductor. This only really works for conductors thicker than the skin depth. 
    # from: http://qucs.sourceforge.net/tech/node86.html#SECTION001316000000000000000
    delta = (1.25*t/np.pi)*(1 + np.log(4*np.pi*w/t)) 
    print(f'delta = {delta:.3f}')
    weff = w + delta
    geff = g - delta
    k1 = weff/(weff+2*geff)
    gnd = weff*gndw #assume gnds = w *gndw
    a = weff/2
    b = (2*geff + 2*a)/2
    c = (2*gnd + 2*b)/2
    k = (c/b)*np.sqrt((b**2 - a**2)/(c**2 - a**2))
    kp = np.sqrt(1-k**2)
    
    k5 =(np.sinh(np.pi*c*0.5/h)/np.sinh(np.pi*b*0.5/h))*np.sqrt((np.sinh(np.pi*b*0.5/h)**2 - np.sinh(np.pi*a*0.5/h)**2)/(np.sinh(np.pi*c*0.5/h)**2 - np.sinh(np.pi*a*0.5/h)**2)) ### ian has 0.25 where 0.5 is, as in Eq 2.3?...check limit
    kp5 = np.sqrt(1-k5**2)
    return 1 + (-1+eps)*0.5*ellipk(k)*ellipk(kp5)/(ellipk(kp)*ellipk(k5))


### Two layer substrate with individual dielectric constants
def epseff_d(eps1, eps2, h1, h2, w, g):
    k0 = w/(w + 2.0*g)
    k0_p = np.sqrt(1-np.square(k0))
    k1 = np.sinh(np.pi*w/(h1 + h2)*0.25)/np.sinh((np.pi*(w + 2.0*g))/(4.0*(h1+h2)))
    k1_p = np.sqrt(1 - np.square(k1))
    k2 = np.sinh(np.pi*w/(h2)*0.25)/np.sinh((np.pi*(w + 2.0*g))/(4.0*(h2)))
    k2_p = np.sqrt(1 - np.square(k2))
    e_eff = 1 + ((eps1 - 1)/2.0)*((ellipk(k1)/ellipk(k1_p))*(ellipk(k0_p)/ellipk(k0))) + ((eps2 - eps2)/2.0)*((ellipk(k2)/ellipk(k2_p))*(ellipk(k0_p)/ellipk(k0)))
    return e_eff

### Conductive Backplane CPW
def epseff_CBCPW(eps, h, w, g):
    k = w/(w + 2.0*g)
    k_p = np.sqrt(1.0 - np.square(k))
    k3 = np.tanh(np.pi*w/(4.0*h))/np.tanh(np.pi*(w + 2*g)/(4.0*h))
    k3_p = np.sqrt(1.0 - np.square(k3))
    e_eff = (1 + eps * (ellipk(k_p)/ellipk(k))*(ellipk(k3)/ellipk(k3_p)))/(1 + (ellipk(k_p)/ellipk(k))*(ellipk(k3)/ellipk(k3_p)))
    return e_eff

def Z0(g, w, h):
    eps = epseff(g, w, h)
    return (np.pi*29.979*ellipk(np.sqrt(1-(w/(w+2*g))**2)))/(np.sqrt(eps)*ellipk(w/(w+2*g)))

def Z0_V1(g, eps, h, w):
    k0 = w/(w+2.0*g)
    k0_p = np.sqrt(1-np.square(k0))    
    return (29.979*np.pi/np.sqrt(epseff(eps, h, w, g)))*(ellipk(k0_p)/ellipk(k0))

def Z0_V2(g, w, h, eps_eff, eps):
    k0 = w/(w+2.0*g)
    k0_p = np.sqrt(1-np.square(k0))
    return (29.979*np.pi/np.sqrt(eps_eff(g, w, h, eps)))*(ellipk(k0_p)/ellipk(k0))

def Z0_V2test(g, w, h, eps_eff, eps):
    k0 = w/(w+2.0*g)
    k0_p = np.sqrt(1-np.square(k0))
    print(f'eps_eff = {eps_eff(g, w, h, eps):.4f}, k0 = {k0:.4f}, ratio = {(ellipk(k0_p)/ellipk(k0)):.4f}')    
    return (29.979*np.pi/np.sqrt(eps_eff(g, w, h, eps)))*(ellipk(k0_p)/ellipk(k0))



def Z0_CBCPW(eps_eff, h, w, g):
    k = w/(w + 2.0*g)
    k_p = np.sqrt(1.0 - np.square(k))
    k3 = np.tanh(np.pi*w/(4.0*h))/np.tanh(np.pi*(w + 2*g)/(4.0*h))
    k3_p = np.sqrt(1.0 - np.square(k3))
    zn = (np.pi*2*29.979/np.sqrt(eps_eff))*1/((ellipk(k)/ellipk(k_p)) + (ellipk(k3)/ellipk(k3_p)))
    return zn

#return necessary gap for 50ohm line given center conductor width
def getgap(w, h, eps_eff, impedance,  eps):
    gval = w/2
    while (np.abs(Z0_V2(gval, w, h, eps_eff, eps) - impedance) > 0.05):
        if Z0_V2(gval, w, h, eps_eff, eps)> impedance:
            gval = gval*0.98
        else:
            gval = gval*1.02
#    print("gap= "+ str(gval) + "; eps = " + str(epseff(gval, w,h,eps)))
    return gval


#### test functions for finite width gnds
def getgaptest(w, h, impedance, eps, gndw):
    gval = w/2
    while (np.abs(Z0_V2_test(gval, w, h, eps,gndw) - impedance) > 0.05):
        if Z0_V2_test(gval, w, h, eps,gndw)> impedance:
            gval = gval*0.98
        else:
            gval = gval*1.02
    print("gap= "+ str(gval) + "; eps = " + str(epseff_g(gval, w,h,eps,gndw=gndw)))
    return gval

def Z0_V2_test(g, w, h, eps,gndw):
    k0 = w/(w+2.0*g)
    k0_p = np.sqrt(1-np.square(k0))    
    return (29.979*np.pi/np.sqrt(epseff_g(g, w, h, eps,gndw=gndw)))*(ellipk(k0_p)/ellipk(k0))


def getpads(spacing, thickness, eps_eff, impedance, eps):
    h = thickness
    w = spacing/2
    gap = spacing - w
    while (np.abs(Z0_V2(gap, w, h, eps_eff, eps) - impedance) > 0.05):
        if Z0_V2(gap, w, h, eps_eff, eps)> impedance:
            gap = gap*0.98
            w = spacing - gap
        else:
            gap = gap*1.02
            w = spacing - gap
            
    return round(w,3), round(gap,3)

# def getgap_2 (w, h = thickness1, eps = eps1, eps_eff = epseff):
#     gap = w/2
#     #eps_eff = epseff(gap, h, w)
#     res = minimize(Z0_V2, gap, args = (w, h, eps_eff, eps))
#     return res.x

def half_wavelength(g, w, h, eps, freq):
    epsval = epseff(g, w, h, eps)
    lam_2 = 0.5*light/(np.sqrt(epsval)*freq)
    return lam_2     

def quarter_wavelength(g, w, h, eps, freq):
    epsval = epseff(g, w, h, eps)
    lam_4 = 0.25*light/(np.sqrt(epsval)*freq)
    return lam_4

def getfreq(g, w, h, eps, lr):
    epsval = epseff(g, w, h, eps)
    freq = 0.5*light/(np.sqrt(epsval)*lr)
    return freq
    
#print(getfreq(23.8, 53, 580, 10.8, 10.8e-3)/1e9)





    
    
    
    
    
    
    