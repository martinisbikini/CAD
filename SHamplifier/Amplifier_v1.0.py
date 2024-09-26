# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 16:16:09 2021

@author: iwh
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 14:28:15 2019

@author: iwh
"""

import phidl as pd
import numpy as np
from scipy.special import ellipk
import CPW_Calc2
import phidl.utilities
import AmpDevices as ad
from sys import exit

# NOTE: some of these functions used these global parameters (every gap calculator will)
#sub_eps = 10.8      #Real part of substrate Dielectric Constant-Sapphire
sub_eps = 11.9      #Real part of substrate Dielectric Constant--Si
thickness = 380     #Substrate Thickness in micron
impedance = 50      #CPW Design impedance
eps_eff = CPW_Calc2.epseff_g      #Effective epsilon calcualtion to use for calculation. CPW_Calc has several functions for multilayer substrates, conductive back plane, etc.


# ybco_layer = pd.Layer(gds_layer=0, name='ybcolayer',color='red')
# gold_layer = pd.Layer(gds_layer=1, name='goldlayer',color='yellow')  
# mag_layer = pd.Layer(gds_layer = 2, name = 'filmlayer', color = 'black')
# via_layer = pd.Layer(gds_layer = 3, name = 'filmlayer', color = 'purple', alpha = .06, dither = 'dotted')
# fast_layer = pd.Layer(gds_layer = 4, name = 'fastlayer', color = 'green')

#ybco_layer = pd.Layer(gds_layer = 2, name='ybcolayer', color='red')     #Called YBCO_layer for historical reasons. This layer defines the CPW gaps
gold_layer = pd.Layer(gds_layer = 1, name='goldlayer', color='yellow')  #Contact pad layer
mag_layer = pd.Layer(gds_layer = 2, name = 'filmlayer', color = 'purple')
via_layer = pd.Layer(gds_layer = 3, name = 'etchlayer', color = 'green')
wafer_layer = pd.Layer(gds_layer = 4, name = 'waferlayer', color = 'red')
mark_layer = pd.Layer(gds_layer = 5, name = 'marklayer', color = 'blue')
#via_layer = pd.Layer(gds_layer = 3, name = 'filmlayer', color = 'green', alpha = .06, dither = 'dotted')



#%%

output_directory = r'Y:\Martina Kiechle\CAD\SHamplifier'
#######################   FULL DEVICE STRUCTURE ##############################

#CPW paprameters
numperiods = 0
width = 0.8 # center conductor width in narrow portion
padlength = 1000
taper_length = 250  # Length of exponential taper
GSGspacing = 150       #GSG contact pad spacing
swidegndw = 100

#Spinwave channel parameters
overlaptrim = 0 #trim the spinwave channel so the CPWs extend past it (in micron)
straightlength = .05 # CPW length or Spinwave channel width
separation = width*8 + 1 # separation between the center conductors in the 2 CPWs, i.e., the length of the active spinwave channel
extend = 10  #disteance to extend the spinwave channel past the CPW
extendpastswchannel = 1
dcleadwidth = 10
dctaperlength = 20

# DC contact parameters and fiducials

ProbeSpacing = 600 # spacing between the two GSG probes
via_extend = 0.4
numCh = 5
ydistCh = 1



#####Single Deivice
ampdevice = pd.Device('myAmp')
temp = ad.AmplifierDevice_v4(width, straightlength,via_extend,1,10,20, numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, separation, extend, ProbeSpacing)
ampdevice.add_ref(temp)
pd.quickplot2(ampdevice)




# ampdevice = pd.Device('myAmp')
# cpw = ad.CPW_pair(width, straightlength, extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation)
# DCpath = ad.DCpath(width, straightlength, via_extend, dcleadwidth, dctaperlength, numperiods, padlength, taper_length, overlaptrim, separation, extend, ProbeSpacing)
# ampdevice.add_ref(DCpath)
# pd.quickplot2(ampdevice)

exit()



#%% 2 um long sw channels with 200 nm waveguides and varying spinwave channel widths (x4)

######Construct Complete Die layout
stride = 1250 # stride to place devs on chip
vstride = 3000 # vertical stride
widths = [0.1,0.15,0.2,0.4,.6,.8,1,2]
#swchannel = 2 # length of spinwave channel
wg_width = [.2,0.4,0.6,0.8] # waveguide width
extendpast = 1 # how far to go past the swchannel
dclw = 5 # dc lead width
numperiods = 0
# straight_lens = [.02,.03,.04, .05,.06, .075, .1,.15, .2,0.3, .4, .6, .8, 1, 2, 3,4,5] # widths of the spin wave channels

straight_lens = [.1,.2,.3,.4, .5,.6, .75,.8,.8, 1, 1, 2, 2, 3, 3, 4, 4, 5] # widths of the spin wave channels

straight_lens = np.sort(straight_lens)

# make a circle to define the wafer
from phidl import quickplot as qp
wafer = pd.Device('wafer')
die_size = 12000
ref_size = 2000
ndierows = 1
ndiecols = 3
waferedge = pd.geometry.circle(radius = 76200/2, layer = wafer_layer)
centerRef = pd.geometry.rectangle((ref_size, ref_size), layer = via_layer).move((-ref_size/2, -ref_size/2))

diearea = pd.geometry.rectangle((ndiecols*die_size,ndiecols*die_size),layer = wafer_layer).move((-ndiecols*die_size/2,-ndiecols*die_size/2))
boolarea = pd.geometry.boolean(waferedge,diearea, 'not', layer = wafer_layer)
wafer.add_ref(boolarea)
#wafer.write_gds(output_directory + 'wafertest0', unit=1e-6, precision=1e-11)



for c in range(len(wg_width)):
    AmpDeviceArray = pd.Device('ampdevicearray')
    
    swchannel = wg_width[c] * 8 
    extend = wg_width[c] * 3
    count = 0
    vcount = 0
    for w in range(6):
        ######for fixed spinwave channel dimensions######
        temp = ad.AmplifierDevice_v4(wg_width[c], straight_lens[w],via_ext,extendpast,dclw, 2*(dclw-straight_lens[w]), numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, swchannel, extend, ProbeSpacing, YIG, YIG_contact_width)
        #temp = AmplifierDevice(wg_width, straightlength, numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, separation, extend, ProbeSpacing, YIG, YIG_contact_width)
        
        ######for dynamic spinwave channel dimensions######
        #temp = AmplifierDevice(width, straightlength, numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, (100*width), extend, ProbeSpacing, YIG, YIG_contact_width)
    
    
        AmpDeviceArray.add_ref(temp).move([stride*count, -vstride*vcount])
        count +=1

    
    
    count = 0
    vcount = 1
    for w in range(6,12):
        temp = ad.AmplifierDevice_v4(wg_width[c], straight_lens[w],via_ext,extendpast,dclw, 2*(dclw-straight_lens[w]), numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, swchannel, extend, ProbeSpacing, YIG, YIG_contact_width)
        AmpDeviceArray.add_ref(temp).move([stride*count, -vstride*vcount])
        count +=1
    
    numperiods = 0
    count = 0
    vcount = 2
    for w in range(12,18):
        temp = ad.AmplifierDevice_v4(wg_width[c], straight_lens[w],via_ext,extendpast,dclw, 2*(dclw-straight_lens[w]), numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, swchannel, extend, ProbeSpacing, YIG, YIG_contact_width)
        AmpDeviceArray.add_ref(temp).move([stride*count, -vstride*vcount])
        count +=1
    
    
    centerx=((AmpDeviceArray.xmin + AmpDeviceArray.xmax)/2) 
    centery= ((AmpDeviceArray.ymin + AmpDeviceArray.ymax)/2)

    AmpDeviceArray.move([-centerx, -centery])

    #verify
    centerx=((AmpDeviceArray.xmin + AmpDeviceArray.xmax)/2) 
    centery= ((AmpDeviceArray.ymin + AmpDeviceArray.ymax)/2)
    print(centerx)
    print(centery)

    #Print GDSII file and plot
    #AmpDeviceArray.write_gds(output_directory + '\\Amplifier-Device-Full-CHIP-'+str(wg_width[c]*1000)+'nm', unit=1e-6, precision=1e-11)
    
    
    if c == 0:
        wafer.add_ref(AmpDeviceArray).move((die_size, die_size))
        wafer.add_ref(AmpDeviceArray).move((0, die_size))
        
    if c == 1:
        wafer.add_ref(AmpDeviceArray).move((-die_size, die_size))
        wafer.add_ref(AmpDeviceArray).move((die_size, 0))
        
    if c == 2:
        wafer.add_ref(AmpDeviceArray).move((-die_size, 0))
        wafer.add_ref(AmpDeviceArray).move((die_size, -die_size))
        
        
    if c == 3:
        wafer.add_ref(AmpDeviceArray).move((0, -die_size))
        wafer.add_ref(AmpDeviceArray).move((-die_size, -die_size))
    
    pd.quickplot2(AmpDeviceArray)


#%%


grid = pd.Device('grid')
cutline = pd.geometry.rectangle((3*die_size,50), layer = mark_layer)

grid.add_ref(cutline).move((-cutline.xmax/2,diearea.ymax-25))
grid.add_ref(cutline).move((-cutline.xmax/2,diearea.ymax-25-12000))
grid.add_ref(cutline).move((-cutline.xmax/2,diearea.ymax-25-24000))
grid.add_ref(cutline).move((-cutline.xmax/2,diearea.ymax-25-36000))

grid.add_ref(cutline).rotate(90).move((-diearea.xmax+25,-cutline.xmax/2))
grid.add_ref(cutline).rotate(90).move((-diearea.xmax+25+12000,-cutline.xmax/2))
grid.add_ref(cutline).rotate(90).move((-diearea.xmax+25+24000,-cutline.xmax/2))
grid.add_ref(cutline).rotate(90).move((-diearea.xmax+25+36000,-cutline.xmax/2))



wafer.add(centerRef)
wafer.add_ref(grid)
qp(wafer)
wafer.write_gds(output_directory + '\\Wisser_AmplifierWafer_v1.0_cutmarks', unit=1e-6, precision=1e-10)