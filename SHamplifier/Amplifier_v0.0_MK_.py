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
# import AmpDevices as ad
# %matplotlib qt
import CADfunctions_MK as ad

from sys import exit

# NOTE: some of these functions used these global parameters (every gap calculator will)
#sub_eps = 10.8      #Real part of substrate Dielectric Constant-Sapphire
sub_eps = 11.9      #Real part of substrate Dielectric Constant--Si
thickness = 380     #Substrate Thickness in micron
impedance = 50      #CPW Design impedance
eps_eff = CPW_Calc2.epseff_g      #Effective epsilon calcualtion to use for calculation. CPW_Calc has several functions for multilayer substrates, conductive back plane, etc.
spacing = 150


gold_layer = pd.Layer(gds_layer = 1, name='goldlayer', color='yellow')  #Contact pad layer
mag_layer = pd.Layer(gds_layer = 2, name = 'filmlayer', color = 'purple')
via_layer = pd.Layer(gds_layer = 3, name = 'etchlayer', color = 'green')
wafer_layer = pd.Layer(gds_layer = 4, name = 'waferlayer', color = 'red')
mark_layer = pd.Layer(gds_layer = 5, name = 'marklayer', color = 'blue')


#%%

output_directory = r'Y:\Martina Kiechle\CAD\SHamplifier'
#######################   FULL DEVICE STRUCTURE ##############################
#Spinwave channel parameters
overlaptrim = 0 #trim the spinwave channel so the CPWs extend past it (in micron)
straightlength = 1 # CPW length or Spinwave channel width
separation = 5 # separation between the center conductors in the 2 CPWs, i.e., the length of the active spinwave channel
extend = 5  #disteance to extend the spinwave channel past the CPW
extendpastswchannel = 0

# DC contact parameters and fiducials
ProbeSpacing = 600 # spacing between the two GSG probes
via_extend_x = 0.5
via_extend_y = 0.4
dc_padlength = 300
dcleadwidth = 10
dctaperlength = 150
constlen = 2
constrwdth = 0.2

numCh = 7
ydistCh = straightlength*4


#CPW paprameters
numperiods = 0
width = 0.5 # center conductor width in narrow portion
padlength = 1000
taper_length = 250  # Length of exponential taper
GSGspacing = 150       #GSG contact pad spacing
swidegndw = 100


ampdevice = pd.Device('myDevice')

if numCh == 1:
    cpw = ad.CPW_pair(width, numCh*straightlength+(numCh)*ydistCh, extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation)
else:
    cpw = ad.CPW_pair(width, numCh*straightlength+(numCh-1)*ydistCh, extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation)


ampdevice.add_ref(cpw)
centerx = ((ampdevice.xmin + ampdevice.xmax)/2) 
centery = ((ampdevice.ymin + ampdevice.ymax)/2)



DCpath = ad.DCpath(width, straightlength, via_extend_x, via_extend_y, dcleadwidth, dctaperlength, numCh, ydistCh, dc_padlength, taper_length, overlaptrim, separation, extend, ProbeSpacing, constlen, constrwdth)

ampdevice.add_ref(DCpath).move([centerx, centery])

pd.quickplot2(ampdevice)

ampdevice.write_gds(output_directory + '\\AmpDev_array_Testdevice', unit=1e-6, precision=1e-10)


from phidl import quickplot as qp

straightlength = 0.5
numCh = 35
numRow = 6

DCtest_device = pd.Device('myDevice')

pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
ydistCh = pad_w+pad_gap

DCtest = ad.DCtest_path(straightlength, via_extend_x, via_extend_y, dcleadwidth, dctaperlength, numCh, numRow, dc_padlength, overlaptrim, separation, extend, ProbeSpacing, constlen, constrwdth)
DCtest_device.add_ref(DCtest).move([centerx, centery])

# pd.quickplot2(DCtest_device)
# DCtest_device.write_gds(output_directory + '\\DC_Testdevice', unit=1e-6, precision=1e-10)


exit()



#%% 2 um long sw channels with 200 nm waveguides and varying spinwave channel widths (x4)

######Construct Complete Die layout
stride = 1250 # stride to place devs on chip
vstride = 3000 # vertical stride
widths = [0.1,0.15,0.2,0.4,.6,.8,1,2]
wg_width = [0.2,0.4,0.6,0.8] # waveguide width
extendpast = 1 # how far to go past the swchannel
dcleadwidth = 10


straight_lens = [.1,.2,.3,.4, .5,.6, .75,.8,.8, 1, 1, 2, 2, 3, 3, 4, 4, 5] # widths of the spin wave channels
straight_lens = np.sort(straight_lens)


#%% Make cells and place on wafer old devices

die_size = 10000

for c in range(len(wg_width)):
    AmpDeviceArray = pd.Device('ampdevicearray')
    
    extend = 10
    swchannel = 50
    count = 0
    vcount = 0
    
    for w in range(6):
        ######for fixed spinwave channel dimensions######
        temp = ad.kdepdamping_v0(straight_lens[w],wg_width[c],via_ext,extendpast,dclw, 2*(dclw-straight_lens[w]), numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, swchannel, extend, ProbeSpacing, YIG, YIG_contact_width)
        
        
        ######for dynamic spinwave channel dimensions######
        #temp = AmplifierDevice(width, straightlength, numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, (100*width), extend, ProbeSpacing, YIG, YIG_contact_width)
    
    
        AmpDeviceArray.add_ref(temp).move([stride*count, -vstride*vcount])
        count +=1
 
    
    count = 0
    vcount = 1
    for w in range(6,12):

        temp = ad.kdepdamping_v0(straight_lens[w],wg_width[c],via_ext,extendpast,dclw, 2*(dclw-straight_lens[w]), numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, swchannel, extend, ProbeSpacing, YIG, YIG_contact_width)
        
        AmpDeviceArray.add_ref(temp).move([stride*count, -vstride*vcount])
        count +=1
    
    numperiods = 0
    count = 0
    vcount = 2
    for w in range(12,18):

        temp = ad.kdepdamping_v0(straight_lens[w],wg_width[c],via_ext,extendpast,dclw, 2*(dclw-straight_lens[w]), numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, swchannel, extend, ProbeSpacing, YIG, YIG_contact_width)
        
        AmpDeviceArray.add_ref(temp).move([stride*count, -vstride*vcount])
        count +=1
    
    
    centerx=((AmpDeviceArray.xmin + AmpDeviceArray.xmax)/2) 
    centery= ((AmpDeviceArray.ymin + AmpDeviceArray.ymax)/2)

    AmpDeviceArray.move([-centerx, -centery])
    
    
    #cutmarks = hd.cutmarks(50,500,die_size)
    #AmpDeviceArray.add_ref(cutmarks)

    #verify
    centerx=((AmpDeviceArray.xmin + AmpDeviceArray.xmax)/2) 
    centery= ((AmpDeviceArray.ymin + AmpDeviceArray.ymax)/2)
    print(centerx)
    print(centery)

    #Print GDSII file and plot
    AmpDeviceArray.write_gds(output_directory + '\\k2_device-CPW-'+str(wg_width[c]*1000)+'nm', unit=1e-6, precision=1e-11)
        
    pd.quickplot2(AmpDeviceArray)
    
    

# Make cells and place on wafer new devices
extendpast = 5
for c in range(len(wg_width)):
    AmpDeviceArray_new = pd.Device('ampdevicearray_new')
    
    #swchannel = wg_width[c] * 8
    #swchannel = wg_width[c] * 5
    extend = 5
    swchannel = 50
    count = 0
    vcount = 0
    for w in range(6):
        ######for fixed spinwave channel dimensions######
        temp = ad.kdepdamping_gs_v01(straight_lens_gs[w],wg_width[c],via_ext,extendpast,dclw, 2*(dclw-straight_lens_gs[w]), numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, swchannel, extend, ProbeSpacing, YIG, YIG_contact_width)
        
        
        ######for dynamic spinwave channel dimensions######
        #temp = AmplifierDevice(width, straightlength, numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, (100*width), extend, ProbeSpacing, YIG, YIG_contact_width)
    
    
        AmpDeviceArray_new.add_ref(temp).move([stride*count, -vstride*vcount])
        count +=1
  
    count = 0
    vcount = 1
    for w in range(6,12):
        temp = ad.kdepdamping_gs_v01(straight_lens_gs[w],wg_width[c],via_ext,extendpast,dclw, 2*(dclw-straight_lens_gs[w]), numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, swchannel, extend, ProbeSpacing, YIG, YIG_contact_width)
        
        AmpDeviceArray_new.add_ref(temp).move([stride*count, -vstride*vcount])
        count +=1
    
    numperiods = 0
    count = 0
    vcount = 2
    for w in range(12,18):

        temp = ad.kdepdamping_gs_v01(straight_lens_gs[w],wg_width[c],via_ext,extendpast,dclw, 2*(dclw-straight_lens_gs[w]), numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, swchannel, extend, ProbeSpacing, YIG, YIG_contact_width)
        
        AmpDeviceArray_new.add_ref(temp).move([stride*count, -vstride*vcount])
        count +=1
    
    
    centerx=((AmpDeviceArray_new.xmin + AmpDeviceArray_new.xmax)/2) 
    centery= ((AmpDeviceArray_new.ymin + AmpDeviceArray_new.ymax)/2)

    AmpDeviceArray_new.move([-centerx, -centery])
    

    centerx=((AmpDeviceArray_new.xmin + AmpDeviceArray_new.xmax)/2) 
    centery= ((AmpDeviceArray_new.ymin + AmpDeviceArray_new.ymax)/2)
    print(centerx)
    print(centery)

    #Print GDSII file and plot
    AmpDeviceArray_new.write_gds(output_directory + '\\k2_device-GS'+str(wg_width[c]*1000)+'nm', unit=1e-6, precision=1e-11)
    
    pd.quickplot2(AmpDeviceArray)




#%% Make quarter

die_size = 10000
QuarterWafer = pd.Device('QuarterWafer')
QuarterWafer.add_ref(AmpDeviceArray)
QuarterWafer.add_ref(AmpDeviceArray).move((die_size,0))
QuarterWafer.add_ref(AmpDeviceArray_new).move((0,-die_size))
QuarterWafer.add_ref(AmpDeviceArray_new).move((die_size,-die_size))

hd.CenterDev(QuarterWafer)
#pd.quickplot(QuarterWafer)


#%% Place on wafer

wafer = pd.Device('wafer')
ref_size = 2000
ndierows = 5
ndiecols = 5

# make a circle to define the wafer

waferedge = pd.geometry.circle(radius = 76200/2, layer = wafer_layer)
centerRef = pd.geometry.rectangle((ref_size, ref_size), layer = {via_layer,mag_layer}).move((-ref_size/2, -ref_size/2))

diearea = pd.geometry.rectangle((ndiecols*die_size,ndiecols*die_size),layer = wafer_layer).move((-ndiecols*die_size/2,-ndiecols*die_size/2))
boolarea = pd.geometry.boolean(waferedge,diearea, 'not', layer = wafer_layer)
wafer.add_ref(boolarea)

waferedge = pd.geometry.circle(radius = 76200/2, layer = wafer_layer)
centerRef = pd.geometry.rectangle((ref_size, ref_size), layer = mag_layer).move((-ref_size/2, -ref_size/2))

diearea = pd.geometry.rectangle((ndiecols*die_size,ndiecols*die_size),layer = wafer_layer).move((-ndiecols*die_size/2,-ndiecols*die_size/2))
boolarea = pd.geometry.boolean(waferedge,diearea, 'not', layer = wafer_layer)
wafer.add_ref(boolarea)


wafer.add_ref(QuarterWafer).move((-1.5*die_size,1.5*die_size))
wafer.add_ref(QuarterWafer).move((1.5*die_size,1.5*die_size))
wafer.add_ref(QuarterWafer).move((-1.5*die_size,-1.5*die_size))
wafer.add_ref(QuarterWafer).move((1.5*die_size,-1.5*die_size))

EBLmark = hd.EBL_marks()
wafer.add_ref(EBLmark)
wafer.add_ref(centerRef)
pd.quickplot(wafer)

wafer.write_gds(output_directory + '\\SHamps_wafer', unit=1e-6, precision=1e-10)






