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

output_directory = r'C:\Users\mnk36\OneDrive - NIST\CAD\SHamplifier\CADfiles'
#######################   FULL DEVICE STRUCTURE ##############################
#Spinwave channel parameters
overlaptrim = 0 # trim the spinwave channel so the CPWs extend past it (in micron)
straightlength = 1 # CPW length or Spinwave channel width
separation = 5 # separation between the center conductors in the 2 CPWs, i.e., the length of the active spinwave channel
extend1 = 5  # disteance to extend the spinwave channel past the CPW
extend2 = 0.8
extendpastswchannel = 0

# DC contact parameters and fiducials
ProbeSpacing = 600 # spacing between the two GSG probes
via_extend_x = 0.4
via_extend_y = 0
dc_padlength = 300
dcleadwidth = 10
dctaperlength = 150
constlen = 2
constrwdth = 0.2
constrtpr = 1.5

numCh = 7
ydistCh = straightlength*4


#CPW paprameters
numperiods = 0
width = 0.5 # center conductor width in narrow portion
gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
cpw_width = 2*gap+3*width 
padlength = 1000
taper_length = 250  # Length of exponential taper
GSGspacing = 150       #GSG contact pad spacing
swidegndw = 100

gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
extralength = ((2*gap)+(3*width)) 


ampdevice = pd.Device('myDevice')

# if numCh == 1:
#     cpw = ad.CPW_pair(width, numCh*straightlength+(numCh)*ydistCh, extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation)
# else:
#     cpw = ad.CPW_pair(width, numCh*straightlength+(numCh-1)*ydistCh, extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation)

if numCh == 1:
    cpw = ad.CPW_single(width, numCh*straightlength+(numCh)*ydistCh, extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation)
else:
    cpw = ad.CPW_single(width, numCh*straightlength+(numCh-1)*ydistCh, extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation)


ampdevice.add_ref(cpw)
centerx = ((ampdevice.xmin + ampdevice.xmax)/2) 
centery = ((ampdevice.ymin + ampdevice.ymax)/2)



DCpath = ad.DCpath_single(width, straightlength, via_extend_x, via_extend_y, dcleadwidth, dctaperlength, numCh, ydistCh, dc_padlength, taper_length, overlaptrim, separation, extend1, extend2, ProbeSpacing, constlen, constrwdth,constrtpr)

ampdevice.add_ref(DCpath).move([centerx+constlen/2+constrtpr+extralength-extend1/2, centery+ampdevice.ymax/2-(numCh*straightlength+(numCh)*ydistCh)/2])
# ampdevice.add_ref(DCpath).move([centerx+2*constlen, centery+ampdevice.ymax/2-(numCh*straightlength+(numCh)*ydistCh)/2])

pd.quickplot2(ampdevice)

ampdevice.write_gds(output_directory + '\\AmpDev_array_Testdevice', unit=1e-6, precision=1e-10)
# exit()



#%%

from phidl import quickplot as qp

straightlength = 0.5
numCol = 5
numRow = 35
constlen = 2
constrwdth = 0.1
separation = 5 # separation between the center conductors in the 2 CPWs, i.e., the length of the active spinwave channel
extend1 = 2.8  
extend2 = 2*via_extend_x #0.8

DCtest_device = pd.Device('myDevice')

pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
ydistCh = pad_w+pad_gap
padlength_dc = 300



for j in range(numCol):
    
    
    for i in range(5):
        straight_lens = [1, 2, 3, 4, 5]
        constrwdth = straight_lens
        DCtest = ad.DCtest_path(straight_lens[j], via_extend_x, via_extend_y, dcleadwidth, dctaperlength, extend1, extend2, padlength_dc, overlaptrim, separation, ProbeSpacing, constlen, constrwdth[j],constrtpr)
        DCtest_device.add_ref(DCtest).move([centerx+j*2.8*ProbeSpacing, centery+i*250])

    
    for i in range(5,10):
        extend1 = 5.8

        DCtest = ad.DCtest_path(straight_lens[j], via_extend_x, via_extend_y, dcleadwidth, dctaperlength, extend1, extend2, padlength_dc, overlaptrim, separation, ProbeSpacing, constlen, constrwdth[j],constrtpr)
        DCtest_device.add_ref(DCtest).move([centerx+j*2.8*ProbeSpacing, centery+i*250])
        
    for i in range(10,15):
        extend2 = 5.8

        DCtest = ad.DCtest_path(straight_lens[j], via_extend_x, via_extend_y, dcleadwidth, dctaperlength, extend1, extend2, padlength_dc, overlaptrim, separation, ProbeSpacing, constlen, constrwdth[j],constrtpr)
        DCtest_device.add_ref(DCtest).move([centerx+j*2.8*ProbeSpacing, centery+i*250])
        
    for i in range(15,20):
        extend2 = 0.8
        straight_lens = [0.2, 0.4, 0.5, 0.6, 0.8]
        constrwdth = straight_lens
        DCtest = ad.DCtest_path(straight_lens[j], via_extend_x, via_extend_y, dcleadwidth, dctaperlength, extend1, extend2, padlength_dc, overlaptrim, separation, ProbeSpacing, constlen, constrwdth[j],constrtpr)
        DCtest_device.add_ref(DCtest).move([centerx+j*2.8*ProbeSpacing, centery+i*250])

    for i in range(20,25):
        extend1 = 2.8
        straight_lens = [0.2, 0.4, 0.5, 0.6, 0.8]
        constrwdth = straight_lens
        DCtest = ad.DCtest_path(straight_lens[j], via_extend_x, via_extend_y, dcleadwidth, dctaperlength, extend1, extend2, padlength_dc, overlaptrim, separation, ProbeSpacing, constlen, constrwdth[j],constrtpr)
        DCtest_device.add_ref(DCtest).move([centerx+j*2.8*ProbeSpacing, centery+i*250])

    for i in range(25,30):
        straight_lens = [1, 2, 3, 4, 5]
        constrwdth = np.array(straight_lens)/10
    
        DCtest = ad.DCtest_path(straight_lens[j], via_extend_x, via_extend_y, dcleadwidth, dctaperlength, extend1, extend2, padlength_dc, overlaptrim, separation, ProbeSpacing, constlen, constrwdth[j],constrtpr)
        DCtest_device.add_ref(DCtest).move([centerx+j*2.8*ProbeSpacing, centery+i*250])
     
        
    for i in range(30,35):
        straight_lens = [5, 5, 5, 5, 5]
        # constrwdth = np.array(straight_lens)/10
        extend1 = 2.8+i-30
        DCtest = ad.DCtest_path(straight_lens[j], via_extend_x, via_extend_y, dcleadwidth, dctaperlength, extend1, extend2, padlength_dc, overlaptrim, separation, ProbeSpacing, constlen, constrwdth[j],constrtpr)
        DCtest_device.add_ref(DCtest).move([centerx+j*2.8*ProbeSpacing, centery+i*250])
        
       
# DCtest = ad.DCtest_path(straightlength, via_extend_x, via_extend_y, dcleadwidth, dctaperlength, extend1, extend2, padlength, overlaptrim, separation, ProbeSpacing, constlen, constrwdth,constrtpr)
# DCtest_device.add_ref(DCtest).move([centerx, centery])

# DCtest_device.move([((DCtest_device.xmin + DCtest_device.xmax)/2) ,((DCtest_device.ymin + DCtest_device.ymax)/2)])

# pd.quickplot2(DCtest_device)
# DCtest_device.write_gds(output_directory + '\\DC_Testdevice', unit=1e-6, precision=1e-10)


# exit()

#%% 2 um long sw channels with 200 nm waveguides and varying spinwave channel widths (x4)

######Construct Complete Die layout
stride = 1250 # stride to place devs on chip
vstride = 3000 # vertical stride
widths = [0.5,1,1.5,2,3,4]
wg_width = [0.5] # waveguide width
# extendpast = 1 # how far to go past the swchannel
dcleadwidth = 10
separation = 5
# ydistCh = straightlength*4



straight_lens = np.array([[5,4,3,2,1,.8,  5,5,5,5,5,5,      5,4,3,2,1,.8],
                          [5,4,3,2,1,.8,  5,5,5,5,5,5,      5,4,3,2,1,.8],
                          [5,4,3,2,1,.8,  5,5,5,5,5,5,      5,4,3,2,1,.8]]) # widths of the spin wave channels

# constrictions_w = straight_lens # use this if no constriction wanted 

constrictions_w = np.array([[5,4,3,2,1,.8, .8,.6,.5,.4,.3,.2, 5,4,3,2,1,.8],
                            [5,4,3,2,1,.8, .8,.6,.5,.4,.3,.2, 5,4,3,2,1,.8],
                            [5,4,3,2,1,.8, .8,.6,.5,.4,.3,.2, 5,4,3,2,1,.8]]) #contriction of the spin wave channels



ch_spc = 2
ydistCh = np.array(straight_lens)+ch_spc

# straight_lens = np.sort(straight_lens)
constlen = 2
constrwdth = 0.2
constrtpr = 1.5
extend1 = 5.8  
extend2 = 2.8


add_cdt_ext = 1 # adding a conduit extention beyond the via


# exit()

#%% Make cells for CPW devices

die_size = 10000

for c in range(len(straight_lens)):
    AmpDeviceArray = pd.Device('ampdevicearray')
    
    extend = 5
    # swchannel = 50
    count = 0
    vcount = 0
    
    for w in range(6):
   
        numCh = 7
        constrwdth = constrictions_w[c,w]
        # cpw = ad.CPW_pair(wg_width[c], numCh*straightlength+(numCh-1)*ydistCh, extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation)
        cpw = ad.CPW_single(wg_width[0], numCh*ydistCh[c,w], extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation)
        
        
        # dcpath = ad.DCpath(width, straightlength, via_extend_x, via_extend_y, dcleadwidth, dctaperlength, numCh, ydistCh, dc_padlength, taper_length, overlaptrim, separation, extend, ProbeSpacing, constlen, constrwdth)
        dcpath = ad.DCpath_single(widths[w], straight_lens[c,w], via_extend_x, via_extend_y, dcleadwidth, dctaperlength, numCh, ydistCh[c,w], dc_padlength, taper_length, add_cdt_ext, separation, extend1, extend2, ProbeSpacing, constlen, constrwdth,constrtpr)
        
        centerx = ((dcpath.xmin + dcpath.xmax)/2) 
        centery = ((dcpath.ymin + dcpath.ymax)/2)
        
        # dcpath.move([centerx+constrtpr, cpw.ymax-(numCh+1)*ydistCh[w]/2-straight_lens[w]/4])
        # dcpath.move([centerx, cpw.ymax-(numCh*straight_lens[c,w]+(numCh-1)*ydistCh[c,w])/2])
        dcpath.move([centerx, centery])
        cpw.move([centerx, centery-padlength-taper_length-numCh*ydistCh[c,w]/2-ch_spc/2])
        
        
        ######for dynamic spinwave channel dimensions######
        #temp = AmplifierDevice(width, straightlength, numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, (100*width), extend, ProbeSpacing, YIG, YIG_contact_width)
        
        ad.fid_squares(dcpath,50,10)
    
        AmpDeviceArray.add_ref(cpw).move([stride*count, -vstride*vcount])
        AmpDeviceArray.add_ref(dcpath).move([stride*count, -vstride*vcount])
        count +=1
 
    
    count = 0
    vcount = 1
    
    
    for w in range(6,12):

        numCh = 7 
        constrwdth = constrictions_w[c,w]
        
        cpw = ad.CPW_single(wg_width[0], (numCh)*ydistCh[c,w], extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation)
        
        dcpath = ad.DCpath_single(widths[w-6], straight_lens[c,w], via_extend_x, via_extend_y, dcleadwidth, dctaperlength, numCh, ydistCh[c,w], dc_padlength, taper_length, add_cdt_ext, separation, extend1, extend2, ProbeSpacing, constlen, constrwdth, constrtpr)
        
        # centerx = ((cpw.xmin + cpw.xmax)/2) 
        # centery = ((cpw.ymin + cpw.ymax)/2) 
        
        centerx = ((dcpath.xmin + dcpath.xmax)/2) 
        centery = ((dcpath.ymin + dcpath.ymax)/2)
        
        dcpath.move([centerx, centery])
        cpw.move([centerx, centery-padlength-taper_length-numCh*ydistCh[c,w]/2-ch_spc/2])
        
        # move([centerx+constlen/2+constrtpr+extralength-extend1/2, centery+ampdevice.ymax/2-(numCh*straightlength+(numCh)*ydistCh)/2])

        ad.fid_squares(dcpath,50,10)

        AmpDeviceArray.add_ref(cpw).move([stride*count, -vstride*vcount])
        AmpDeviceArray.add_ref(dcpath).move([stride*count, -vstride*vcount])

        
        count +=1
    
    numperiods = 0
    count = 0
    vcount = 2
    
    for w in range(12,18):

        numCh = 1
        constrwdth = constrictions_w[c,w]
        cpw = ad.CPW_single(wg_width[0], numCh*straight_lens[c,w]+(numCh)*ydistCh[c,w], extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation) 
        
        
        dcpath = ad.DCpath_single(widths[w-12], straight_lens[c,w], via_extend_x, via_extend_y, dcleadwidth, dctaperlength, numCh, ydistCh[c,w], dc_padlength, taper_length, add_cdt_ext, separation, extend1, extend2, ProbeSpacing, constlen, constrwdth,constrtpr)
        centery = ((cpw.ymin + cpw.ymax)/2) 
        
        centerx = ((dcpath.xmin + dcpath.xmax)/2) 
        centery = ((dcpath.ymin + dcpath.ymax)/2)
        
        dcpath.move([centerx, centery])
        cpw.move([centerx, centery-padlength-taper_length-numCh*ydistCh[c,w]/2-ch_spc/2])
        
        # dcpath.move([centerx, cpw.ymax-(numCh*straight_lens[c,w]+numCh*ydistCh[c,w])])
        
        ad.fid_squares(dcpath,50,10)


        AmpDeviceArray.add_ref(cpw).move([stride*count, -(vstride-500)*vcount])
        AmpDeviceArray.add_ref(dcpath).move([stride*count, -(vstride-500)*vcount])
        


        
        count +=1
    
    
    
    
    
    centerx=((AmpDeviceArray.xmin + AmpDeviceArray.xmax)/2) 
    centery= ((AmpDeviceArray.ymin + AmpDeviceArray.ymax)/2)

    AmpDeviceArray.move([-centerx, -centery])
    
    
    # cutmarks = ad.cutmarks(50,500,die_size)
    # AmpDeviceArray.add_ref(cutmarks)

    #verify
    centerx=((AmpDeviceArray.xmin + AmpDeviceArray.xmax)/2) 
    centery= ((AmpDeviceArray.ymin + AmpDeviceArray.ymax)/2)

    #Print GDSII file and plot
    AmpDeviceArray.write_gds(output_directory + '\\Amp_chip-'+str(c+1), unit=1e-6, precision=1e-11)
    
        

pd.quickplot2(AmpDeviceArray)

exit()
    

#%% Make quarter

die_size = 10000
QuarterWafer = pd.Device('QuarterWafer')

QuarterWafer.add_ref(DCtest_device).move([-((DCtest_device.xmin + DCtest_device.xmax)/2) ,-((DCtest_device.ymin + DCtest_device.ymax)/2)])

# QuarterWafer.add_ref(AmpDeviceArray)
QuarterWafer.add_ref(AmpDeviceArray).move((die_size,0))
QuarterWafer.add_ref(AmpDeviceArray).move((0,-die_size))
QuarterWafer.add_ref(AmpDeviceArray).move((die_size,-die_size))

ad.CenterDev(QuarterWafer)
pd.quickplot(QuarterWafer)

# exit()
    

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

EBLmark = ad.EBL_marks()
wafer.add_ref(EBLmark)
wafer.add_ref(centerRef)
pd.quickplot(wafer)

wafer.write_gds(output_directory + '\\SH_Amps_wafer', unit=1e-6, precision=1e-10)






