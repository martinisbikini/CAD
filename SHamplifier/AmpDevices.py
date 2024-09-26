# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 13:21:01 2023

@author: jjw7
"""

import phidl as pd
import numpy as np
from scipy.special import ellipk
import CPW_Calc2
import phidl.utilities
import phidl.geometry as pg


#%%
# NOTE: some of these functions used these global parameters (every gap calculator will)
sub_eps = 11.9      #Real part of substrate Dielectric Constant--Si
thickness = 380     #Substrate Thickness in micron
impedance = 50      #CPW Design impedance
eps_eff = CPW_Calc2.epseff_g      #Effective epsilon calcualtion to use for calculation. CPW_Calc has several functions for multilayer substrates, conductive back plane, etc.
c_pad_length = 200  #Center Contact Pad Length. This should be shorter than ground pad length for a resonator
g_pad_length = 600  #Ground contact pad length
segments = 12       # Number of segments for exponential taper from contact pads to CPW
taper_length = 250  # Length of exponential taper
spacing = 150       #GSG contact pad spacing

#ybco_layer = pd.Layer(gds_layer = 2, name='ybcolayer', color='red')     #Called YBCO_layer for historical reasons. This layer defines the CPW gaps
gold_layer = pd.Layer(gds_layer = 1, name='goldlayer', color='yellow')  #Contact pad layer
mag_layer = pd.Layer(gds_layer = 2, name = 'filmlayer', color = 'purple')
via_layer = pd.Layer(gds_layer = 3, name = 'etchlayer', color = 'green')
wafer_layer = pd.Layer(gds_layer = 4, name = 'waferlayer', color = 'red')
mark_layer = pd.Layer(gds_layer = 5, name = 'marklayer', color = 'blue')

### Function used to define the exponential taper from Contact Pads to main CPW
def taper(w, ccw, segments):
    x_taper = []
    y_init = np.round(np.linspace(0, taper_length, segments),3)
    y_taper = y_init.tolist()
    for i in y_taper:
        x_taper.append((ccw - w)/2*np.exp(-((segments/2)*i/taper_length)) + w/2)   
    x_taper[segments -1] = w/2
    y_new = y_taper.copy()
    for i in x_taper[::-1]:
        x_taper.append(CPW_Calc2.getgap(2*i, thickness, eps_eff, impedance ,sub_eps) + i)
    #x_taper[segments*2 -1] = w/2 + gap
    for i in range(len(y_new)):
        y_taper.append(y_new[-(1+i)])
    return x_taper, y_taper

### Function used to define Linear taper from Contact Pads to main CPW
def lintaper(w, ccw, segments):
    x_taper = []
    y_init = np.round(np.linspace(0, taper_length, segments),3)
    y_taper = y_init.tolist()
    for i in y_taper:
        x_taper.append(((ccw - w)/2)*(1-i/taper_length) + w/2)   
    x_taper[segments -1] = w/2
    y_new = y_taper.copy()
    for i in x_taper[::-1]:
        x_taper.append(CPW_Calc2.getgap(2*i, thickness, eps_eff, impedance ,sub_eps) + i)
    #x_taper[segments*2 -1] = w/2 + gap
    for i in range(len(y_new)):
        y_taper.append(y_new[-(1+i)])
    return x_taper, y_taper

### define function based on taper to define outer edges of ground plane. Pass in ground-edge-to-ground edge widths. 
def groundtaper(w, ccw, segments):
    x_taper = []
    y_init = np.round(np.linspace(0, taper_length, segments),3)
    y_taper = y_init.tolist()
    for i in y_taper:
        x_taper.append((ccw - w)/2*np.exp(-((segments/2)*i/taper_length)) + w/2)   
    x_taper[segments -1] = w/2
 #   y_new = y_taper.copy()
    nx_taper = [-x for x in x_taper]
    nx_taper.reverse()
    ny_taper = y_taper.copy()
    ny_taper.reverse()
    x_taper.extend(nx_taper)
    y_taper.extend(ny_taper)
#    print(x_taper)
    return x_taper, y_taper

### define function based on taper to define outer edges of ground plane. Pass in ground-edge-to-ground edge widths. 
def lingroundtaper(w, ccw, segments):
    x_taper = []
    y_init = np.round(np.linspace(0, taper_length, segments),3)
    y_taper = y_init.tolist()
    for i in y_taper:
        x_taper.append(((ccw - w)/2)*(1-i/taper_length) + w/2)   
    x_taper[segments -1] = w/2
 #   y_new = y_taper.copy()
    nx_taper = [-x for x in x_taper]
    nx_taper.reverse()
    ny_taper = y_taper.copy()
    ny_taper.reverse()
    x_taper.extend(nx_taper)
    y_taper.extend(ny_taper)
#    print(x_taper)
    return x_taper, y_taper

def makeplane(w,ccw,segments = segments, exptaper = True):
    # makes a plane of conductor to "not" with tapered gaps made by "maketapers" 
    taperplane = pd.Device("Taperplane")
    if exptaper is True:
        xp_taper,yp_taper = groundtaper(w,ccw,segments)
    else:
        xp_taper,yp_taper = lingroundtaper(w,ccw,segments)        
    taperplane.add_polygon([tuple(xp_taper), tuple(yp_taper)], layer = gold_layer)
    return taperplane

def makeCPWseg(w, gndw, length, epseff = eps_eff,subeps = sub_eps):
    # given a center conductor width, ground-edge to ground-edge width, and length, make a cpw. maketapers should probably be redone to allow for diff eps...
    # note 11/21: all structures are currently in gold_layer
    CPWseg = pd.geometry.rectangle(size=(gndw,length),layer = gold_layer)
    CPWseg.move([-gndw/2,0])
    
    gap = CPW_Calc2.getgap(w, thickness, epseff, impedance, subeps) #note that some params are global: imp, sub_eps 
    gaps = pd.Device('gaps')
    gaprect1 = pd.geometry.rectangle(size=(gap,length), layer = gold_layer)
    gaprect2 = pd.geometry.rectangle(size=(gap,length), layer = gold_layer)
    gap1=gaps.add_ref(gaprect1)
    gap1.move([w/2,0])
    gap2=gaps.add_ref(gaprect2)
    gap2.move([-(w/2+gap),0])
    CPWseg = pd.geometry.boolean(CPWseg,gaps, 'not', layer = gold_layer)
#    pd.quickplot(CPWseg)
    return CPWseg

def makeCPWseg_fixg(w, sgndw, length, epseff = eps_eff,subeps = sub_eps):
    # wrapper for making a CPW segment, in which the gnd width is set explicitly.
    # takes width, single gnd width, length
    gap = CPW_Calc2.getgap(w, thickness, epseff, impedance, subeps) #note that some params are global: imp, sub_eps 
    totalwidth = w + 2*gap + 2*sgndw
    CPWseg = makeCPWseg(w, totalwidth, length)
    return CPWseg

def maketapergaps(w,spacing,segments = segments,exptaper = True):
    # takes width, probe spacing to taper from to make exponentially tapered gaps for CPW
    tapers = pd.Device("Tapers")
    length = 0
    pad_w, pad_g = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    if exptaper is True:
        x_taper, y_taper = taper(w, pad_w, segments)
    else:
        x_taper, y_taper = lintaper(w, pad_w, segments)
        
    x_taper[-1] = pad_w/2 + pad_g
#    for i in range(len(y_taper)):
#        y_taper[i] = y_taper[i] - (length/2 + taper_length)
    TAPER1 = pd.Device("Taper1")
    TAPER1.add_polygon([tuple(x_taper), tuple(y_taper)], layer = gold_layer)
    #TAPER2 = pd.Device("Taper2")
    for i in range(len(x_taper)):
        x_taper[i] = -1*x_taper[i]
    TAPER1.add_polygon([tuple(x_taper), tuple(y_taper)], layer = gold_layer)
#    TAPER3 = pd.geometry.copy(TAPER1).rotate(180)
    tapers.add_ref(TAPER1)
#    pd.quickplot2(tapers)
    return tapers

def maketaperCPW(nw, ngndw,spacing,widegndw,exptaper = True):
    # makes tapered CPW. uses hard-coded taper_length above :( 
    # takes width, probe spacing to taper from to make exponentially tapered gaps for CPW
    D = pd.Device('taperCPW')
    # make tapered CPW from plane and gaps
    plane = makeplane(ngndw,widegndw,segments=20,exptaper = exptaper)
    tapers = maketapergaps(nw,spacing,segments=20, exptaper = exptaper)
    D << pd.geometry.boolean(plane,tapers,'not', layer = gold_layer)
    return D

def maketaperCPW_fixg(nw, nsgndw,spacing,widesgndw,exptaper = True):
    # makes tapered CPW with defined SINGLE-SIDED ground widths(...sgndw). uses hard-coded taper_length above :( 
    # takes width, probe spacing to taper from to make exponentially or linearly tapered gaps for CPW
    # get gap sizes to calculate total widths
    bigw, biggap =CPW_Calc2.getpads(spacing, thickness, eps_eff,impedance, sub_eps)
    smallgap = CPW_Calc2.getgap(nw, thickness, eps_eff, impedance, sub_eps)
    # make total widths
    ngndw = nw + 2*smallgap + 2* nsgndw
    widegndw = bigw + 2*biggap + 2*widesgndw
    # make tapered CPW from plane and gaps
    D = maketaperCPW(nw, ngndw,spacing,widegndw, exptaper=exptaper)
    return D

def makeshort(width, height):
    # makes a short of total width width. could parametrize height
    D = pd.Device('blank')
    short = pd.geometry.rectangle(size = (width, height), layer = gold_layer)
    short.move([-width/2,0])
    D.add_ref(short)
    return D

def makemeander_seg(width, straightlength, numperiods):
    # same as regular, with double grounds between to make uniform current period
    # makes a meander using width. assumes that grounds = widths, calculates gaps
    # makes a through. Takes probe spacing, wide single-side gnd width, pad lengths
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
# make a set of arcs. this will double the grounds to make essentially parallel wires
    curveCPW = pd.Device('curveCPW')
    centgnd = pd.geometry.arc(radius = (width+gap)/2, width=width, theta = 180,layer = gold_layer)
    outgnd = pd.geometry.arc(radius = 3*(width+gap)/2, width=width, theta = 180,layer = gold_layer)
    centcond = pd.geometry.arc(radius = 5*(width+gap)/2, width = width, theta = 180,layer = gold_layer)
    curveCPW.add_ref(centgnd)
    curveCPW.add_ref(outgnd)
    curveCPW.add_ref(centcond)
# make a straight section of CPW
    straightCPW = makeCPWseg_fixg(width,width,straightlength)
# make a period, consisting of curve, straight, topcurve, straight
    oneperiod = pd.Device('oneperiod')
    bottcurve = oneperiod.add_ref(curveCPW).mirror(p1=(0,0),p2=(1,0)).movex(destination=(width+gap)/2)
    upCPW = oneperiod.add_ref(straightCPW).movex(destination=2*(width+gap))
    topcurve = oneperiod.add_ref(curveCPW).move(destination=(7*(width+gap)/2,straightlength))
    downCPW = oneperiod.add_ref(straightCPW).movex(destination=(5*(width+gap)))
    periodlen = (5+1)*(width+gap)
    oneperiod.movex(width+gap)
# make meander
    themeander = pd.Device('themeander')
    firstlen = themeander.add_ref(straightCPW)
    submeanders = []
    for step in range(numperiods):
        theshift = (step)*periodlen
        themeander.add_ref(oneperiod).movex(destination=theshift)
#    pd.quickplot(themeander)
# add bits at start and end to make import into sonnet simpler
    bridgeCPWseg = makeCPWseg_fixg(width,width,(bottcurve.ymax-bottcurve.ymin)+0.1)
    # assemble. Origins are all wacky
# make bottom connector
    topbridgesegref = themeander.add_ref(bridgeCPWseg).move([0,straightlength])
    bottbridgesegref = themeander.add_ref(bridgeCPWseg).move([numperiods*periodlen,-(bridgeCPWseg.ymax-bridgeCPWseg.ymin)])
    return(themeander)

def makemeanderthroughCPW(width, straightlength, numperiods, spacing, swidegndw, padlength,exptaperflag = True):
    # makes a meander using width. assumes that grounds = widths, calculates gaps
    # makes a through. Takes probe spacing, wide single-side gnd width, pad lengths
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
# make a set of arcs
    curveCPW = pd.Device('curveCPW')
    centgnd = pd.geometry.arc(radius = width/4, width=width/2, theta = 180,layer = gold_layer)
    outgnd = pd.geometry.arc(radius = 2*(width+gap), width=width, theta = 180,layer = gold_layer)
    centcond = pd.geometry.arc(radius = (width+gap), width = width, theta = 180,layer = gold_layer)
    curveCPW.add_ref(centgnd)
    curveCPW.add_ref(outgnd)
    curveCPW.add_ref(centcond)
# make a straight section of CPW
    straightCPW = makeCPWseg_fixg(width,width,straightlength)
# make a period, consisting of curve, straight, topcurve, straight
    oneperiod = pd.Device('oneperiod')
    bottcurve = oneperiod.add_ref(curveCPW).mirror(p1=(0,0),p2=(1,0))
    upCPW = oneperiod.add_ref(straightCPW).movex(destination=(width+gap))
    topcurve = oneperiod.add_ref(curveCPW).move(destination=(2*(width+gap),straightlength))
    downCPW = oneperiod.add_ref(straightCPW).movex(destination=(3*width+3*gap))
    periodlen = 4*width+4*gap
    oneperiod.movex(width+gap)
# make meander
    themeander = pd.Device('themeander')
    firstlen = themeander.add_ref(straightCPW)
    submeanders = []
    for step in range(numperiods):
        theshift = (step)*periodlen
        themeander.add_ref(oneperiod).movex(destination=theshift)
    #pd.quickplot2(themeander)
# make full cpw through 
    through = pd.Device('throughCPW')
    # get pads based on probe spacing
    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    # use these to make a wide CPW
    gndw = swidegndw # should make gnd metal width widegndw
    sgndwidth =width # assuming gnds are CC width
    wideCPWseg = makeCPWseg_fixg(pad_w,gndw,padlength)
    # make taper
    widesegheight = wideCPWseg.ymax
    taperCPW = maketaperCPW_fixg(width,sgndwidth,spacing, swidegndw, exptaper = exptaperflag)
    # make narrow CPW seg to extend past meander arcs
    bridgeCPWseg = makeCPWseg_fixg(width,sgndwidth,(bottcurve.ymax-bottcurve.ymin))
    # assemble. Origins are all wacky
# make bottom connector
    bottConnect = pd.Device('bottConn')
    bottpadref = bottConnect.add_ref(wideCPWseg)
    botttaperref = bottConnect.add_ref(taperCPW).move([0,widesegheight])
    bottbridgesegref = bottConnect.add_ref(bridgeCPWseg).move([0,widesegheight+(taperCPW.ymax-taperCPW.ymin)])
# now add it to the through
    bottConnref =through.add_ref(bottConnect).movex(destination=numperiods*periodlen)
    meanderCPWref = through.add_ref(themeander).move([0,bottConnect.ymax])
#    toptaperref = through.add_ref(taperCPW).rotate(180).move([0,widesegheight+2*taperCPW.ymax+(themeander.ymax-themeander.ymin)]) # account for 0 ymax on rotation of taper. 
#    toppadref = through.add_ref(wideCPWseg).move([0,widesegheight+2*taperCPW.ymax+(themeander.ymax-themeander.ymin)])
    topConnref = through.add_ref(bottConnect).mirror(p1=(0,0),p2=(1,0))
    print('bottConnect ymax='+str(bottConnect.ymax))
    topConnref.move((0,straightlength+2*bottConnect.ymax))
    textlabelref = pd.geometry.text("{} nm".format(round(width*1000, 0)), size = 20, layer = gold_layer, justify = 'center').move([0, -50])
    through.add_ref(textlabelref)
   # pd.quickplot2(through)
    return through

def makemeanderthroughCPW2(width, straightlength, numperiods, spacing, swidegndw, padlength,exptaperflag = True, textlabel = False):
    # same as regualr, with double grounds between to make uniform current period
    # makes a meander using width. assumes that grounds = widths, calculates gaps
    # makes a through. Takes probe spacing, wide single-side gnd width, pad lengths
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
# make a set of arcs. this will double the grounds to make essentially parallel wires
    curveCPW = pd.Device('curveCPW')
    centgnd = pd.geometry.arc(radius = (width+gap)/2, width=width, theta = 180,layer = gold_layer)
    outgnd = pd.geometry.arc(radius = 3*(width+gap)/2, width=width, theta = 180,layer = gold_layer)
    centcond = pd.geometry.arc(radius = 5*(width+gap)/2, width = width, theta = 180,layer = gold_layer)
    curveCPW.add_ref(centgnd)
    curveCPW.add_ref(outgnd)
    curveCPW.add_ref(centcond)
# make a straight section of CPW
    straightCPW = makeCPWseg_fixg(width,width,straightlength)
# make a period, consisting of curve, straight, topcurve, straight
    oneperiod = pd.Device('oneperiod')
    bottcurve = oneperiod.add_ref(curveCPW).mirror(p1=(0,0),p2=(1,0)).movex(destination=(width+gap)/2)
    upCPW = oneperiod.add_ref(straightCPW).movex(destination=2*(width+gap))
    topcurve = oneperiod.add_ref(curveCPW).move(destination=(7*(width+gap)/2,straightlength))
    downCPW = oneperiod.add_ref(straightCPW).movex(destination=(5*(width+gap)))
    periodlen = (5+1)*(width+gap)
    oneperiod.movex(width+gap)
# make meander
    themeander = pd.Device('themeander')
    firstlen = themeander.add_ref(straightCPW)
    submeanders = []
    for step in range(numperiods):
        theshift = (step)*periodlen
        themeander.add_ref(oneperiod).movex(destination=theshift)
   # pd.quickplot2(themeander)
# make full cpw through 
    through = pd.Device('throughCPW')
    # get pads based on probe spacing
    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    # use these to make a wide CPW
    gndw = swidegndw # should make gnd metal width widegndw
    sgndwidth =width # assuming gnds are CC width
    wideCPWseg = makeCPWseg_fixg(pad_w,gndw,padlength)
    # make taper
    widesegheight = wideCPWseg.ymax
    taperCPW = maketaperCPW_fixg(width,sgndwidth,spacing, swidegndw, exptaper = exptaperflag)
    # make narrow CPW seg to extend past meander arcs
    bridgeCPWseg = makeCPWseg_fixg(width,sgndwidth,(bottcurve.ymax-bottcurve.ymin))
    # assemble. Origins are all wacky
# make bottom connector
    bottConnect = pd.Device('bottConn')
    bottpadref = bottConnect.add_ref(wideCPWseg)
    botttaperref = bottConnect.add_ref(taperCPW).move([0,widesegheight])
    bottbridgesegref = bottConnect.add_ref(bridgeCPWseg).move([0,widesegheight+(taperCPW.ymax-taperCPW.ymin)])
# now add it to the through
    bottConnref =through.add_ref(bottConnect).movex(destination=numperiods*periodlen)
    meanderCPWref = through.add_ref(themeander).move([0,bottConnect.ymax])
#    toptaperref = through.add_ref(taperCPW).rotate(180).move([0,widesegheight+2*taperCPW.ymax+(themeander.ymax-themeander.ymin)]) # account for 0 ymax on rotation of taper. 
#    toppadref = through.add_ref(wideCPWseg).move([0,widesegheight+2*taperCPW.ymax+(themeander.ymax-themeander.ymin)])
    topConnref = through.add_ref(bottConnect).mirror(p1=(0,0),p2=(1,0))
    print('bottConnect ymax='+str(bottConnect.ymax))
    topConnref.move((0,straightlength+2*bottConnect.ymax))
    textlabelref = pd.geometry.text("{} nm".format(round(width*1000, 0)), size = 20, layer = gold_layer, justify = 'center').move([0, -50])
    if textlabel is True:
        through.add_ref(textlabelref)
   # pd.quickplot2(through)
    return through

def makemeandershortCPW(width, straightlength, numperiods, spacing, swidegndw, padlength,exptaperflag = True):
    # makes a meander using width. assumes that grounds = widths, calculates gaps
    # makes a SHORT. Takes probe spacing, wide single-side gnd width, pad lengths
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
# make a set of arcs
    curveCPW = pd.Device('curveCPW')
    centgnd = pd.geometry.arc(radius = width/4, width=width/2, theta = 180,layer = gold_layer)
    outgnd = pd.geometry.arc(radius = 2*(width+gap), width=width, theta = 180,layer = gold_layer)
    centcond = pd.geometry.arc(radius = (width+gap), width = width, theta = 180,layer = gold_layer)
    curveCPW.add_ref(centgnd)
    curveCPW.add_ref(outgnd)
    curveCPW.add_ref(centcond)
# make a straight section of CPW
    straightCPW = makeCPWseg_fixg(width,width,straightlength)
# make a period, consisting of curve, straight, topcurve, straight
    oneperiod = pd.Device('oneperiod')
    bottcurve = oneperiod.add_ref(curveCPW).mirror(p1=(0,0),p2=(1,0))
    upCPW = oneperiod.add_ref(straightCPW).movex(destination=(width+gap))
    topcurve = oneperiod.add_ref(curveCPW).move(destination=(2*(width+gap),straightlength))
    downCPW = oneperiod.add_ref(straightCPW).movex(destination=(3*width+3*gap))
    periodlen = 4*width+4*gap
    oneperiod.movex(width+gap)
# make meander
    themeander = pd.Device('themeander')
    firstlen = themeander.add_ref(straightCPW)
    submeanders = []
    for step in range(numperiods):
        theshift = (step)*periodlen
        themeander.add_ref(oneperiod).movex(destination=theshift)
   # pd.quickplot2(themeander)
# make full shorted cpw 
    short = pd.Device('shortCPW')
    # get pads based on probe spacing
    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    # use these to make a wide CPW
    gndw = swidegndw # should make gnd metal width widegndw
    sgndwidth =width # assuming gnds are CC width
    wideCPWseg = makeCPWseg_fixg(pad_w,gndw,padlength)
    # make taper
    widesegheight = wideCPWseg.ymax
    taperCPW = maketaperCPW_fixg(width,sgndwidth,spacing, swidegndw, exptaper = exptaperflag)
# make bridge segment
    bridgeCPWseg = makeCPWseg_fixg(width,sgndwidth,(bottcurve.ymax-bottcurve.ymin))
# assemble. 
# make bottom connector
    bottConnect = pd.Device('bottConn')
    bottpadref = bottConnect.add_ref(wideCPWseg)
    botttaperref = bottConnect.add_ref(taperCPW).move([0,widesegheight])
    bottbridgesegref = bottConnect.add_ref(bridgeCPWseg).move([0,widesegheight+(taperCPW.ymax-taperCPW.ymin)])
# now add it to the through
    bottConnref = short.add_ref(bottConnect).movex(destination=numperiods*periodlen)
    meanderCPWref = short.add_ref(themeander).move([0,bottConnect.ymax])
#make a shorting rect, move it into place
    topConnref = short.add_ref(pd.geometry.rectangle(size=(2*gap+3*width,3*width), layer = gold_layer)).move((-(gap+3*width/2),straightlength+bottConnect.ymax))
    textlabelref = pd.geometry.text("{} nm".format(round(width*1000, 0)), size = 20, layer = gold_layer, justify = 'center').move([0, -50])
    short.add_ref(textlabelref)
  #  pd.quickplot2(short)
    return short
       
def makemeandershortCPW2(width, straightlength, numperiods, spacing, swidegndw, padlength,exptaperflag = True, textlabel = False):
    # same as other version but with double grounds
    # makes a meander using width. assumes that grounds = widths, calculates gaps
    # makes a SHORT. Takes probe spacing, wide single-side gnd width, pad lengths
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
# make a set of arcs. this will double the grounds to make essentially parallel wires
    curveCPW = pd.Device('curveCPW')
    centgnd = pd.geometry.arc(radius = (width+gap)/2, width=width, theta = 180,layer = gold_layer)
    outgnd = pd.geometry.arc(radius = 3*(width+gap)/2, width=width, theta = 180,layer = gold_layer)
    centcond = pd.geometry.arc(radius = 5*(width+gap)/2, width = width, theta = 180,layer = gold_layer)
    curveCPW.add_ref(centgnd)
    curveCPW.add_ref(outgnd)
    curveCPW.add_ref(centcond)
# make a straight section of CPW
    straightCPW = makeCPWseg_fixg(width,width,straightlength)
# make a period, consisting of curve, straight, topcurve, straight
    oneperiod = pd.Device('oneperiod')
    bottcurve = oneperiod.add_ref(curveCPW).mirror(p1=(0,0),p2=(1,0)).movex(destination=(width+gap)/2)
    upCPW = oneperiod.add_ref(straightCPW).movex(destination=2*(width+gap))
    topcurve = oneperiod.add_ref(curveCPW).move(destination=(7*(width+gap)/2,straightlength))
    downCPW = oneperiod.add_ref(straightCPW).movex(destination=(5*(width+gap)))
    periodlen = (5+1)*(width+gap)
    oneperiod.movex(width+gap)
# make meander
    themeander = pd.Device('themeander')
    firstlen = themeander.add_ref(straightCPW)
    submeanders = []
    for step in range(numperiods):
        theshift = (step)*periodlen
        themeander.add_ref(oneperiod).movex(destination=theshift)
  #  pd.quickplot2(themeander)
# make full shorted cpw 
    short = pd.Device('shortCPW')
    # get pads based on probe spacing
    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    # use these to make a wide CPW
    gndw = swidegndw # should make gnd metal width widegndw
    sgndwidth =width # assuming gnds are CC width
    wideCPWseg = makeCPWseg_fixg(pad_w,gndw,padlength)
    # make taper
    widesegheight = wideCPWseg.ymax
    taperCPW = maketaperCPW_fixg(width,sgndwidth,spacing, swidegndw, exptaper = exptaperflag)
# make bridge segment
    bridgeCPWseg = makeCPWseg_fixg(width,sgndwidth,(bottcurve.ymax-bottcurve.ymin))
# assemble. 
# make bottom connector
    bottConnect = pd.Device('bottConn')
    bottpadref = bottConnect.add_ref(wideCPWseg)
    botttaperref = bottConnect.add_ref(taperCPW).move([0,widesegheight])
    bottbridgesegref = bottConnect.add_ref(bridgeCPWseg).move([0,widesegheight+(taperCPW.ymax-taperCPW.ymin)])
# now add it to the through
    bottConnref = short.add_ref(bottConnect).movex(destination=numperiods*periodlen)
    meanderCPWref = short.add_ref(themeander).move([0,bottConnect.ymax])
#make a shorting rect, move it into place
    topConnref = short.add_ref(pd.geometry.rectangle(size=(2*gap+3*width,3*width), layer = gold_layer)).move((-(gap+3*width/2),straightlength+bottConnect.ymax))
    textlabelref = pd.geometry.text("{} nm".format(round(width*1000, 0)), size = 20, layer = gold_layer, justify = 'center').move([0, -50])
    if textlabel is True:
        short.add_ref(textlabelref)
    #pd.quickplot2(short)
    return short

def makethroughCPW(width, gndwidth, spacing, widegndw, length, padlength,exptaperflag = True):
    # makes a through. Takes narrow width, narrow gnd width (just the metal), probe spacing, wide gnd widths, length of narrow section, pad lengths
    through = pd.Device('throughCPW')
    # get pads based on probe spacing
    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    # use these to make a wide CPW
    gndw = widegndw + pad_gap+ pad_w/2 # should make gnd metal width widegndw
    wideCPWseg = makeCPWseg(pad_w,gndw,padlength)
    # make taper
    widesegheight = wideCPWseg.ymax
    taperCPW = maketaperCPW(width,gndwidth,spacing, (widegndw+pad_gap+pad_w/2), exptaper = exptaperflag)
    # make narrow CPW
    narrowCPWseg = makeCPWseg(width,gndwidth,length)
    # assemble
    bottpadref = through.add_ref(wideCPWseg)
    botttaperref = through.add_ref(taperCPW).move([0,widesegheight])
    narrowCPWref = through.add_ref(narrowCPWseg).move([0,widesegheight+taperCPW.ymax])
    toptaperref = through.add_ref(taperCPW).rotate(180).move([0,widesegheight+2*taperCPW.ymax+(narrowCPWseg.ymax-narrowCPWseg.ymin)]) # account for 0 ymax on rotation of taper. 
    toppadref = through.add_ref(wideCPWseg).move([0,widesegheight+2*taperCPW.ymax+(narrowCPWseg.ymax-narrowCPWseg.ymin)])
    textlabelref = pd.geometry.text("{} nm".format(round(width*1000, 0)), size = 20, layer = gold_layer, justify = 'center').move([0, -50])
    through.add_ref(textlabelref)
#    pd.quickplot2(through)
    return through

def makethroughCPW_fixg(width, sgndwidth, spacing, swidegndw, length, padlength,exptaperflag = True, labelflag = True):
    # makes a through. Takes narrow width, narrow gnd width (single side metal width), probe spacing, wide single-side gnd width, length of narrow section, pad lengths
    through = pd.Device('throughCPW')
    # get pads based on probe spacing
    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    # use these to make a wide CPW
    gndw = swidegndw # should make gnd metal width widegndw
    wideCPWseg = makeCPWseg_fixg(pad_w,gndw,padlength)
    # make taper
    widesegheight = wideCPWseg.ymax
    taperCPW = maketaperCPW_fixg(width,sgndwidth,spacing, swidegndw, exptaper = exptaperflag)
    # make narrow CPW
    narrowCPWseg = makeCPWseg_fixg(width,sgndwidth,length)
    # assemble
    bottpadref = through.add_ref(wideCPWseg)
    botttaperref = through.add_ref(taperCPW).move([0,widesegheight])
    narrowCPWref = through.add_ref(narrowCPWseg).move([0,widesegheight+taperCPW.ymax])
    toptaperref = through.add_ref(taperCPW).rotate(180).move([0,widesegheight+2*taperCPW.ymax+(narrowCPWseg.ymax-narrowCPWseg.ymin)]) # account for 0 ymax on rotation of taper. 
    toppadref = through.add_ref(wideCPWseg).move([0,widesegheight+2*taperCPW.ymax+(narrowCPWseg.ymax-narrowCPWseg.ymin)])
    if labelflag is True:
        textlabelref = pd.geometry.text("{} nm".format(round(width*1000, 0)), size = 20, layer = gold_layer, justify = 'center').move([0, -50])
        through.add_ref(textlabelref)
#    pd.quickplot2(through)
    return through

def makeshortedCPW(width, gndwidth, spacing, widegndw, length, padlength,exptaperflag = True):
    # makes a Short. Takes narrow width, narrow gnd width (just the metal), probe spacing, wide gnd widths, length of narrow section, pad lengths
    theDev = pd.Device('throughCPW')
    # get pads based on probe spacing
    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    # use these to make a wide CPW
    gndw = widegndw + pad_gap+ pad_w/2 # should make gnd metal width widegndw
    wideCPWseg = makeCPWseg(pad_w,gndw,padlength)
    # make taper
    widesegheight = wideCPWseg.ymax
    taperCPW = maketaperCPW(width,gndwidth,spacing, (widegndw+pad_gap+pad_w/2), exptaper = exptaperflag)
    # make narrow CPW
    narrowCPWseg = makeCPWseg(width,gndwidth,length)
    # make short based on width of narrow CPW element
    CPWwidth = np.abs(narrowCPWseg.xmax-narrowCPWseg.xmin)
    theshort = makeshort(CPWwidth,CPWwidth)
    # assemble
    bottpadref = theDev.add_ref(wideCPWseg)
    botttaperref = theDev.add_ref(taperCPW).move([0,widesegheight])
    narrowCPWref = theDev.add_ref(narrowCPWseg).move([0,widesegheight+taperCPW.ymax])
    topshortref = theDev.add_ref(theshort).move([0,widesegheight+taperCPW.ymax+(narrowCPWseg.ymax-narrowCPWseg.ymin)]) # account for 0 ymax on rotation of taper. 
    textlabelref = pd.geometry.text("{} nm".format(round(width*1000, 0)), size = 20, layer = gold_layer, justify = 'center').move([0, -50])
    theDev.add_ref(textlabelref)
#    pd.quickplot2(theDev)
    return theDev

def makeshortedCPW_fixg(width, sgndwidth, spacing, swidegndw, length, padlength,exptaperflag = True):
    # makes a Short. Takes narrow width, single-side gnd width (just the metal), probe spacing, wide gnd widths, length of narrow section, pad lengths
    theDev = pd.Device('throughCPW')
    # get pads based on probe spacing
    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    # use these to make a wide CPW
    gndw = swidegndw
    wideCPWseg = makeCPWseg_fixg(pad_w,gndw,padlength)
    # make taper
    widesegheight = wideCPWseg.ymax
    taperCPW = maketaperCPW_fixg(width,sgndwidth,spacing, swidegndw, exptaper = exptaperflag)
    # make narrow CPW
    narrowCPWseg = makeCPWseg_fixg(width,sgndwidth,length)
    # make short based on width of narrow CPW element
    CPWwidth = np.abs(narrowCPWseg.xmax-narrowCPWseg.xmin)
    theshort = makeshort(CPWwidth,CPWwidth)
    # assemble
    bottpadref = theDev.add_ref(wideCPWseg)
    botttaperref = theDev.add_ref(taperCPW).move([0,widesegheight])
    narrowCPWref = theDev.add_ref(narrowCPWseg).move([0,widesegheight+taperCPW.ymax])
    topshortref = theDev.add_ref(theshort).move([0,widesegheight+taperCPW.ymax+(narrowCPWseg.ymax-narrowCPWseg.ymin)]) # account for 0 ymax on rotation of taper. 
    textlabelref = pd.geometry.text("{} nm".format(round(width*1000, 0)), size = 20, layer = gold_layer, justify = 'center').move([0, -50])
    theDev.add_ref(textlabelref)
#    pd.quickplot2(theDev)
    return theDev

def makeshortedsegment_fixg(width, sgndwidth, length):
    # makes a Short. Takes narrow width, single-side gnd width (just the metal), probe spacing, wide gnd widths, length of narrow section, pad lengths
    theDev = pd.Device('throughCPW')
    narrowCPWseg = makeCPWseg_fixg(width,sgndwidth,length)
    # make short based on width of narrow CPW element
    CPWwidth = np.abs(narrowCPWseg.xmax-narrowCPWseg.xmin)
    theshort = makeshort(CPWwidth,CPWwidth)
    # assemble
    narrowCPWref = theDev.add_ref(narrowCPWseg)
    topshortref = theDev.add_ref(theshort).move([0,(narrowCPWseg.ymax-narrowCPWseg.ymin)]) # account for 0 ymax on rotation of taper. 
#    textlabelref = pd.geometry.text("{} nm".format(round(width*1000, 0)), size = 20, layer = gold_layer, justify = 'center').move([0, -50])
#    theDev.add_ref(textlabelref)
 #   pd.quickplot2(theDev)
    return theDev

def taperCPW(w, spacing, taper_length, freq):
    # This is Ian's function to make his resonators. Different layering
    gap = CPW_Calc2.getgap(w, thickness, eps_eff, impedance, sub_eps)
    #print (gap)
    length_t = (CPW_Calc2.half_wavelength(gap, w, thickness, sub_eps, freq))*1e6
    print("Ideal Len:{}".format(length_t))
    length = (CPW_Calc2.half_wavelength(gap, w, thickness, sub_eps, freq))*1e6 - 2*taper_length - 2*c_pad_length
    #print (length + taper_length + c_pad_length)
    pad_w, pad_g = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    #print (pad_w)
    
    
    #Add gaps for CPW
    CPW = pd.Device("CPW")
    CPW.add_polygon([(-(w/2 + gap), -length/2), (-(w/2 + gap), length/2), (-w/2, length/2), (-w/2,-length/2)], layer = ybco_layer )
    CPW.add_polygon([((w/2 + gap), -length/2), ((w/2 + gap), length/2), (w/2, length/2), (w/2,-length/2)], layer = ybco_layer )
    
    #Add gaps for Contact Pads
    g_len = length + 2*taper_length
    BOX1 = pd.Device("box")
    BOX1.add_polygon([(-(pad_w/2 + pad_g), -g_len/2), (-(pad_w/2 + pad_g), -(g_len/2 + g_pad_length)), (-pad_w/2, -(g_len/2 + g_pad_length)), (-pad_w/2, -g_len/2)], layer = ybco_layer )
    BOX2 = pd.geometry.copy(BOX1).move([pad_w+pad_g,0])
    BOX3 = pd.geometry.copy(BOX1).move([pad_w+pad_g, g_len + g_pad_length])
    BOX4 = pd.geometry.copy(BOX1).move([0, g_len + g_pad_length])
    CPW.add_ref([BOX1, BOX2, BOX3, BOX4])
    
    #Add tapers to complete CPW
    x_taper, y_taper = taper(w, pad_w, segments)
    x_taper[-1] = pad_w/2 + pad_g
    for i in range(len(y_taper)):
        y_taper[i] = y_taper[i] - (length/2 + taper_length)
    TAPER1 = pd.Device("Taper1")
    TAPER1.add_polygon([tuple(x_taper), tuple(y_taper)], layer = ybco_layer)
    #TAPER2 = pd.Device("Taper2")
    for i in range(len(x_taper)):
        x_taper[i] = -1*x_taper[i]
    TAPER1.add_polygon([tuple(x_taper), tuple(y_taper)], layer = ybco_layer)
    TAPER3 = pd.geometry.copy(TAPER1).rotate(180)
    CPW.add_ref([TAPER1, TAPER3])
    
    
    #Remove area behind Center Contact
    CBOX1 = pd.Device("Cbox")
    b_len = g_len
    CBOX1.add_polygon([(-(pad_w/2), -(b_len/2 + c_pad_length)), (-(pad_w/2),-(b_len/2 + g_pad_length)), ((pad_w/2),-(b_len/2 + g_pad_length)),\
                       ((pad_w/2), -(b_len/2 + c_pad_length))], layer = ybco_layer)
    CBOX2 = pd.geometry.copy(CBOX1).rotate(180)
    CPW.add_ref([CBOX1, CBOX2])
    
    #Add contact pads
    GSGBOX = pd.Device("Contacts")
    GSGBOX.add_polygon([(-pad_w/2 , -b_len/2 ), (-pad_w/2, -(b_len/2 + c_pad_length)), (pad_w/2, -(b_len/2 + c_pad_length)), (pad_w/2, -b_len/2)], layer = gold_layer)
    GSGBOX.add_polygon([(-(pad_w/2+ pad_g), -b_len/2), (-(pad_w/2+ pad_g + pad_w), -b_len/2), (-(pad_w/2+ pad_g + pad_w), -(b_len/2 + g_pad_length)),\
                        (-(pad_w/2+ pad_g), -(b_len/2 + g_pad_length))], layer = gold_layer)
    GSGBOX.add_polygon([((pad_w/2+ pad_g), -b_len/2), ((pad_w/2+ pad_g + pad_w), -b_len/2), ((pad_w/2+ pad_g + pad_w), -(b_len/2 + g_pad_length)),\
                        ((pad_w/2+ pad_g), -(b_len/2 + g_pad_length))], layer = gold_layer)
    GSGBOX2 = pd.geometry.copy(GSGBOX).rotate(180)
    CPW.add_ref([GSGBOX, GSGBOX2])
    
    #Add center mark
    CM = pd.geometry.rectangle(size = (100, 20), layer = gold_layer)
    CM.move([ w + gap + 100, -10])
    CM2 = pd.geometry.copy(CM).rotate(180)
    CPW.add_ref([CM, CM2])
    
    #Add marks to contacts
    mark_sep = 100
    GM = []
    m_w = 100
    m_h = 20
    m_g = 50
    for i in range(int(g_pad_length/mark_sep + 1)):
        GM.append(pd.geometry.rectangle(size= (m_w, m_h), layer = gold_layer).move([pad_w/2 + pad_g + pad_w + m_g, length/2 + taper_length - 10 + mark_sep*i]))
        GM.append(pd.geometry.rectangle(size= (m_w, m_h), layer = gold_layer).move([-(pad_w/2 + pad_g + pad_w + m_g + m_w), length/2 + taper_length - 10 + mark_sep*i]))
        GM.append(pd.geometry.rectangle(size= (m_w, m_h), layer = gold_layer).move([pad_w/2 + pad_g + pad_w + m_g, -(length/2 + taper_length - 10 + mark_sep*i + m_h)]))
        GM.append(pd.geometry.rectangle(size= (m_w, m_h), layer = gold_layer).move([-(pad_w/2 + pad_g + pad_w + m_g + m_w), -(length/2 + taper_length - 10 + mark_sep*i + m_h)]))
    CPW.add_ref(GM)
    
    #Add CPW text
    TXT1 = pd.geometry.text("{} GHz".format(round(freq/1.0e9, 3)), size = 100, layer = gold_layer, justify = 'center').move([0, -(length/2 + taper_length + g_pad_length + 250)])
    TXT2 = pd.geometry.text("{}/{}".format(round(w, 1),round(gap, 1)), size = 100, layer = gold_layer, justify = 'center').move([0, -(length/2 + taper_length + g_pad_length + 400)])
    CPW.add_ref([TXT1, TXT2])
    
    return CPW


def CPW_pair(width, straightlength, extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation):
    mDevice = pd.Device('cpw')
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
    extralength = ((2.5*gap)+(3*width)) # this is theaightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False,  extra bit to account for the radius from the meander as it turns back around
    
    CPW1 = makemeandershortCPW2(width,straightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False, textlabel = False)
    CPW2 = makemeandershortCPW2(width, straightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False, textlabel = False).rotate(180, center = [((CPW1.xmax+CPW1.xmin)/2),CPW1.ymin]).move([separation,(2*(padlength+taper_length))+straightlength+(2*extralength)])
    
    mDevice.add_ref(CPW1)
    mDevice.add_ref(CPW2)

    return(mDevice)

def DCpath(width, straightlength, via_extend, dcleadwidth, dctaperlength, numperiods, padlength, taper_length, overlaptrim, separation, extend, ProbeSpacing):

    mDevice = pd.Device('DCpath')


    SWchan = pd.Device('SWchannel')
    xpoints = [-((separation/2)+extend+via_extend+(straightlength)), ((separation/2)+extend+via_extend+(straightlength)), ((separation/2)+via_extend+extend), -((separation/2)+via_extend+extend)]  # define the corners of the polygon
    ypoints = [-((straightlength/2)-overlaptrim), -((straightlength/2)-overlaptrim), ((straightlength/2)-overlaptrim), ((straightlength/2)-overlaptrim)]

    SWchan.add_polygon([xpoints, ypoints], layer = mag_layer)
    
    Via1 = pd.Device('via1')
    Via2 = pd.Device('via2')

    #calculate dimensions and placement locations
    via_enlarge = 0.4
    xvia = [-(straightlength/2)-via_enlarge, (straightlength/2), (straightlength/2), -(straightlength/2)-via_enlarge]
    yvia1 = [-(straightlength/2)-via_enlarge, -(straightlength/2)-via_enlarge, (straightlength/2)+via_enlarge, (straightlength/2+via_enlarge)]
    yvia2 = yvia1
    viashift = (separation/2)+extend+(straightlength/2)
    #Make via polygons
    Via1.add_polygon([xvia, yvia1], layer = via_layer)
    Via2.add_polygon([xvia, yvia2], layer = via_layer)

    # place vias
    Via1.move([-viashift, 0])
    Via2.move([viashift+via_enlarge, 0])

    Lead1 = pd.Device('lead1')
    Lead2 = pd.Device('lead2')

    extralength = 2
    dcleadlength = ProbeSpacing - (viashift - straightlength/2) + separation/2 - (taper_length+ extralength + overlaptrim+straightlength/2)
    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    xlead = [-(dcleadlength/2), (dcleadlength/2), (dcleadlength/2), -(dcleadlength/2)]
    ylead = [-(dcleadwidth/2), -(dcleadwidth/2), (dcleadwidth/2), (dcleadwidth/2)]
    xpad = [-(pad_w/1), (pad_w/1), (pad_w/1), -(pad_w/1)] #doubled the width per Matt's request
    ypad = [-(padlength/2), -(padlength/2), (padlength/2), (padlength/2)]
    
    #Assemble and place both lead+pad structures
    Lead1.add_polygon([xlead, ylead], layer = gold_layer).move([centerx - viashift + straightlength/2 -(dcleadlength/2), centery])
    Lead2.add_polygon([xlead, ylead], layer = gold_layer).move([centerx + viashift - straightlength/2 + (dcleadlength/2), centery])
    Lead1.add_polygon([xpad, ypad], layer = gold_layer).move([-ProbeSpacing, padlength/2])
    Lead2.add_polygon([xpad, ypad], layer = gold_layer).move([centerx+ProbeSpacing+separation/2, 2*centery-padlength/2])
    
    #assmeble and place the connecting curved sections of the DC leads
    
    connect1 = pd.geometry.arc(radius = taper_length + extralength + overlaptrim+ straightlength/2, width = dcleadwidth, theta = 90, layer = gold_layer)
    Lead1.add_ref(connect1).rotate(angle = 90, center = (0,0)).move([centerx - viashift + straightlength/2 -(dcleadlength), centery-(taper_length+ extralength + overlaptrim+straightlength/2)])
    Lead2.add_ref(connect1).rotate(angle = 90, center = (0,0)).move([centerx - viashift + straightlength/2 -(dcleadlength), centery-(taper_length+ extralength + overlaptrim+straightlength/2)]).rotate(angle = 180, center = (centerx, centery))
    
    Lead1.move([-dctaperlength,0])
    Lead2.move([dctaperlength,0])
    
    xleadtaper1 = [Lead1.xmax, Lead1.xmax + dctaperlength, Lead1.xmax + dctaperlength, Lead1.xmax]
    yleadtaper1 = [Lead1.ymax, SWchan.ymax, SWchan.ymin, Lead1.ymax - dcleadwidth]
    Lead1.add_polygon([xleadtaper1, yleadtaper1], layer = gold_layer)
    
    xleadtaper2 = [Lead2.xmin, Lead2.xmin - dctaperlength, Lead2.xmin - dctaperlength, Lead2.xmin]
    yleadtaper2 = [Lead2.ymin, SWchan.ymin, SWchan.ymax, Lead2.ymin + dcleadwidth]
    Lead2.add_polygon([xleadtaper2, yleadtaper2], layer = gold_layer)
        
        
    # Construct  device
    mDevice.add_ref(SWchan)
    mDevice.add_ref(Via1)
    mDevice.add_ref(Via2)
    mDevice.add_ref(Lead1)
    mDevice.add_ref(Lead2)

    ##Add Label
    textlabel = f"{numperiods}x {int(width*1000)} nm ({separation}um x {straightlength}um)"
    print(textlabel)                         
    textlabelref = pd.geometry.text(textlabel, size = 25, layer = gold_layer, justify = 'center').move([0, -50])
    mDevice.add_ref(textlabelref)
    
    #Add squares in the corners for ALL layers.  This is needed to prevent overlay issues with BEAMER
    CornerSquare = pd.Device('cornerSquare')
    
    corn_offset = 50  # offset of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
    corn_width = 10   # width of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
    
    corn_xmin = mDevice.xmin 
    corn_xmax = mDevice.xmax
    corn_ymin = mDevice.ymin
    corn_ymax = mDevice.ymax
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    
    mDevice.add_ref(CornerSquare)
    
    
    return(mDevice)
    


def AmplifierDevice_v4(width, straightlength, via_extend, extendpastswchannel,dcleadwidth,dctaperlength, numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, separation, extend, ProbeSpacing):
    mDevice = pd.Device('MyDevice')
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
    extralength = ((2.5*gap)+(3*width)) # this is theaightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False,  extra bit to account for the radius from the meander as it turns back around
    
    CPW1 = makemeandershortCPW2(width,straightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False, textlabel = False)
    CPW2 = makemeandershortCPW2(width, straightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False, textlabel = False).rotate(180, center = [((CPW1.xmax+CPW1.xmin)/2),CPW1.ymin]).move([separation,(2*(padlength+taper_length))+straightlength+(2*extralength)])
    
    mDevice.add_ref(CPW1)
    mDevice.add_ref(CPW2)

    #calculate the center of the 2 CPWs to place the SW channel (i.e., the nominal origin), should be a symmetric structure at this point
    centerx = ((mDevice.xmin + mDevice.xmax)/2) 
    centery = ((mDevice.ymin + mDevice.ymax)/2)
    
    ##### Make spinwave channel - polygon
    
    SWchan = pd.Device('SWchannel')
    #define corners of the polygon
    #for i in range(numCh):
    xpoints = [-((separation/2)+extend+via_extend+(straightlength)), ((separation/2)+extend+via_extend+(straightlength)), ((separation/2)+via_extend+extend), -((separation/2)+via_extend+extend)]  # define the corners of the polygon
    ypoints = [-((straightlength/2)-overlaptrim), -((straightlength/2)-overlaptrim), ((straightlength/2)-overlaptrim), ((straightlength/2)-overlaptrim)]
    
    #make and place spinwave channel
    #for i in range(numCh):
    SWchan.add_polygon([xpoints, ypoints], layer = mag_layer)
    SWchan.move([centerx, centery])
    
    ####### Make Vias
    Via1 = pd.Device('via1')
    Via2 = pd.Device('via2')
    
    #calculate dimensions and placement locations

    via_enlarge = 0.4
    xvia = [-(straightlength/2)-via_enlarge, (straightlength/2), (straightlength/2), -(straightlength/2)-via_enlarge]
    yvia1 = [-(straightlength/2)-via_enlarge, -(straightlength/2)-via_enlarge, (straightlength/2)+via_enlarge, (straightlength/2+via_enlarge)]
    yvia2 = yvia1
    viashift = (separation/2)+extend+(straightlength/2)
    #Make via polygons
    Via1.add_polygon([xvia, yvia1], layer = via_layer)
    Via2.add_polygon([xvia, yvia2], layer = via_layer)
    # place vias
    Via1.move([centerx-viashift, centery])
    Via2.move([centerx+viashift+via_enlarge, centery])
    
    ########  Make dc leads and contacts
    Lead1 = pd.Device('lead1')
    Lead2 = pd.Device('lead2')
    

    dcleadlength = ProbeSpacing - (viashift - straightlength/2) + separation/2 - (taper_length+ extralength + overlaptrim+straightlength/2)
    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    xlead = [-(dcleadlength/2), (dcleadlength/2), (dcleadlength/2), -(dcleadlength/2)]
    ylead = [-(dcleadwidth/2), -(dcleadwidth/2), (dcleadwidth/2), (dcleadwidth/2)]
    xpad = [-(pad_w/1), (pad_w/1), (pad_w/1), -(pad_w/1)] #doubled the width per Matt's request
    ypad = [-(padlength/2), -(padlength/2), (padlength/2), (padlength/2)]
    
    #Assemble and place both lead+pad structures
    Lead1.add_polygon([xlead, ylead], layer = gold_layer).move([centerx - viashift + straightlength/2 -(dcleadlength/2), centery])
    Lead2.add_polygon([xlead, ylead], layer = gold_layer).move([centerx + viashift - straightlength/2 + (dcleadlength/2), centery])
    Lead1.add_polygon([xpad, ypad], layer = gold_layer).move([-ProbeSpacing, padlength/2])
    Lead2.add_polygon([xpad, ypad], layer = gold_layer).move([centerx+ProbeSpacing+separation/2, 2*centery-padlength/2])
    
    #assmeble and place the connecting curved sections of the DC leads
    
    connect1 = pd.geometry.arc(radius = taper_length + extralength + overlaptrim+ straightlength/2, width = dcleadwidth, theta = 90, layer = gold_layer)
    Lead1.add_ref(connect1).rotate(angle = 90, center = (0,0)).move([centerx - viashift + straightlength/2 -(dcleadlength), centery-(taper_length+ extralength + overlaptrim+straightlength/2)])
    Lead2.add_ref(connect1).rotate(angle = 90, center = (0,0)).move([centerx - viashift + straightlength/2 -(dcleadlength), centery-(taper_length+ extralength + overlaptrim+straightlength/2)]).rotate(angle = 180, center = (centerx, centery))
    
    Lead1.move([-dctaperlength,0])
    Lead2.move([dctaperlength,0])
    
    xleadtaper1 = [Lead1.xmax, Lead1.xmax + dctaperlength, Lead1.xmax + dctaperlength, Lead1.xmax]
    yleadtaper1 = [Lead1.ymax, SWchan.ymax, SWchan.ymin, Lead1.ymax - dcleadwidth]
    Lead1.add_polygon([xleadtaper1, yleadtaper1], layer = gold_layer)
    
    
    xleadtaper2 = [Lead2.xmin, Lead2.xmin - dctaperlength, Lead2.xmin - dctaperlength, Lead2.xmin]
    yleadtaper2 = [Lead2.ymin, SWchan.ymin, SWchan.ymax, Lead2.ymin + dcleadwidth]
    Lead2.add_polygon([xleadtaper2, yleadtaper2], layer = gold_layer)
        
        
    # Construct  device
    mDevice.add_ref(SWchan)
    mDevice.add_ref(Via1)
    mDevice.add_ref(Via2)
    mDevice.add_ref(Lead1)
    mDevice.add_ref(Lead2)


    
    ##Add Label
    #textlabel2 = f"{int(width*1000)} nm {numperiods: 0.3f}x"
    textlabel = f"{numperiods}x {int(width*1000)} nm ({separation}um x {straightlength}um)"
    print(textlabel)                         
    textlabelref = pd.geometry.text(textlabel, size = 25, layer = gold_layer, justify = 'center').move([0, -50])
    mDevice.add_ref(textlabelref)
    

    #Add squares in the corners for ALL layers.  This is needed to prevent overlay issues with BEAMER
    CornerSquare = pd.Device('cornerSquare')
    
    corn_offset = 50  # offset of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
    corn_width = 10   # width of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
    
    corn_xmin = mDevice.xmin 
    corn_xmax = mDevice.xmax
    corn_ymin = mDevice.ymin
    corn_ymax = mDevice.ymax
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    
    mDevice.add_ref(CornerSquare)
    
    
    
    return(mDevice)




def SpinChannelTest_v1p0(straightlength, via_extend, extendpastswchannel,dcleadwidth,dctaperlength, numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, separation, extend, ProbeSpacing, YIG, YIG_contact_width):
    mDevice = pd.Device('MyDevice')
    
    ##### Make spinwave channel - polygon
    
    SWchan = pd.Device('SWchannel')
    #define corners of the polygon
    xpoints = [-((separation/2)+extend+via_extend+(straightlength)), ((separation/2)+extend+via_extend+(straightlength)), ((separation/2)+via_extend+extend), -((separation/2)+via_extend+extend)]  # define the corners of the polygon
    ypoints = [-((straightlength/2)-overlaptrim), -((straightlength/2)-overlaptrim), ((straightlength/2)-overlaptrim), ((straightlength/2)-overlaptrim)]
    
    #calculate the center of the 2 CPWs to place the SW channel (i.e., the nominal origin), should be a symmetric structure at this point
    centerx=((mDevice.xmin + mDevice.xmax)/2) 
    centery= ((mDevice.ymin + mDevice.ymax)/2)
    
    #make and place spinwave channel
    
    for i in range(5):
        
        
        SWchan.add_polygon([xpoints, ypoints], layer = mag_layer).move([0,-3*i])
        SWchan.move([centerx, centery])
        
        if i == 0:
            xvia_left = [SWchan.xmin-1,SWchan.xmin+1,SWchan.xmin+1,SWchan.xmin-1]
            xvia_right = [SWchan.xmax-1,SWchan.xmax+1,SWchan.xmax+1,SWchan.xmax-1]
            y_via = [SWchan.ymax+1,SWchan.ymax+1,SWchan.ymin-1,SWchan.ymin-1]
        
        SWchan.add_polygon([xvia_left,y_via], layer = via_layer).move([0,-3*i])
        SWchan.add_polygon([xvia_right,y_via], layer = via_layer).move([0,-3*i])
  
    mDevice.add_ref(SWchan)   
    # Construct  device
    
    # for i in range(10):
        
    #     mDevice.add_ref(SWchan).move([0,-3*i])
    
    centerx=((mDevice.xmin + mDevice.xmax)/2) 
    centery= ((mDevice.ymin + mDevice.ymax)/2)
    
    mDevice.move([-centerx, -centery])
    ##Add Label
    #textlabel2 = f"{int(width*1000)} nm {numperiods: 0.3f}x"
    textlabel = f"{np.round(straightlength,decimals = 3)}um"
    print(textlabel)                    
    textlabelref = pd.geometry.text(textlabel, size = 50, layer = mark_layer, justify = 'center').move([-50,0])
    mDevice.add_ref(textlabelref).move([SWchan.xmin+textlabelref.xmin,-SWchan.ymax/2-textlabelref.ymax/2])
    

    
    
    # Add marks for the spinwave channels to find in SEM
    
   
    
    #Add squares in the corners for ALL layers.  This is needed to prevent overlay issues with BEAMER
    CornerSquare = pd.Device('cornerSquare')
    
    corn_offset = 50  # offset of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
    corn_width = 10   # width of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
    
    corn_xmin = mDevice.xmin 
    corn_xmax = mDevice.xmax
    corn_ymin = mDevice.ymin
    corn_ymax = mDevice.ymax
    
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mark_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mark_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mark_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mark_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    

    
    mDevice.add_ref(CornerSquare)
    
    centerx=((mDevice.xmin + mDevice.xmax)/2) 
    centery= ((mDevice.ymin + mDevice.ymax)/2)
    
    mDevice.move([-centerx, -centery])
    
    return(mDevice)


def kdepdamping_v0(width, straightlength, via_extend, extendpastswchannel,dcleadwidth,dctaperlength, numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, separation, extend, ProbeSpacing, YIG, YIG_contact_width):
    mDevice = pd.Device('MyDevice')
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
    extralength = ((2.5*gap)+(3*width)) # this is theaightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False,  extra bit to account for the radius from the meander as it turns back around
    
    CPW1 = makemeandershortCPW2(width,straightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False, textlabel = False)
    CPW2 = makemeandershortCPW2(width, straightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False, textlabel = False).rotate(180, center = [((CPW1.xmax+CPW1.xmin)/2),CPW1.ymin]).move([separation,(2*(padlength+taper_length))+straightlength+(2*extralength)])
    
    mDevice.add_ref(CPW1)
    mDevice.add_ref(CPW2)
    
    ##### Make spinwave channel - polygon
    
    SWchan = pd.Device('SWchannel')
    #define corners of the polygon
    xpoints = [-((separation/2)+extend+via_extend+(straightlength)), ((separation/2)+extend+via_extend+(straightlength)), ((separation/2)+extend+via_extend+(straightlength)), -((separation/2)+extend+via_extend+(straightlength))]  # define the corners of the polygon
    ypoints = [-((straightlength/2)-overlaptrim), -((straightlength/2)-overlaptrim), ((straightlength/2)-overlaptrim), ((straightlength/2)-overlaptrim)]
    
    #calculate the center of the 2 CPWs to place the SW channel (i.e., the nominal origin), should be a symmetric structure at this point
    centerx=((mDevice.xmin + mDevice.xmax)/2) 
    centery= ((mDevice.ymin + mDevice.ymax)/2)
    
    #make and place spinwave channel
    SWchan.add_polygon([xpoints, ypoints], layer = mag_layer)
    SWchan.move([centerx, centery])
    
    ####### Make Vias
    Via1 = pd.Device('via1')
    Via2 = pd.Device('via2')
    
    #calculate dimensions and placement locations
    if YIG is True:
        xvia = [-(separation/2-YIG_contact_width/2), (separation/2-YIG_contact_width/2), (separation/2-YIG_contact_width/2), -(separation/2-YIG_contact_width/2)]    
        yvia1 = [-(straightlength/2), -(straightlength/2), (straightlength/2), (straightlength/2)]
        viashift = 0  # keeping the separation between the CPW and contact the same as YIG_contact_width for now to simplify
        #Make via polygons
        Via1.add_polygon([xvia, yvia1], layer = via_layer)
        # place vias
        Via1.move([centerx-viashift, centery])
            
    else:
        via_enlarge = 0.4
        xvia = [-(straightlength/2)-via_enlarge, (straightlength/2), (straightlength/2), -(straightlength/2)-via_enlarge]
        yvia1 = [-(straightlength/2)-via_enlarge, -(straightlength/2)-via_enlarge, (straightlength/2)+via_enlarge, (straightlength/2+via_enlarge)]
        yvia2 = yvia1
        viashift = (separation/2)+extend+(straightlength/2)
        #Make via polygons
        Via1.add_polygon([xvia, yvia1], layer = via_layer)
        Via2.add_polygon([xvia, yvia2], layer = via_layer)
        # place vias
        Via1.move([centerx-viashift, centery])
        Via2.move([centerx+viashift+via_enlarge, centery])
    
    ########  Make dc leads and contacts
    Lead1 = pd.Device('lead1')
    Lead2 = pd.Device('lead2')
    
    if YIG is True:
        #Extra bit to go down to the Pt in between the CPWs.  The use of 'via' is leftover from legacy since this was copied from the via script
        xvia = [-(YIG_contact_width/2), (YIG_contact_width/2), (YIG_contact_width/2), -(YIG_contact_width/2)]    
        yvia1 = [-(straightlength/2), -(straightlength/2), (straightlength/0.5), (straightlength/0.5)]
        yvia2 = [-(straightlength/0.5), -(straightlength/0.5), (straightlength/2), (straightlength/2)]    
        viashift = (separation/2)-YIG_contact_width  # keeping the separation between the CPW and contact the same as YIG_contact_width for now to simplify
        #Make  polygons
        Lead1.add_polygon([xvia, yvia1], layer = gold_layer)
        Lead2.add_polygon([xvia, yvia2], layer = gold_layer)
        # place items
        Lead1.move([centerx-viashift, centery])
        Lead2.move([centerx+viashift, centery])
        
        fudge = 6 #fudge factor - for some reason I cannot get the domension to work out Grrrrr
        R = taper_length + extralength + overlaptrim + straightlength/2 # radius of curved section
        dcleadlength = ProbeSpacing - R  + 2*YIG_contact_width 
        pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
        xlead = [-(dcleadlength/2), (dcleadlength/2)+ extralength+ 2*YIG_contact_width-fudge, (dcleadlength/2)+ extralength+ 2*YIG_contact_width-fudge, -(dcleadlength/2)]
        ylead = [-(straightlength/2), -(straightlength/2), (straightlength/2), (straightlength/2)]
        xpad = [-(pad_w/1), (pad_w/1), (pad_w/1), -(pad_w/1)] #doubled the width per Matt's request
        ypad = [-(padlength/2), -(padlength/2), (padlength/2)+ (1.5*straightlength), (padlength/2) + (1.5*straightlength)]
        
     
        
        #Assemble and place both lead+pad structures
        Lead1.add_polygon([xlead, ylead], layer = gold_layer).move([-dcleadlength/2 + YIG_contact_width, centery + (1.5*straightlength)])
        Lead2.add_polygon([xlead, ylead], layer = gold_layer).move([-dcleadlength/2 + YIG_contact_width, centery + (1.5*straightlength)]).rotate(angle = 180, center = (centerx, centery))
        Lead1.add_polygon([xpad, ypad], layer = gold_layer).move([-ProbeSpacing, padlength/2])
        Lead2.add_polygon([xpad, ypad], layer = gold_layer).move([-ProbeSpacing, padlength/2]).rotate(angle = 180, center = (centerx, centery))
        
        #assmeble and place the connecting curved sections of the DC leads
        connect1 = pd.geometry.arc(radius = R, width = straightlength, theta = 90, layer = gold_layer)
        Lead1.add_ref(connect1).rotate(angle = 90, center = (0,0)).move([-dcleadlength + YIG_contact_width, centery-(taper_length+ extralength + overlaptrim+straightlength/2) + (1.5*straightlength)])
        Lead1.add_ref(connect1).rotate(angle = 90, center = (0,0)).move([-dcleadlength + YIG_contact_width, centery-(taper_length+ extralength + overlaptrim+straightlength/2) + (1.5*straightlength)]).rotate(angle = 180, center = (centerx, centery))    
        
            
    else:
        dcleadlength = ProbeSpacing - (viashift - straightlength/2) + separation/2 - (taper_length+ extralength + overlaptrim+straightlength/2)
        pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
        xlead = [-(dcleadlength/2), (dcleadlength/2), (dcleadlength/2), -(dcleadlength/2)]
        ylead = [-(dcleadwidth/2), -(dcleadwidth/2), (dcleadwidth/2), (dcleadwidth/2)]
        xpad = [-(pad_w/1), (pad_w/1), (pad_w/1), -(pad_w/1)] #doubled the width per Matt's request
        ypad = [-(padlength/2), -(padlength/2), (padlength/2), (padlength/2)]
        
        #Assemble and place both lead+pad structures
        Lead1.add_polygon([xlead, ylead], layer = gold_layer).move([centerx - viashift + straightlength/2 -(dcleadlength/2), centery])
        Lead2.add_polygon([xlead, ylead], layer = gold_layer).move([centerx + viashift - straightlength/2 + (dcleadlength/2), centery])
        Lead1.add_polygon([xpad, ypad], layer = gold_layer).move([-ProbeSpacing, padlength/2])
        Lead2.add_polygon([xpad, ypad], layer = gold_layer).move([centerx+ProbeSpacing+separation/2, 2*centery-padlength/2])
        
        #assmeble and place the connecting curved sections of the DC leads
        
        connect1 = pd.geometry.arc(radius = taper_length + extralength + overlaptrim+ straightlength/2, width = dcleadwidth, theta = 90, layer = gold_layer)
        Lead1.add_ref(connect1).rotate(angle = 90, center = (0,0)).move([centerx - viashift + straightlength/2 -(dcleadlength), centery-(taper_length+ extralength + overlaptrim+straightlength/2)])
        Lead2.add_ref(connect1).rotate(angle = 90, center = (0,0)).move([centerx - viashift + straightlength/2 -(dcleadlength), centery-(taper_length+ extralength + overlaptrim+straightlength/2)]).rotate(angle = 180, center = (centerx, centery))
        
        Lead1.move([-dctaperlength,0])
        Lead2.move([dctaperlength,0])
        
        
        xleadtaper1 = [Lead1.xmax, Lead1.xmax + dctaperlength, Lead1.xmax + dctaperlength, Lead1.xmax]
        yleadtaper1 = [Lead1.ymax, SWchan.ymax, SWchan.ymin, Lead1.ymax - dcleadwidth]
        Lead1.add_polygon([xleadtaper1, yleadtaper1], layer = gold_layer)
        
        
        xleadtaper2 = [Lead2.xmin, Lead2.xmin - dctaperlength, Lead2.xmin - dctaperlength, Lead2.xmin]
        yleadtaper2 = [Lead2.ymin, SWchan.ymin, SWchan.ymax, Lead2.ymin + dcleadwidth]
        Lead2.add_polygon([xleadtaper2, yleadtaper2], layer = gold_layer)
        
        
    # Construct  device
    mDevice.add_ref(SWchan)
    #mDevice.add_ref(Via1)
    #mDevice.add_ref(Via2)
    #mDevice.add_ref(Lead1)
    #mDevice.add_ref(Lead2)
    
    ##Add Label
    #textlabel2 = f"{int(width*1000)} nm {numperiods: 0.3f}x"
    textlabel = f"{numperiods}x {int(width*1000)} nm ({separation}um x {straightlength}um)"
    print(textlabel)                         
    textlabelref = pd.geometry.text(textlabel, size = 25, layer = gold_layer, justify = 'center').move([0, -50])
    mDevice.add_ref(textlabelref)
    
    #Add squares in the corners for ALL layers.  This is needed to prevent overlay issues with BEAMER
    CornerSquare = pd.Device('cornerSquare')
    
    corn_offset = 50  # offset of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
    corn_width = 10   # width of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
    
    corn_xmin = mDevice.xmin 
    corn_xmax = mDevice.xmax
    corn_ymin = mDevice.ymin
    corn_ymax = mDevice.ymax
    
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    

    
    mDevice.add_ref(CornerSquare)
    
    
    
    return(mDevice)











################### Devices with just ground and signal ########################################



def kdepdamping_gs(width, straightlength, via_extend, extendpastswchannel,dcleadwidth,dctaperlength, numperiods, GSGspacing, swidegndw, padlength, taper_length, overlaptrim, separation, extend, ProbeSpacing, YIG, YIG_contact_width):
    mDevice = pd.Device('MyDevice')
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
    extralength = ((2.5*gap)+(3*width)) # this is theaightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False,  extra bit to account for the radius from the meander as it turns back around
    
    CPW1 = makemeandershortCPW2_gs(width,straightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False, textlabel = False)
    CPW2 = makemeandershortCPW2_gs(width, straightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False, textlabel = False).rotate(180, center = [((CPW1.xmax+CPW1.xmin)/2),CPW1.ymin]).move([separation,(2*(padlength+taper_length))+straightlength+(2*extralength)])
    
    mDevice.add_ref(CPW1)
    mDevice.add_ref(CPW2)
    
    ##### Make spinwave channel - polygon
    
    SWchan = pd.Device('SWchannel')
    #define corners of the polygon
    # xpoints = [-((separation/2)+extend+via_extend+(straightlength)), ((separation/2)+extend+via_extend+(straightlength)), ((separation/2)+extend+via_extend+(straightlength)), -((separation/2)+extend+via_extend+(straightlength))]  # define the corners of the polygon
    # ypoints = [-((straightlength/2)-overlaptrim), -((straightlength/2)-overlaptrim), ((straightlength/2)-overlaptrim), ((straightlength/2)-overlaptrim)]
    
    xpoints = [-((separation/2)+extend+via_extend-width), ((separation/2)+extend+via_extend-width), ((separation/2)+extend+via_extend-width), -((separation/2)+extend+via_extend-width)]  # define the corners of the polygon
    ypoints = [-((straightlength/2)-overlaptrim), -((straightlength/2)-overlaptrim), ((straightlength/2)-overlaptrim), ((straightlength/2)-overlaptrim)]
    
    
    #calculate the center of the 2 CPWs to place the SW channel (i.e., the nominal origin), should be a symmetric structure at this point
    centerx=((mDevice.xmin + mDevice.xmax)/2) 
    centery= ((mDevice.ymin + mDevice.ymax)/2)
    
    #make and place spinwave channel
    SWchan.add_polygon([xpoints, ypoints], layer = mag_layer)
    SWchan.move([centerx, centery])
    
    mDevice.add_ref(SWchan)
    
    
    
    textlabel = f"{numperiods}x {int(width*1000)} nm ({separation}um x {straightlength}um)"
    print(textlabel)                         
    textlabelref = pd.geometry.text(textlabel, size = 25, layer = gold_layer, justify = 'center').move([0, -50])
    mDevice.add_ref(textlabelref)
    
    #Add squares in the corners for ALL layers.  This is needed to prevent overlay issues with BEAMER
    CornerSquare = pd.Device('cornerSquare')
    
    corn_offset = 50  # offset of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
    corn_width = 10   # width of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
    
    corn_xmin = mDevice.xmin 
    corn_xmax = mDevice.xmax
    corn_ymin = mDevice.ymin
    corn_ymax = mDevice.ymax
    
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = via_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mag_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    

    
    mDevice.add_ref(CornerSquare)
    
    
    
    return(mDevice)



def makemeandershortCPW2_gs(width, straightlength, numperiods, spacing, swidegndw, padlength,exptaperflag = True, textlabel = True):
    # same as other version but with double grounds
    # makes a meander using width. assumes that grounds = widths, calculates gaps
    # makes a SHORT. Takes probe spacing, wide single-side gnd width, pad lengths
    
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
# make a set of arcs. this will double the grounds to make essentially parallel wires
    curveCPW = pd.Device('curveCPW')
    centgnd = pd.geometry.arc(radius = (width+gap)/2, width=width, theta = 180,layer = gold_layer)
    outgnd = pd.geometry.arc(radius = 3*(width+gap)/2, width=width, theta = 180,layer = gold_layer)
    centcond = pd.geometry.arc(radius = 5*(width+gap)/2, width = width, theta = 180,layer = gold_layer)
    curveCPW.add_ref(centgnd)
    curveCPW.add_ref(outgnd)
    curveCPW.add_ref(centcond)
# make a straight section of CPW
    straightCPW = makeCPWseg_fixg_gs_fixgap(width,width,straightlength)
# make a period, consisting of curve, straight, topcurve, straight
    oneperiod = pd.Device('oneperiod')
    bottcurve = oneperiod.add_ref(curveCPW).mirror(p1=(0,0),p2=(1,0)).movex(destination=(width+gap)/2)
    upCPW = oneperiod.add_ref(straightCPW).movex(destination=2*(width+gap))
    topcurve = oneperiod.add_ref(curveCPW).move(destination=(7*(width+gap)/2,straightlength))
    downCPW = oneperiod.add_ref(straightCPW).movex(destination=(5*(width+gap)))
    periodlen = (5+1)*(width+gap)
    oneperiod.movex(width+gap)
# make meander
    themeander = pd.Device('themeander')
    firstlen = themeander.add_ref(straightCPW)
    submeanders = []
    for step in range(numperiods):
        theshift = (step)*periodlen
        themeander.add_ref(oneperiod).movex(destination=theshift)
  #  pd.quickplot2(themeander)
# make full shorted cpw 
    short = pd.Device('shortCPW')
    # get pads based on probe spacing
    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    # use these to make a wide CPW
    gndw = swidegndw # should make gnd metal width widegndw
    sgndwidth =width # assuming gnds are CC width
    wideCPWseg = makeCPWseg_fixg_gs(pad_w,gndw,padlength)
    #pd.quickplot2(wideCPWseg)
    # make taper
    widesegheight = wideCPWseg.ymax
    #taperCPW = maketaperCPW_fixg_gs(width,sgndwidth,spacing, swidegndw, exptaper = exptaperflag)
# make bridge segment
    bridgeCPWseg = makeCPWseg_fixg_gs_fixgap(width,sgndwidth,(bottcurve.ymax-bottcurve.ymin))
    #bridgeCPWseg = makeCPWseg_fixg_gs_fixgap(width,sgndwidth,(3))
# assemble. 
# make bottom connector
    bottConnect = pd.Device('bottConn')
    bottpadref = bottConnect.add_ref(wideCPWseg)
    
    #botttaperref = bottConnect.add_ref(taperCPW).move([0,widesegheight])
    bottbridgesegref = bottConnect.add_ref(bridgeCPWseg).move([-19,widesegheight+(taper_length)])
 
    taper1_xpts = [bottConnect.xmin,bottbridgesegref.xmin,bottbridgesegref.xmin+20+2*width,bottConnect.xmin+250]
    taper1_ypts = [widesegheight,widesegheight+taper_length,widesegheight+taper_length,widesegheight]
    
    bottConnect.add_polygon([taper1_xpts,taper1_ypts],layer = gold_layer)


    taper2_xpts = [bottConnect.xmax,bottbridgesegref.xmax,bottbridgesegref.xmax-width,bottConnect.xmax-100]
    taper2_ypts = [widesegheight,widesegheight+taper_length,widesegheight+taper_length,widesegheight]
    
    bottConnect.add_polygon([taper2_xpts,taper2_ypts],layer = gold_layer)

#now add it to the through
    bottConnref = short.add_ref(bottConnect).movex(destination=numperiods*periodlen)
    meanderCPWref = short.add_ref(themeander).move([-19,bottConnect.ymax])
#make a shorting rect, move it into place
    #topConnref = short.add_ref(pd.geometry.rectangle(size=(2*gap+3*width,3*width), layer = gold_layer)).move((-(gap+3*width/2),straightlength+bottConnect.ymax))
    topConnref = short.add_ref(pd.geometry.rectangle(size=(2*20+3*width,20*width), layer = gold_layer)).move((-(20+3*width/2)-19,straightlength+bottConnect.ymax))
    
    
    
    textlabelref = pd.geometry.text("{} nm".format(round(width*1000, 0)), size = 20, layer = gold_layer, justify = 'center').move([0, -50])
    if textlabel is True:
        short.add_ref(textlabelref)
   # pd.quickplot2(short)
    return short


def makeCPWseg_gs(w, gndw, length, epseff = eps_eff,subeps = sub_eps):
    # given a center conductor width, ground-edge to ground-edge width, and length, make a cpw. maketapers should probably be redone to allow for diff eps...
    # note 11/21: all structures are currently in gold_layer
    CPWseg = pd.geometry.rectangle(size=(gndw,length),layer = gold_layer)
    CPWseg.move([-gndw/2,0])
    
    gap = CPW_Calc2.getgap(w, thickness, epseff, impedance, subeps) #note that some params are global: imp, sub_eps 
    gaps = pd.Device('gaps')
    gaprect1 = pd.geometry.rectangle(size=(gap,length), layer = gold_layer)
    #gaprect2 = pd.geometry.rectangle(size=(gap,length), layer = gold_layer)
    gap1=gaps.add_ref(gaprect1)
    gap1.move([w/2,0])
    #gap2=gaps.add_ref(gaprect2)
    #gap2.move([-(w/2+gap),0])
    CPWseg = pd.geometry.boolean(CPWseg,gaps, 'not', layer = gold_layer)
#    pd.quickplot(CPWseg)
    return CPWseg


def makeCPWseg_gs_fixgap(w, gndw, length, epseff = eps_eff,subeps = sub_eps):
    # given a center conductor width, ground-edge to ground-edge width, and length, make a cpw. maketapers should probably be redone to allow for diff eps...
    # note 11/21: all structures are currently in gold_layer
    CPWseg = pd.geometry.rectangle(size=(gndw,length),layer = gold_layer)
    CPWseg.move([-gndw/2,0])
    
    gap = 20
    #gap = CPW_Calc2.getgap(w, thickness, epseff, impedance, subeps) #note that some params are global: imp, sub_eps 
    gaps = pd.Device('gaps')
    gaprect1 = pd.geometry.rectangle(size=(gap,length), layer = gold_layer)
    #gaprect2 = pd.geometry.rectangle(size=(gap,length), layer = gold_layer)
    gap1=gaps.add_ref(gaprect1)
    gap1.move([w/2,0])
    #gap2=gaps.add_ref(gaprect2)
    #gap2.move([-(w/2+gap),0])
    CPWseg = pd.geometry.boolean(CPWseg,gaps, 'not', layer = gold_layer)
#    pd.quickplot(CPWseg)
    return CPWseg

def makeCPWseg_fixg_gs(w, sgndw, length, epseff = eps_eff,subeps = sub_eps):
    # wrapper for making a CPW segment, in which the gnd width is set explicitly.
    # takes width, single gnd width, length
    
    
    gap = CPW_Calc2.getgap(w, thickness, epseff, impedance, subeps) #note that some params are global: imp, sub_eps 
    totalwidth = w + 2*gap + 2*sgndw
    CPWseg = makeCPWseg_gs(w, totalwidth, length)
    return CPWseg

def makeCPWseg_fixg_gs_fixgap(w, sgndw, length, epseff = eps_eff,subeps = sub_eps):
    # wrapper for making a CPW segment, in which the gnd width is set explicitly.
    # takes width, single gnd width, length
    
    
    gap = 20
    #gap = CPW_Calc2.getgap(w, thickness, epseff, impedance, subeps) #note that some params are global: imp, sub_eps 
    totalwidth = w + 2*gap + 2*sgndw
    CPWseg = makeCPWseg_gs_fixgap(w, totalwidth, length)
    return CPWseg




def maketapergaps_gs(w,spacing,segments = segments,exptaper = True):
    # takes width, probe spacing to taper from to make exponentially tapered gaps for CPW
    tapers = pd.Device("Tapers")
    length = 0
    pad_w, pad_g = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    if exptaper is True:
        x_taper, y_taper = taper(w, pad_w, segments)
    else:
        x_taper, y_taper = lintaper_gs(w, pad_w, segments)
        
    x_taper[-1] = pad_w/2 + pad_g
#    for i in range(len(y_taper)):
#        y_taper[i] = y_taper[i] - (length/2 + taper_length)
    TAPER1 = pd.Device("Taper1")
    TAPER1.add_polygon([tuple(x_taper), tuple(y_taper)], layer = gold_layer)
    #TAPER2 = pd.Device("Taper2")
    for i in range(len(x_taper)):
        x_taper[i] = -1*x_taper[i]
    #TAPER1.add_polygon([tuple(x_taper), tuple(y_taper)], layer = gold_layer)
#    TAPER3 = pd.geometry.copy(TAPER1).rotate(180)
    tapers.add_ref(TAPER1)
#    pd.quickplot2(tapers)
    return tapers

def maketaperCPW_gs(nw, ngndw,spacing,widegndw,exptaper = True):
    # makes tapered CPW. uses hard-coded taper_length above :( 
    # takes width, probe spacing to taper from to make exponentially tapered gaps for CPW
    D = pd.Device('taperCPW')
    # make tapered CPW from plane and gaps
    plane = makeplane(ngndw,widegndw,segments=20,exptaper = exptaper)
    tapers = maketapergaps_gs(nw,spacing,segments=20, exptaper = exptaper)
    D << pd.geometry.boolean(plane,tapers,'not', layer = gold_layer)
    return D

def maketaperCPW_fixg_gs(nw, nsgndw,spacing,widesgndw,exptaper = True):
    # makes tapered CPW with defined SINGLE-SIDED ground widths(...sgndw). uses hard-coded taper_length above :( 
    # takes width, probe spacing to taper from to make exponentially or linearly tapered gaps for CPW
    # get gap sizes to calculate total widths
    bigw, biggap =CPW_Calc2.getpads(spacing, thickness, eps_eff,impedance, sub_eps)
    smallgap = CPW_Calc2.getgap(nw, thickness, eps_eff, impedance, sub_eps)
    # make total widths
    ngndw = nw + 2*smallgap + 2* nsgndw
    widegndw = bigw + 2*biggap + 2*widesgndw
    # make tapered CPW from plane and gaps
    D = maketaperCPW_gs(nw, ngndw,spacing,widegndw, exptaper=exptaper)
    return D

def makeshort_gs(width, height):
    # makes a short of total width width. could parametrize height
    D = pd.Device('blank')
    short = pd.geometry.rectangle(size = (width, height), layer = gold_layer)
    short.move([-width/2,0])
    D.add_ref(short)
    return D

def lintaper_gs(w, ccw, segments):
    x_taper = []
    y_init = np.round(np.linspace(0, taper_length, segments),3)
    y_taper = y_init.tolist()
    for i in y_taper:
        x_taper.append(((ccw - w)/2)*(1-i/taper_length) + w/2)   
    x_taper[segments -1] = w/2
    y_new = y_taper.copy()
    for i in x_taper[::-1]:
        x_taper.append(CPW_Calc2.getgap(2*i, thickness, eps_eff, impedance ,sub_eps) + i)
        
    #x_taper[segments*2 -1] = w/2 + gap
    for i in range(len(y_new)):
        y_taper.append(y_new[-(1+i)])
    return x_taper, y_taper



#### Testing CPWs with different PMMA baking ##########

def CPWtest_v0(width, length):
    
    testCPW = pd.Device('testCPW')
    temp = makeCPWseg_fixg(width,width, length)
    gegewidth = temp.xmax - temp.xmin
    sh = makeshort(gegewidth,gegewidth)
    testCPW.add_ref(sh).move([0,temp.ymax-temp.ymin])
    testCPW.add_ref(temp)
    
    return testCPW
    
def CPWtest_v0_array(width,length,num_wgs):
    
    FinalDev = pd.Device('FinalDev')
    CPWarr = pd.Device('testCPWarr')
    temp = CPWtest_v0(width,length)
    spacing = 3*(temp.xmax-temp.xmin)
    
    for n in range(num_wgs):
        
        CPWarr.add_ref(temp).move([spacing*n,0])
    
    
    centerx=((CPWarr.xmin + CPWarr.xmax)/2) 
    centery= ((CPWarr.ymin + CPWarr.ymax)/2)
    
    CPWarr.move([-centerx, -centery])
    FinalDev.add_ref(CPWarr)
    ##Add Label
    #textlabel2 = f"{int(width*1000)} nm {numperiods: 0.3f}x"
    textlabel = f"{np.round(width,decimals = 3)}um"
    print(textlabel)                    
    textlabelref = pd.geometry.text(textlabel, size = 50, layer = mark_layer, justify = 'center').move([-20,0])
    
    FinalDev.add_ref(textlabelref).move([-(textlabelref.xmax-textlabelref.xmin)/2-(CPWarr.xmax-CPWarr.xmin)/2,-(textlabelref.ymax-textlabelref.ymin)/2])



    #Add squares in the corners for ALL layers.  This is needed to prevent overlay issues with BEAMER
    CornerSquare = pd.Device('cornerSquare')
    
    corn_offset = 50  # offset of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
    corn_width = 10   # width of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
    
    corn_xmin = FinalDev.xmin 
    corn_xmax = FinalDev.xmax
    corn_ymin = FinalDev.ymin
    corn_ymax = FinalDev.ymax
    
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mark_layer). move([corn_xmin-corn_offset, corn_ymin-corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mark_layer). move([corn_xmax+corn_offset, corn_ymin-corn_offset])
    
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mark_layer). move([corn_xmax+corn_offset, corn_ymax+corn_offset])
    
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = gold_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    CornerSquare.add_polygon([(-corn_width/2, corn_width/2, corn_width/2, -corn_width/2), (corn_width/2, corn_width/2, -corn_width/2, -corn_width/2)], layer = mark_layer). move([corn_xmin-corn_offset, corn_ymax+corn_offset])
    

    
    FinalDev.add_ref(CornerSquare)
    
    centerx=((FinalDev.xmin + FinalDev.xmax)/2) 
    centery= ((FinalDev.ymin + FinalDev.ymax)/2)
    
    FinalDev.move([-centerx, -centery])
    
    return FinalDev


###### HMOMM Devices

def makecross(width, length):
    
    cross1 = pd.Device('HMOMMDev')
    cross2 = pd.Device('HMOMMDev')
    
    cross1.add_polygon([(-length/2,-length/2,length/2,length/2),(-width/2,width/2,width/2,-width/2)],layer = gold_layer)
    cross2.add_polygon([(-width/2,width/2,width/2,-width/2),(-length/2,-length/2,length/2,length/2)],layer = gold_layer)
    
    cross  = pg.boolean(A = cross1, B = cross2, operation = 'or',  precision = 1e-6,
               num_divisions = [1,1], layer = gold_layer)
    return cross

def vertrec(width, length):
    
    cross2 = pd.Device('HMOMMDev')
    
    cross2.add_polygon([(-width/2,width/2,width/2,-width/2),(-length/2,-length/2,length/2,length/2)],layer = gold_layer)
    
    return cross2


def makeHMOMMDot_v0(width, straightlength, numperiods, spacing, swidegndw, padlength,dot_diam,numdots,exptaperflag = False, textlabel = True):
    
    HMOMMDev = pd.Device('HMOMMDev')
    
    # make the CPW structure #
    CPWtaper = makemeandershortCPW2(width, straightlength, numperiods, spacing, swidegndw, padlength,exptaperflag , textlabel)
    
    # get dimensions to place crosses/dots
    narrowgap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
    crosslen = narrowgap/2
    crosswid = crosslen/100
    
    
    fidmark = makecross(crosswid,crosslen)
    dot = pg.circle(dot_diam/2)
    fidmark_arr = pd.Device('fidmarkarray')
    #fidmark_arr.add_ref(dot)
    fidmark_arr.add_array(fidmark,rows = numdots+1, columns = 4, spacing =(width+narrowgap,straightlength/(numdots+1)))
    fidmark_arr.add_array(dot,rows = numdots, columns = 4, spacing =(width+narrowgap,straightlength/(numdots+1))).move((0,(straightlength/(numdots+1))/2))
    
    center_arr_x = (fidmark_arr.xmax + fidmark_arr.xmin) / 2
    center_arr_y = (fidmark_arr.ymax + fidmark_arr.ymin) / 2
    fidmark_arr.move((-center_arr_x,-center_arr_y))
    
    
    HMOMMDev.add_ref(CPWtaper)
    HMOMMDev.add_ref(fidmark_arr).move((0,CPWtaper.ymax-fidmark_arr.ymax-4*width))
    
    return HMOMMDev


def maketapershortCPW(width, straightlength, spacing, swidegndw, padlength,exptaperflag = False, textlabel = True):
    # same as other version but with double grounds
    # makes a meander using width. assumes that grounds = widths, calculates gaps
    # makes a SHORT. Takes probe spacing, wide single-side gnd width, pad lengths
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)

# make a straight section of CPW
    straightCPW = makeCPWseg_fixg(width,width,straightlength)

# make meander
    themeander = pd.Device('themeander')
    themeander.add_ref(straightCPW)
  #  pd.quickplot2(themeander)
# make full shorted cpw 
    short = pd.Device('shortCPW')
    # get pads based on probe spacing
    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    # use these to make a wide CPW
    gndw = swidegndw # should make gnd metal width widegndw
    sgndwidth =width # assuming gnds are CC width
    wideCPWseg = makeCPWseg_fixg(pad_w,gndw,padlength)
    # make taper
    widesegheight = wideCPWseg.ymax
    taperCPW = maketaperCPW_fixg(width,sgndwidth,spacing, swidegndw, exptaper = exptaperflag)
# assemble. 
# make bottom connector
    bottConnect = pd.Device('bottConn')
    bottpadref = bottConnect.add_ref(wideCPWseg)
    botttaperref = bottConnect.add_ref(taperCPW).move([0,widesegheight])
# now add it to the through
    bottConnref = short.add_ref(bottConnect).movex(destination=0)
    meanderCPWref = short.add_ref(themeander).move([0,bottConnect.ymax])
#make a shorting rect, move it into place
    topConnref = short.add_ref(pd.geometry.rectangle(size=(2*gap+3*width,3*width), layer = gold_layer)).move((-(gap+3*width/2),straightlength+bottConnect.ymax))
    textlabelref = pd.geometry.text("{} nm".format(round(width*1000, 0)), size = 20, layer = gold_layer, justify = 'center').move([0, -50])
    if textlabel is True:
        short.add_ref(textlabelref)
    #pd.quickplot2(short)
    return short


def makeHMOMMDot_v1(width, straightlength, spacing, swidegndw, padlength,dot_diam,exptaperflag = False, textlabel = True):
    
    HMOMMDev = pd.Device('HMOMMDev')
    
    # make the CPW structure #
    CPWtaper = maketapershortCPW(width, straightlength, spacing, swidegndw, padlength,exptaperflag , textlabel)
    
    # get dimensions to place crosses/dots
    narrowgap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
    # crosslen = max(narrowgap/2,3)
    # crosswid = max(crosslen/5,1)
    
    crosslen = 5
    crosswid = 2    

    
    HMOMMDev.add_ref(CPWtaper)

    fidmark = makecross(crosswid,crosslen)
    fidmarkvert = vertrec(crosswid,crosslen)
    
    for i in range(len(dot_diam)):
        dot_diam_wk = dot_diam[i]
        dot = pg.circle(dot_diam_wk/2,layer = mag_layer)
        fidmark_dot_line = pd.Device('fidmarkdot')
        #fidmark_arr.add_ref(dot)
        fidmark_dot_line.add_ref(fidmark)
        fidmark_dot_line.add_ref(dot).move((0,-(3-fidmark_dot_line.ymin+dot.ymax)))
        fidmark_dot_line.add_ref(fidmark).move((0,-(3-fidmark_dot_line.ymin+fidmark.ymax)))
        fidmark_dot_line.add_ref(dot).move((0,-(3-fidmark_dot_line.ymin+dot.ymax)))
        fidmark_dot_line.add_ref(fidmark).move((0,-(3-fidmark_dot_line.ymin+fidmark.ymax)))
        fidmark_dot_line.add_ref(dot).move((0,-(3-fidmark_dot_line.ymin+dot.ymax)))
        fidmark_dot_line.add_ref(fidmark).move((0,-(3-fidmark_dot_line.ymin+fidmark.ymax)))
        
        
        fidmark_arr = pd.Device('fidmarkarr')
        fidmark_arr.add_array(fidmark_dot_line,rows = 1, columns = 4, spacing =(width+narrowgap,0))
        
        
        
        center_arr_x = (fidmark_arr.xmax + fidmark_arr.xmin) / 2
        center_arr_y = (fidmark_arr.ymax + fidmark_arr.ymin) / 2
        fidmark_arr.move((-center_arr_x,-center_arr_y))
       
         
        textlabel = f"{np.round(dot_diam_wk,decimals = 3)}um"
        print(textlabel)                    
        textlabelref = pd.geometry.text(textlabel, size = 5, justify = 'center',layer = gold_layer)
        textlabelref.y = 0
        
        fidmark_arr.add_ref(textlabelref).move([-10-2.5*width-textlabelref.xmax,0])
    
        HMOMMDev.add_ref(fidmark_arr).move((0,CPWtaper.ymax-50-fidmark_arr.ymax-3*width-60*i))

    return HMOMMDev


def cutmarks(width, length,diesize):
    
    cross1 = pd.Device('cross1')
    cross2 = pd.Device('cross2')
    
    cross1.add_polygon([(0,0,length,length),(-width/2,width/2,width/2,-width/2)],layer = gold_layer).move([0,width/2])
    cross2.add_polygon([(-width/2,width/2,width/2,-width/2),(0,0,length,length)],layer = gold_layer).move([width/2,0])
    
    
    cross  = pg.boolean(A = cross1, B = cross2, operation = 'or',  precision = 1e-6,
               num_divisions = [1,1], layer = gold_layer)
    
    marks = pd.Device('marks')
    marks.add_ref(cross).rotate(270).move([-diesize/2,diesize/2])
    marks.add_ref(cross).rotate(180).move([diesize/2,diesize/2])
    marks.add_ref(cross).rotate(90).move([diesize/2,-diesize/2])
    marks.add_ref(cross).move([-diesize/2,-diesize/2])
    
    return marks
    


































