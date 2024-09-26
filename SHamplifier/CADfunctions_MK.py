import phidl as pd
import phidl.routing as pr
import numpy as np
from scipy.special import ellipk
import CPW_Calc2
#import phidl.utilities
import phidl.geometry as pg
import phidl.path as pp

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


gold_layer = pd.Layer(gds_layer = 1, name='goldlayer', color='yellow')  # Contact pad layer
mag_layer = pd.Layer(gds_layer = 2, name = 'filmlayer', color = 'purple')
via_layer = pd.Layer(gds_layer = 3, name = 'etchlayer', color = 'green')
wafer_layer = pd.Layer(gds_layer = 4, name = 'waferlayer', color = 'red')
mark_layer = pd.Layer(gds_layer = 5, name = 'marklayer', color = 'blue')

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

def CPW_pair(width, straightlength, extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation):
    mDevice = pd.Device('cpw')
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
    extralength = ((2.5*gap)+(3*width)) 
    
    CPW1 = makemeandershortCPW2(width,straightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False, textlabel = False)
    CPW2 = makemeandershortCPW2(width,straightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False, textlabel = False).rotate(180, center = [((CPW1.xmax+CPW1.xmin)/2),CPW1.ymin]).move([separation,(2*(padlength+taper_length))+straightlength+(2*extralength)])
    
    mDevice.add_ref(CPW1)
    mDevice.add_ref(CPW2)

    return(mDevice)


def CPW_single(width, straightlength, extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength, taper_length, separation):
    mDevice = pd.Device('cpw')
    gap = CPW_Calc2.getgap(width, thickness, eps_eff, impedance, sub_eps)
    extralength = ((2.5*gap)+(3*width)) 
    
    CPW1 = makemeandershortCPW2(width,straightlength+extendpastswchannel, numperiods, GSGspacing, swidegndw, padlength,exptaperflag = False, textlabel = False)
    
    mDevice.add_ref(CPW1)

    return(mDevice)


def DCpath(width, straightlength, via_extend_x, via_extend_y, dcleadwidth, dctaperlength, numCh, ydistCh, padlength, taper_length, overlaptrim, separation, extend1, extend2, ProbeSpacing, constlen, constrwdth):

    mDevice = pd.Device('DCpath')
    SWchan = pd.Device('SWchannel')

    viashift = (separation+extend1+extend2+2*via_extend_x)

    len_box1 = separation/2+extend1+via_extend_x-1.5*constlen
    len_box2 = separation/2+extend2+via_extend_x-1.5*constlen

    
    box1 = SWchan << pg.compass([len_box1,straightlength], layer = mag_layer).move([-viashift/4-constlen*0.75,0])
    box2 = SWchan << pg.compass([len_box2,straightlength], layer = mag_layer).move([+viashift/4+constlen*0.75,0])
    box3 = SWchan << pg.compass([constlen,constrwdth], layer = mag_layer).move([0,0])


    constriction = pr.route_smooth(port1 = box1.ports['E'], port2 = box3.ports['W'], smooth_options=  {'corner_fun': pp.euler, 'use_eff': True}, layer = mag_layer,)

    SWchan.add_ref(constriction)
    SWchan.add_ref(constriction).rotate(angle = 180, center = (0,0))

    Via1 = pd.Device('via1')
    Via2 = pd.Device('via2')
    
    xvia = [-via_extend_x, via_extend_x, via_extend_x, -via_extend_x]
    yvia = [-via_extend_y, -via_extend_y, straightlength+via_extend_y, straightlength+via_extend_y]

    # Make via polygons
    Via1.add_polygon([xvia, yvia], layer = via_layer)
    Via2.add_polygon([xvia, yvia], layer = via_layer)


    # place vias
    Via1.move([-viashift/2+via_extend_x, -straightlength/2])
    Via2.move([viashift/2-via_extend_x, -straightlength/2])
    #Via1.move([-viashift/2, -straightlength/2])
    #Via2.move([viashift/2, -straightlength/2])


    Lead1 = pd.Device('lead1')
    Lead2 = pd.Device('lead2')

    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    dcleadlength = (pad_w*3+pad_gap*2)/2-separation-dctaperlength+20 

    xlead = [-(dcleadlength/2), -(dcleadlength/2), (dcleadlength/2), (dcleadlength/2)]
    ylead = [-(dcleadwidth/2), (dcleadwidth/2), (dcleadwidth/2), -(dcleadwidth/2)]

    xpad = [-(padlength/2), -(padlength/2), (padlength/2), (padlength/2)]
    ypad = [-(pad_w/2), (pad_w/2), (pad_w/2), -(pad_w/2)] 


    for i in range(numCh):
        mDevice.add_ref(SWchan).move([0, -ydistCh*i])
        mDevice.add_ref(Via1).move([0, -ydistCh*i])
        mDevice.add_ref(Via2).move([0, -ydistCh*i])


    centerx = ((mDevice.xmin + mDevice.xmax)/2) 
    centery = ((mDevice.ymin + mDevice.ymax)/2)
    mDevice.move([centerx, -centery])


    yleadtaper = []

    for i in range(numCh):
        Lead1.add_polygon([xlead, ylead], layer = gold_layer).move([-dctaperlength-viashift/2+2*via_extend_x, dcleadwidth*(2*i-numCh+1)+(numCh-1)*(ydistCh)/2])
        Lead1.add_polygon([xpad, ypad], layer = gold_layer).move([-ProbeSpacing, taper_length+numCh*ydistCh+pad_w/2+(pad_w+pad_gap)*i])

        Lead1.add_polygon([xvia, yvia], layer = gold_layer).move([-viashift/2+via_extend_x, -straightlength/2+ydistCh*i])
        xleadtaper = [Lead1.xmax-2*via_extend_x, Lead1.xmax-2*via_extend_x, Lead1.xmax-dctaperlength+dcleadlength/2, Lead1.xmax-dctaperlength+dcleadlength/2] 
        yleadtaper.append([Via1.ymax+ydistCh*i, Via1.ymin+ydistCh*i, Lead1.ymin+2*dcleadwidth*i, Lead1.ymin+2*dcleadwidth*i+dcleadwidth])


    for i in range(numCh):
        Lead1.add_polygon([xleadtaper, yleadtaper[i]], layer = gold_layer)

    Lead1.move([centerx, centery])

    for i in range(numCh):
        port1 = Lead1.add_port(name='U'+str(i*2+1), midpoint=(Lead1.xmax-dctaperlength-dcleadlength/2, 2*dcleadwidth*i-(numCh-1)*dcleadwidth), width=dcleadwidth, orientation=180)
        port2 = Lead1.add_port(name='U'+str(i*2+2), midpoint=(Lead1.xmin+padlength, taper_length+(numCh-1)*ydistCh+pad_w/2+(pad_w+pad_gap)*i), width=dcleadwidth, orientation=0)
        connect1 = pr.route_smooth(port1, port2, path_type = 'U', length1 = 200-20*i, layer = gold_layer)
        Lead1.add_ref(connect1).move([0, 0])

    Lead2.add_ref(Lead1).rotate(angle = 180, center = (0,0)).move([0, 0])
    mDevice.add_ref(Lead1)
    mDevice.add_ref(Lead2)
    mDevice.mirror((0, 0))

    return(mDevice)

def DCpath_single(width, straightlength, via_extend_x, via_extend_y, dcleadwidth, dctaperlength, numCh, ydistCh, padlength, taper_length, add_cdt_ext, separation, extend1, extend2, ProbeSpacing, constlen, constrwdth,constrtpr):

    mDevice = pd.Device('DCpath')
    SWchan = pd.Device('SWchannel')

    viashift = extend1+extend2+constlen+2*constrtpr+2*via_extend_x
    separation = constlen+2*constrtpr


    box1 = SWchan << pg.compass([extend1,straightlength], layer = mag_layer).move([-constlen/2-constrtpr-extend1/2,0])
    box2 = SWchan << pg.compass([extend2,straightlength], layer = mag_layer).move([extend2/2+constlen/2+constrtpr,0])
    box3 = SWchan << pg.compass([constlen,constrwdth], layer = mag_layer).move([0,0])

    constriction = pr.route_smooth(port1 = box1.ports['E'], port2 = box3.ports['W'], smooth_options=  {'corner_fun': pp.euler, 'use_eff': True}, layer = mag_layer,)

    SWchan.add_ref(constriction)
    SWchan.add_ref(constriction).rotate(angle = 180, center = (0,0))
    #SWchan.move([0,0])

    # hardcoded extention past the via to avoid SW reflection
    ext_len = 5
    if add_cdt_ext == 1:
        #SWchan.add_polygon([[-ext_len/2,ext_len/2,ext_len/2,-ext_len/2-straightlength], [0,0,straightlength,straightlength]], layer = mag_layer).move([-extend1-constrtpr-constlen/2-ext_len/2,-straightlength/2])
        #SWchan.add_polygon([[-ext_len/2,ext_len/2-straightlength,ext_len/2,-ext_len/2], [0,0,straightlength,straightlength]], layer = mag_layer).move([extend2+constrtpr+constlen/2+ext_len/2,-straightlength/2])
        SWchan.add_polygon([[-ext_len/2,ext_len/2,ext_len/2,-ext_len/2-straightlength], [0,0,straightlength,straightlength]], layer = mag_layer).move([-extend1-constrtpr-constlen/2-ext_len/2,-straightlength/2])
        SWchan.add_polygon([[-ext_len/2,ext_len/2,ext_len/2+straightlength,-ext_len/2], [0,0,straightlength,straightlength]], layer = mag_layer).move([extend2+constrtpr+constlen/2+ext_len/2,-straightlength/2])

    #box4 = SWchan << pg.compass([5,straightlength], layer = mag_layer).move([-extend1-constrtpr-constlen/2-2.5,0])
    #box5 = SWchan << pg.compass([5,straightlength], layer = mag_layer).move([extend2+constrtpr+constlen/2+2.5,0])


    Via1 = pd.Device('via1')
    Via2 = pd.Device('via2')
    
    xvia = [-via_extend_x, via_extend_x, via_extend_x, -via_extend_x]
    yvia = [-via_extend_y, -via_extend_y, straightlength+via_extend_y, straightlength+via_extend_y]

    # Make via polygons
    Via1.add_polygon([xvia, yvia], layer = via_layer)
    Via2.add_polygon([xvia, yvia], layer = via_layer)

    # place vias
    Via1.move([-extend1-constrtpr-constlen/2+via_extend_x, -straightlength/2])
    Via2.move([extend2+constrtpr+constlen/2-via_extend_x, -straightlength/2])


    Lead1 = pd.Device('lead1')
    Lead2 = pd.Device('lead2')

    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    dcleadlength = (pad_w*3+pad_gap*2)/2-separation-dctaperlength+20 

    xlead = [-(dcleadlength/2), -(dcleadlength/2), (dcleadlength/2), (dcleadlength/2)]
    ylead = [-(dcleadwidth/2), (dcleadwidth/2), (dcleadwidth/2), -(dcleadwidth/2)]

    xpad = [-(padlength/2), -(padlength/2), (padlength/2), (padlength/2)]
    ypad = [-(pad_w/2), (pad_w/2), (pad_w/2), -(pad_w/2)] 


    for i in range(numCh):
        mDevice.add_ref(SWchan).move([0, -ydistCh*i])
        mDevice.add_ref(Via1).move([0, -ydistCh*i])
        mDevice.add_ref(Via2).move([0, -ydistCh*i])


    centerx = ((mDevice.xmin + mDevice.xmax)/2) 
    centery = ((mDevice.ymin + mDevice.ymax)/2)


    mDevice.move([-centerx, -centery])


    yleadtaper = []

    for i in range(numCh):
        Lead1.add_polygon([xlead, ylead], layer = gold_layer).move([-centerx-dctaperlength-viashift/2+3*via_extend_x, dcleadwidth*(2*i-numCh+1)+(numCh-1)*(ydistCh)/2])
        Lead1.add_polygon([xpad, ypad], layer = gold_layer).move([-ProbeSpacing, taper_length+numCh*ydistCh+pad_w/2+(pad_w+pad_gap)*i])

        Lead1.add_polygon([xvia, yvia], layer = gold_layer).move([-centerx-viashift/2+2*via_extend_x, -straightlength/2+ydistCh*i])
        xleadtaper = [Lead1.xmax-2*via_extend_x, Lead1.xmax-2*via_extend_x, Lead1.xmax-dctaperlength+dcleadlength/2, Lead1.xmax-dctaperlength+dcleadlength/2] 
        yleadtaper.append([Via1.ymax+ydistCh*i, Via1.ymin+ydistCh*i, Lead1.ymin+2*dcleadwidth*i, Lead1.ymin+2*dcleadwidth*i+dcleadwidth])


    for i in range(numCh):
        Lead1.add_polygon([xleadtaper, yleadtaper[i]], layer = gold_layer)

    Lead1.move([centerx, centery])

    for i in range(numCh):
        port1 = Lead1.add_port(name='U'+str(i*2+1), midpoint=(Lead1.xmax-dctaperlength-dcleadlength/2, 2*dcleadwidth*i-(numCh-1)*dcleadwidth), width=dcleadwidth, orientation=180)
        port2 = Lead1.add_port(name='U'+str(i*2+2), midpoint=(Lead1.xmin+padlength, taper_length+(numCh-1)*ydistCh+pad_w/2+(pad_w+pad_gap)*i), width=dcleadwidth, orientation=0)
        connect1 = pr.route_smooth(port1, port2, path_type = 'U', length1 = 200-20*i, layer = gold_layer)
        Lead1.add_ref(connect1).move([0, 0])

    Lead2.add_ref(Lead1).rotate(angle = 180, center = (0,0)).move([0, 0])
    mDevice.add_ref(Lead1)
    mDevice.add_ref(Lead2)
    mDevice.mirror((0, 0))


    if constrwdth == straightlength:
        textlabel = f"{np.round(viashift-6*via_extend_x,1)} um x {straightlength} um"
    
    else:
        textlabel = f"{np.round(viashift-6*via_extend_x,1)} um x {straightlength} um ({constlen}um x {constrwdth}um)"

    textlabelref = pd.geometry.text(textlabel, size = 25, layer = gold_layer, justify = 'center').move([centerx, centery-1300])
    mDevice.add_ref(textlabelref)

    return(mDevice)


def DCtest_path(straightlength, via_extend_x, via_extend_y, dcleadwidth, dctaperlength, extend1, extend2, padlength, overlaptrim, separation, ProbeSpacing, constlen, constrwdth,constrtpr):

    mDevice = pd.Device('DCpath')
    SWchan = pd.Device('SWchannel')

    viashift = (extend1+extend2+constlen+2*constrtpr)
    separation = constlen+2*constrtpr


    box1 = SWchan << pg.compass([extend1,straightlength], layer = mag_layer).move([-constlen/2-constrtpr-extend1/2,0])
    box2 = SWchan << pg.compass([extend2,straightlength], layer = mag_layer).move([extend2/2+constlen/2+constrtpr,0])
    box3 = SWchan << pg.compass([constlen,constrwdth], layer = mag_layer).move([0,0])

    constriction = pr.route_smooth(port1 = box1.ports['E'], port2 = box3.ports['W'], smooth_options=  {'corner_fun': pp.euler, 'use_eff': True}, layer = mag_layer,)

    SWchan.add_ref(constriction)
    SWchan.add_ref(constriction).rotate(angle = 180, center = (0,0))


    Via1 = pd.Device('via1')
    Via2 = pd.Device('via2')
    
    xvia = [-via_extend_x, via_extend_x, via_extend_x, -via_extend_x]
    yvia = [-via_extend_y, -via_extend_y, straightlength+via_extend_y, straightlength+via_extend_y]

    # Make via polygons
    Via1.add_polygon([xvia, yvia], layer = via_layer)
    Via2.add_polygon([xvia, yvia], layer = via_layer)


    # place vias
    Via1.move([-extend1-constrtpr-constlen/2+via_extend_x, -straightlength/2])
    Via2.move([extend2+constrtpr+constlen/2-via_extend_x, -straightlength/2])

    Lead1 = pd.Device('lead1')
    Lead2 = pd.Device('lead2')

    pad_w, pad_gap = CPW_Calc2.getpads(spacing, thickness, eps_eff, impedance, sub_eps)
    dcleadlength = (pad_w*3+pad_gap*2)/2-separation-dctaperlength+20 

    xlead = [-(dcleadlength/2), -(dcleadlength/2), (dcleadlength/2), (dcleadlength/2)]
    ylead = [-(dcleadwidth/2), (dcleadwidth/2), (dcleadwidth/2), -(dcleadwidth/2)]

    xpad = [-(padlength/2), -(padlength/2), (padlength/2), (padlength/2)]
    ypad = [-(pad_w/2), (pad_w/2), (pad_w/2), -(pad_w/2)] 

    mDevice.add_ref(SWchan)
    mDevice.add_ref(Via1)
    mDevice.add_ref(Via2)

    centerx = ((mDevice.xmin + mDevice.xmax)/2) 
    centery = ((mDevice.ymin + mDevice.ymax)/2)
    mDevice.move([-centerx, -centery])

    


    Lead1.add_polygon([xlead, ylead], layer = gold_layer).move([-centerx-dctaperlength-viashift/2-dcleadlength/2, 0])
    Lead1.add_polygon([xpad, ypad], layer = gold_layer).move([-centerx-ProbeSpacing-viashift/2, -centery])
    Lead1.add_polygon([xvia, yvia], layer = gold_layer).move([viashift/2-via_extend_x, -straightlength/2])


    xleadtaper1 = [-viashift/2, -viashift/2, -centerx-viashift/2-dctaperlength, -centerx-viashift/2-dctaperlength] 
    yleadtaper1 = [Via1.ymax, Via1.ymin, Lead1.ymax-pad_w/2-dcleadwidth/2, Lead1.ymax-pad_w/2+dcleadwidth/2]

    xleadtaper2 = [Lead1.xmin+padlength,Lead1.xmin+padlength,-centerx-viashift/2-dctaperlength-dcleadlength,-centerx-viashift/2-dctaperlength-dcleadlength] 
    yleadtaper2 = [-pad_w/2,pad_w/2,dcleadwidth/2,-dcleadwidth/2]

    Lead1.add_polygon([xleadtaper1, yleadtaper1], layer = gold_layer)
    Lead1.add_polygon([xleadtaper2, yleadtaper2], layer = gold_layer)

    Lead2.add_ref(Lead1).rotate(angle = 180, center = (0,0)).move([0, 0])
    mDevice.add_ref(Lead1)
    mDevice.add_ref(Lead2)

    DeviceArray = pd.Device('Testarray')
    DeviceArray.add_ref(mDevice)

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
    
    DeviceArray.add_ref(CornerSquare)



    if constrwdth == straightlength:
        textlabel = f"{np.round(viashift-4*via_extend_x,1)} um x {straightlength} um"
    
    else:
        textlabel = f"{np.round(viashift-4*via_extend_x,1)} um x {straightlength} um ({constlen}um x {constrwdth}um)"

    #print(textlabel)                         
    textlabelref = pd.geometry.text(textlabel, size = 25, layer = gold_layer, justify = 'center').move([-centerx-ProbeSpacing+100, pad_w/2+28])
    mDevice.add_ref(textlabelref)

    #for j in range(numRow):
    #    for i in range(numCh):
    #        DeviceArray.add_ref(mDevice).move([j*2*(Lead1.xmin-pad_w/2), (2*pad_w+pad_gap)*i])


    return(DeviceArray)


def fid_squares(dev, corn_offset, corn_width):
        CornerSquare = pd.Device('cornerSquare')
    
        #corn_offset = 50  # offset of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
        #corn_width = 10   # width of fiducial squares in corners to elimiate bulk/sleeve overlap issues with BEAMER
        
        corn_xmin = dev.xmin 
        corn_xmax = dev.xmax
        corn_ymin = dev.ymin
        corn_ymax = dev.ymax
        
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
        
        dev.add_ref(CornerSquare)
    



def EBL_marks():
    
    EBL_mark = pd.Device('EBLMark')
    cross1 = pd.Device('cross1')
    cross2 = pd.Device('cross2')
    
    cross1.add_polygon([(-1000,-1000,1000,1000),(-1.5,1.5,1.5,-1.5)],layer = mark_layer)
    cross2.add_polygon([(-1.5,1.5,1.5,-1.5),(-1000,-1000,1000,1000)],layer = mark_layer)
    
    cross  = pg.boolean(A = cross1, B = cross2, operation = 'or',  precision = 1e-6,
               num_divisions = [1,1], layer = mark_layer)
    
    
    EBL_mark.add_ref(cross).move((-32500,2000))
    EBL_mark.add_ref(cross).move((32500,2000))
    
    
    return EBL_mark

def CenterDev(dev):
    
    dev.x = 0
    dev.y = 0
    
    return dev

def cutmarks(width, length,diesize):
    
    cross1 = pd.Device('cross1')
    cross2 = pd.Device('cross2')
    
    cross1.add_polygon([(0,0,length,length),(-width/2,width/2,width/2,-width/2)],layer = mark_layer).move([0,width/2])
    cross2.add_polygon([(-width/2,width/2,width/2,-width/2),(0,0,length,length)],layer = mark_layer).move([width/2,0])
    
    
    cross  = pg.boolean(A = cross1, B = cross2, operation = 'or',  precision = 1e-6,
               num_divisions = [1,1], layer = mark_layer)
    
    marks = pd.Device('marks')
    marks.add_ref(cross).rotate(270).move([-diesize/2,diesize/2])
    marks.add_ref(cross).rotate(180).move([diesize/2,diesize/2])
    marks.add_ref(cross).rotate(90).move([diesize/2,-diesize/2])
    marks.add_ref(cross).move([-diesize/2,-diesize/2])
    
    return marks


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
    xpad = [-(pad_w), (pad_w), (pad_w), -(pad_w)] #doubled the width per Matt's request
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

