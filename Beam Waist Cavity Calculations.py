# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 12:30:18 2022

@author: rubin
"""
from numpy import *
import sympy as sp
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PlotUtilities import *
from matplotlib.widgets import Slider, Button


plt.close('all')

#Radii of Curvature in M


#Starting beam parameters
#r = sp.Symbol('r')
#theta = sp.Symbol('theta')

Pump_current = 14 #A 
pump_power = (.0757*Pump_current+19.014)*Pump_current  #Pump Power Curve Fit
thermal_lens_focal_length = .85#438.2*(pump_power**-1.229) #Thermal Lens Curve fot from paper


#Define or test distances here for beam waist calc
#d1 = sp.Symbol('d1')
#d2 = sp.Symbol('d2')

yag_index = 1.83
rod_length = .08 * yag_index
rod_height = 3

d1 = .267 + rod_length/2
d2 = .0205 + rod_length/2
d3 = sp.Symbol('d3')

#Thermal lens taken from paper in YAG Thermal Lens Folder (A Compact Side)

 #What range to search for values in


#Negative lens distance from rod face

neg_focal = inf

#Starting beam Parameters


#wavelength
y = 1064e-9

def gen_matrix_np (a,b,c,d):
    return(array([[a,b],[c,d]]))

def mm (mat1,mat2):
    return(matmul(mat1,mat2))

def slider_wrapper(d1,d2,d3,thermal_lens_focal_length,neg_focal,y):
    cavity_length = 2*d1+2*d2
    r = .0015
    R1 = inf
    R2 = .5
    R3 = .3
    theta = -arctan(r/thermal_lens_focal_length)
    Beam = array([[r],[theta]])
    
    
    d2NL = d2
    d1NL = 0
    
    
    
    S2M1 = gen_matrix_np(1, d1, 0, 1)
    M1 = gen_matrix_np(1, 0, -2/R1, 1)
    M12S = gen_matrix_np(1, d1, 0 , 1)
    S2M2 = gen_matrix_np(1, d2, 0 , 1)
    M2 = gen_matrix_np(1, 0, -2/R2, 1)
    M22S = gen_matrix_np(1, d2, 0, 1)
    TL = gen_matrix_np(1,0,-1/thermal_lens_focal_length,1)
    NL = gen_matrix_np(1,0,-1/neg_focal,1)
    TL2NL = gen_matrix_np(1,d1NL,0,1)
    NL2M2 = gen_matrix_np(1,d2NL,0,1)
    M22NL = gen_matrix_np(1,d2NL,0,1)
    NL2TL = gen_matrix_np(1,d1NL,0,1)
    
    #Mirror 3
    
    
    
    M3 = gen_matrix_np(1,0,-2/R3,1)
    M22M3 = gen_matrix_np(1,d3,0,1)
     
    m1 = mm(M1,S2M1)
    #m = mm(S2M1,Beam)
    m2 = mm(M1,m1)
    m3 = mm(M12S,m2)
    m4 = mm(TL,m3)
    m5 = mm(NL2M2,m4)
    m6 = mm(M2,m5)
    m7 = mm (M22M3,m6)
    beforeM3 = m7
    m8 = mm(M3,m7)
    M3mat = mm(M3,m7)
    afterM3 = m8
    m9 = mm(M22M3,m8)
    m10 = mm(M2,m9)
    m11 = mm(M22NL,m10)
    m12 = mm(NL,m11)
    m13 = mm(NL2TL,m12)
    m14 = mm(TL,m13)
    
    
    m = m14
    
    
    
    #inbeamM3 = beforeM3[1][0]
    #outbeamM3 = afterM3[1][0]
    
    
    
    
    inbeamM3 = mm(M22M3,mm(M2,mm(NL2M2,mm(TL,mm(M12S,mm(M1,mm(S2M1,Beam)))))))
    outbeamM3 = mm(M3,mm(M22M3,mm(M2,mm(NL2M2,mm(TL,mm(M12S,mm(M1,mm(S2M1,Beam))))))))
    inbeamM3 = inbeamM3[1][0]
    outbeamM3 = outbeamM3[1][0]
    d3val = sp.solve(outbeamM3 + inbeamM3,d3)
    inter = m
    inter = array([[inter[0][0].subs(d3,d3val[0]),inter[0][1].subs(d3,d3val[0])],[inter[1][0].subs(d3,d3val[0]),inter[1][1].subs(d3,d3val[0])]])
    m = inter
    d3 = d3val[0]
    
    S2M1 = gen_matrix_np(1, d1, 0, 1)
    M1 = gen_matrix_np(1, 0, -2/R1, 1)
    M12S = gen_matrix_np(1, d1, 0 , 1)
    S2M2 = gen_matrix_np(1, d2, 0 , 1)
    M2 = gen_matrix_np(1, 0, -2/R2, 1)
    M22S = gen_matrix_np(1, d2, 0, 1)
    TL = gen_matrix_np(1,0,-1/thermal_lens_focal_length,1)
    NL = gen_matrix_np(1,0,-1/neg_focal,1)
    TL2NL = gen_matrix_np(1,d1NL,0,1)
    NL2M2 = gen_matrix_np(1,d2NL,0,1)
    M22NL = gen_matrix_np(1,d2NL,0,1)
    NL2TL = gen_matrix_np(1,d1NL,0,1)
    
    #Mirror 3
    
    
    
    M3 = gen_matrix_np(1,0,-2/R3,1)
    M22M3 = gen_matrix_np(1,d3,0,1)
    
   
    
   
    
   
    
   
    #dOutput = TL*NL2TL*NL*M22NL*M2*M22M3*M3*M22M3*M2*NL2M2*TL*M12S*M1*S2M1
    
    
    
    
    dA = m[0][0]
    dB = m[0][1]
    dC = m[1][0]
    dD = m[1][1]
    
    
    
    
    
    
    #Multiply one by one and check the radius and divergence at each point
    OutputList = [TL,M22NL,M2,M22M3,M3,M22M3,M2,NL2M2,TL,M12S,M1,S2M1]
    beam_radius_pos = [Beam[0][0]]
    beam_radius_neg = [-1*Beam[0][0]]
    beam_divergence = [Beam[1][0]]
    
    for i in range(len(OutputList)):
        Beam = mm(OutputList[len(OutputList)-1-i],Beam)
        beam_radius_pos.append(Beam[0][0])
        beam_radius_neg.append(-1*Beam[0][0])
        beam_divergence.append(Beam[1][0])
        
    
    
    
    bposmm = []
    bnegmm = []
    for i in range(len(beam_radius_pos)):
        bposmm.append(beam_radius_pos[i]*1000)
        bnegmm.append(beam_radius_neg[i]*1000)
        
    RL = rod_length/2
    
    #distance_traveled = [0,d1,d1,2*d1,2*d1,2*d1+d3,2*d1+d3,2*d1+d3+d4,2*d1+d3+d4,2*d1+d3+2*d4,2*d1+d3+2*d4,2*d1+2*d3+2*d4,2*d1+2*d3+2*d4]
    distance_traveled = [0,d1,d1,2*d1,2*d1,2*d1+d2,2*d1+d2,2*d1+d2+d3,2*d1+d2+d3,2*d1+d2+2*d3,2*d1+d2+2*d3,2*d1+2*d2+2*d3,2*d1+2*d2+2*d3]
    return(distance_traveled,bposmm,bnegmm,dA,dD,d2,dB,R1,R2, cavity_length,R3,d3val[0])



fig,axs = plt.subplots(1,2)
#plt.plot(distance_traveled,bnegmm, color ='r')
#plt.plot(distance_traveled,bposmm, color='r')


mirror_width = .001
mirror1_radius = 25.4/2 
mirror2_radius = 25.4/4
mirror3_radius = 25.4/2


d3init = 0.5

plt.ylabel('Beam Radius(mm)')
plt.xlabel('Cavity Length (m)')



#line, = plt.plot(t, f(t, init_amplitude, init_frequency), lw=2)
#line1, = 
#line2, = 
axs[0].set_xlabel('Distance (mm)')

axcolor = 'lightgoldenrodyellow'
axs[0].margins(x=0)

# adjust the main plot to make room for the sliders
plt.subplots_adjust(left=0.25, bottom=0.25)



axfreq = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
d1_slider = Slider(
    ax=axfreq,
    label='D1',
    valmin=0,
    valmax=3,
    valinit=d1,
)

# Make a vertically oriented slider to control the amplitude
axamp = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
d2_slider = Slider(
    ax=axamp,
    label="D2",
    valmin=0,
    valmax=3,
    valinit=d2,
)

"""
axd3 = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=axcolor)
d3_slider = Slider(
    ax=axd3,
    label="D3",
    valmin=0,
    valmax=3,
    valinit=d3init,
)
"""


axd3 = plt.axes([0.25, 0.05, 0.65, 0.03])

axtherm = plt.axes([0.1, 0.25, 0.0225, 0.63], facecolor=axcolor)
therm_slider = Slider(
    ax=axtherm,
    label="Thermal Lens \n Focal Length",
    valmin=0,
    valmax=5,
    valinit=thermal_lens_focal_length,
    orientation="vertical"
)


  
r1ax = plt.axes([0.05, 0.15, 0.15, 0.05])
r2ax = plt.axes([0.05, 0.1, 0.15, 0.05])
r3ax = plt.axes([0.05,0.05,0.15,0.05])
r4ax = plt.axes([0.05,0,0.15,0.05])


resetax = plt.axes([0.05, 0.15, 0.15, 0.05])
resetax.tick_params(axis = 'x',
                    which = 'both',
                    bottom=False,
                    top=False,
                    labelbottom=False)
resetax.tick_params(axis = 'y',
                    which = 'both',
                    left=False,
                    right=False,
                    labelleft=False)
resetax.text(.075,.25,'Stability')



def update(val):
    
    
    plt.sca(axs[0])
    axs[0].clear()
    plt.title('Cavity Length vs. Beam Radius')
    
    
    k = slider_wrapper(d1_slider.val, d2_slider.val,d3, therm_slider.val,neg_focal, y)
    line1, = axs[0].plot(k[0],k[1], color ='c')
    line1, = axs[0].plot(k[0],k[2], color ='c')
    axs[0].set_xlabel('Distance (mm)')
    d3val = k[11]
    
    
    
    rect1 = MakeRectangle(d1_slider.val,-mirror1_radius/2,mirror_width,mirror1_radius,'blue',False)
    AddPatch(rect1)

    rect2 = MakeRectangle(2*d1_slider.val-rod_length/2,-rod_height/2,rod_length,rod_height,'blue',False)
    AddPatch(rect2)

    rect3 = MakeRectangle(2*d1_slider.val+d2_slider.val,-mirror2_radius/2,mirror_width,mirror2_radius,'blue',False)
    AddPatch(rect3)

    rect4 = MakeRectangle(2*d1_slider.val+2*d2_slider.val+2*d3val-rod_length/2,-rod_height/2,rod_length,rod_height,'blue',False)
    AddPatch(rect4)

    rect5 = MakeRectangle(-rod_length/2,-rod_height/2,rod_length,rod_height,'blue',False)
    AddPatch(rect5)
    
    rect3 = MakeRectangle(2*d1_slider.val+d2_slider.val+d3val,-mirror3_radius/2,mirror_width,mirror3_radius,'blue',False)
    AddPatch(rect3)
    
    rect3 = MakeRectangle(2*d1_slider.val+d2_slider.val+2*d3val,-mirror2_radius/2,mirror_width,mirror2_radius,'blue',False)
    AddPatch(rect3)
    
    
    Stability = (k[3]+k[4]+2)/4
    print('Stability = ',Stability)


    
    if Stability < 0 or Stability > 1:
        print('Cavity is Unstable')
        resetax.set_facecolor('red')
        
    else:
        print('Stable')
        resetax.set_facecolor('green')
        
    beam_waist = (abs(sp.sqrt((y*abs(k[6]))/(pi*sp.sqrt(1-((k[3]+k[4])/2)**2)))*10**6))
    #print('Beam waist = ', beam_waist,' microns')
    
    r1ax.clear()
    r2ax.clear()
    r3ax.clear()
    r4ax.clear()
    axd3.clear()
    
    axd3.text(.075,.25,'d3 = '+str(d3val))
    r4str = 'Beam Waist = ', int(beam_waist),' microns'
    r1str = 'R1 = '+str(k[7])+' m'
    r2str = 'R2 = '+str(k[8])+' m'
    r3str = 'R3 = ' +str(k[10])+' m'
    r1ax.text(.075,.25,r1str)
    r2ax.text(.075,.25,r2str) 
    r3ax.text(.075,.25,r3str)
    r4ax.text(.075,.25,r4str)
    
    
    
    r1ax.tick_params(axis = 'x',
                    which = 'both',
                    bottom=False,
                    top=False,
                    labelbottom=False)
    r1ax.tick_params(axis = 'y',
                        which = 'both',
                        left=False,
                        right=False,
                        labelleft=False)
    r2ax.tick_params(axis = 'x',
                        which = 'both',
                        bottom=False,
                        top=False,
                        labelbottom=False)
    r2ax.tick_params(axis = 'y',
                        which = 'both',
                        left=False,
                        right=False,
                        labelleft=False)
    r3ax.tick_params(axis = 'x',
                        which = 'both',
                        bottom=False,
                        top=False,
                        labelbottom=False)
    r3ax.tick_params(axis = 'y',
                        which = 'both',
                        left=False,
                        right=False,
                        labelleft=False)
    r4ax.tick_params(axis = 'x',
                        which = 'both',
                        bottom=False,
                        top=False,
                        labelbottom=False)
    r4ax.tick_params(axis = 'y',
                        which = 'both',
                        left=False,
                        right=False,
                        labelleft=False)
    
    
    
    
    fig.canvas.draw_idle()
    
    

def update2(val):
    
    plt.sca(axs[1])
    axs[1].clear()
    k = slider_wrapper(d1_slider.val, d2_slider.val,d3,therm_slider.val,neg_focal, y)
    plt.title('Cavity Stability')
    d3val = k[11]
    A = k[3]
    D = k[4]
    stab = (A+D+2)/4
    x = linspace(-3,3,1000)
    y1=[]
    y2=[]
    y3=[]
    y4=[]
    
    
    for i in range(len(x)):
            
        y1.append(1)
        y2.append(3)
        y3.append(0)
        y4.append(-1)
    
    plt.axhline(y=1)
    plt.axhline(y=0)
    plt.plot(0,stab, marker='o',markersize=10,color='b')
    plt.xlim(-1,1)
    plt.ylim(-1,2)
    plt.fill_between(x,y1,y2,color='r')
    plt.fill_between(x,y3,y4,color='r')
    plt.fill_between(x,y1,y3,color='g')
    
    print(k[9])
    
    """
    print('g1 = ',g1)
    print('g2 = ',g2)
    print('g1g2 = ', g1*g2)
    """
    
    
    fig.canvas.draw_idle()
    
    
# register the update function with each slider
d1_slider.on_changed(update)
d2_slider.on_changed(update)
d1_slider.on_changed(update2)
d2_slider.on_changed(update2)
therm_slider.on_changed(update)
therm_slider.on_changed(update2)
#d3_slider.on_changed(update)
#d3_slider.on_changed(update2)
