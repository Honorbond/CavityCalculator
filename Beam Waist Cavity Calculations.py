# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 22:28:54 2021

Transfer Matrices Calculator

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
thermal_lens_focal_length = 438.2*(pump_power**-1.229) #Thermal Lens Curve fot from paper



#Define or test distances here for beam waist calc
#d1 = sp.Symbol('d1')
#d2 = sp.Symbol('d2')

yag_index = 1.83
rod_length = .08 * yag_index
rod_height = 3

d1 = .267 + rod_length/2
d2 = .0205 + rod_length/2

#Thermal lens taken from paper in YAG Thermal Lens Folder (A Compact Side)





 #What range to search for values in


#Negative lens distance from rod face

neg_focal = inf

#Starting beam Parameters


#wavelength
y = 1064e-9





def gen_matrix(a,b,c,d):
    return(sp.Matrix([[a,b],[c,d]]))


#Ray Transfer Matrices

#Generic distance Matrices



def slider_wrapper(d1,d2,thermal_lens_focal_length,neg_focal,y):
    cavity_length = 2*d1+2*d2
    r = .0015
    R1 = inf
    R2 = inf
    theta = -arctan(r/thermal_lens_focal_length)
    Beam = sp.Matrix([[r],[theta]])
    
    d1NL = 0
    d2NL = d2 - d1NL

    
    
    
    S2M1 = gen_matrix(1, d1, 0, 1)
    M1 = gen_matrix(1, 0, -2/R1, 1)
    M12S = gen_matrix(1, d1, 0 , 1)
    S2M2 = gen_matrix(1, d2, 0 , 1)
    M2 = gen_matrix(1, 0, -2/R2, 1)
    M22S = gen_matrix(1, d2, 0, 1)
    TL = gen_matrix(1,0,-1/thermal_lens_focal_length,1)
    NL = gen_matrix(1,0,-1/neg_focal,1)
    TL2NL = gen_matrix(1,d1NL,0,1)
    NL2M2 = gen_matrix(1,d2NL,0,1)
    M22NL = gen_matrix(1,d2NL,0,1)
    NL2TL = gen_matrix(1,d1NL,0,1)
    
    
    
    
    
    
    
    dOutput = TL*NL2TL*NL*M22NL*M2*NL2M2*NL*TL2NL*TL*M12S*M1*S2M1
    
    dA = dOutput[0]
    dB = dOutput[1]
    dC = dOutput[2]
    dD = dOutput[3]
    
    
    #Multiply one by one and check the radius and divergence at each point
    OutputList = [TL,NL2TL,NL,M22NL,M2,NL2M2,NL,TL2NL,TL,M12S,M1,S2M1]
    beam_radius_pos = [Beam[0]]
    beam_radius_neg = [-1*Beam[0]]
    beam_divergence = [Beam[1]]
    
    for i in range(len(OutputList)):
        Beam = OutputList[len(OutputList)-1-i]*Beam
        beam_radius_pos.append(Beam[0])
        beam_radius_neg.append(-1*Beam[0])
        beam_divergence.append(Beam[1])
        
    
    d3 = d1NL
    d4 = d2NL
    
    bposmm = []
    bnegmm = []
    for i in range(len(beam_radius_pos)):
        bposmm.append(beam_radius_pos[i]*1000)
        bnegmm.append(beam_radius_neg[i]*1000)
        
    RL = rod_length/2
    
    distance_traveled = [0,d1,d1,2*d1,2*d1,2*d1+d3,2*d1+d3,2*d1+d3+d4,2*d1+d3+d4,2*d1+d3+2*d4,2*d1+d3+2*d4,2*d1+2*d3+2*d4,2*d1+2*d3+2*d4]
    return(distance_traveled,bposmm,bnegmm,dA,dD,d4,dB,R1,R2, cavity_length)



fig,axs = plt.subplots(1,2)
#plt.plot(distance_traveled,bnegmm, color ='r')
#plt.plot(distance_traveled,bposmm, color='r')


mirror_width = .001
mirror1_radius = 25.4/2 
mirror2_radius = 25.4/4








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






axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
d1_slider = Slider(
    ax=axfreq,
    label='D1',
    valmin=0,
    valmax=3,
    valinit=d1,
)

# Make a vertically oriented slider to control the amplitude
axamp = plt.axes([0.25, 0.05, 0.65, 0.03], facecolor=axcolor)
d2_slider = Slider(
    ax=axamp,
    label="D2",
    valmin=0,
    valmax=3,
    valinit=d2,
)


axtherm = plt.axes([0.1, 0.25, 0.0225, 0.63], facecolor=axcolor)
therm_slider = Slider(
    ax=axtherm,
    label="Thermal Lens \n Focal Length",
    valmin=0,
    valmax=5,
    valinit=thermal_lens_focal_length,
    orientation="vertical"
)


  





r1ax = plt.axes([0.05, 0.1, 0.15, 0.05])
r2ax = plt.axes([0.05, 0.05, 0.15, 0.05])


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
    
    
    k = slider_wrapper(d1_slider.val, d2_slider.val, therm_slider.val,neg_focal, y)
    line1, = axs[0].plot(k[0],k[1], color ='c')
    line1, = axs[0].plot(k[0],k[2], color ='c')
    axs[0].set_xlabel('Distance (mm)')
    
    
    
    rect1 = MakeRectangle(d1_slider.val,-mirror1_radius/2,mirror_width,mirror1_radius,'blue',False)
    AddPatch(rect1)

    rect2 = MakeRectangle(2*d1_slider.val-rod_length/2,-rod_height/2,rod_length,rod_height,'blue',False)
    AddPatch(rect2)

    rect3 = MakeRectangle(2*d1_slider.val+d2_slider.val,-mirror2_radius/2,mirror_width,mirror2_radius,'blue',False)
    AddPatch(rect3)

    rect4 = MakeRectangle(2*d1_slider.val+2*d2_slider.val-rod_length/2,-rod_height/2,rod_length,rod_height,'blue',False)
    AddPatch(rect4)

    rect5 = MakeRectangle(-rod_length/2,-rod_height/2,rod_length,rod_height,'blue',False)
    AddPatch(rect5)
    
    
    
    Stability = (k[3]+k[4]+2)/4
    print('Stability = ',Stability)


    
    
    
    if Stability < 0 or Stability > 1:
        print('Cavity is Unstable')
        resetax.set_facecolor('red')
        
    else:
        print('Stable')
        resetax.set_facecolor('green')
        
    beam_waist = (abs(sp.sqrt((y*abs(k[6]))/(pi*sp.sqrt(1-((k[3]+k[4])/2)**2)))*10**6))
    print('Beam waist = ', beam_waist,' microns')
    
    r1ax.clear()
    r2ax.clear()
    r1str = 'R1 = '+str(k[7])+' m'
    r2str = 'R2 = '+str(k[8])+' m'
    r1ax.text(.075,.25,r1str)
    r2ax.text(.075,.25,r2str)    
    
    fig.canvas.draw_idle()
    
    
    
    

def update2(val):
    
    plt.sca(axs[1])
    axs[1].clear()
    k = slider_wrapper(d1_slider.val, d2_slider.val, therm_slider.val,neg_focal, y)
    plt.title('Cavity Stability')

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
# Create a `matplotlib.widgets.Button` to reset the sliders to initial values.













"""
A=[]
B=[]
C=[]
D=[]
waist = []
output = []
length = linspace(0,cavity_length,1000)
for i in range(len(length)):
    
    output.append(Output.subs(d,length[i]))
    A.append(output[i][0])
    B.append(output[i][1])
    C.append(output[i][2])
    D.append(output[i][3])
    waist.append(abs(sp.sqrt((y*abs(B[i]))/(pi*sp.sqrt(1-((A[i]+D[i])/2)**2)))*10**6))

plt.figure()
plt.plot(100*length,waist)
plt.title('Length vs. Beam Waist')
plt.ylabel('Beam Waist (um)')
plt.xlabel('Arm length (Sample to Either End Mirror) (cm)')

"""
"""
for i in range(len(d)):
    #Sample to Mirror 1
    S2M1 = gen_matrix(1, d[i], 0, 1)

    #Mirror 1

    M1 = gen_matrix(1, 0, -2/R1, 1)

    #Mirror 1 to Sample
    M12S = gen_matrix(1, d[i], 0 , 1)

    #Sample to Mirror 2
    S2M2 = gen_matrix(1, d[i], 0 , 1)

    #Mirror 2
    M2 = gen_matrix(1, 0, -2/R2, 1)

    #Mirror 2 to Sample
    M22S = gen_matrix(1, d[i], 0, 1)

    #Multiply Matrices

    Output.append(M22S*M2*S2M2*M12S*M1*S2M1*Beam)
""" 


