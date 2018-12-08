# -*- coding: utf-8 -*-
"""
Created on Fri Dec 07 11:53:26 2018

@author: Erfan
"""

import CoolProp.CoolProp as cp
import numpy as np
import matplotlib.pyplot as plt
P,T,s,h,hf,hd,Ts,ss={},{},{},{},{},{},{},{}

Pamb=101325 
Tamb=15+273.15

Mexh=289           #Exhuast air Flow

eff=0.365          #Efficiency

Pel=97.7           #Electrical output

P[1]=Pamb
T[1]=Tamb
P[4]=Pamb
T[4]=538+273.15
P[2]=17*P[1]
P[3]=P[2]*0.95   #5% pressure loss in Combustion Chamber

eta_pC=0.80      #Polytropic Efficiency
eta_pT=0.75

gas='Air'
ch4='CH4'

Cp=cp.PropsSI('C','T',T[1],'P',P[1],gas)
Cv=cp.PropsSI ('CVMASS','T',T[1],'P',P[1],gas)

gamma=Cp/Cv
gamma1=(gamma-1)/gamma
gamma2=gamma1/eta_pC
gamma3=gamma1*eta_pT

LHV=50000

T[2]=T[1]*(17**gamma2)
T[3]=T[4]*((P[3]/P[4])**gamma3)
h[2]=cp.PropsSI('H','T',T[2],'P',P[2],gas)
h[3]=cp.PropsSI('H','T',T[3],'P',P[3],gas)
h[4]=cp.PropsSI('H','T',T[4],'P',P[4],gas)

deltaHAir=h[3]-h[2]
hf[3]=cp.PropsSI('H','T',T[3],'P',P[3],ch4)
hf[1]=cp.PropsSI('H','T',T[1],'P',P[1],ch4)
deltaHFuel=hf[3]-hf[1]
x=((deltaHAir+deltaHFuel)/1000)/LHV
M_air=Mexh/(1+x)
hd[3]=h[3]+x*hf[3]

h[1]=cp.PropsSI('H','T',T[1],'P',P[1],gas)
P_comp=M_air*((h[2]-h[1])/1000)
P_turb=Mexh*((hd[3]-h[4])/1000)
P_el=(P_turb-P_comp)*0.36*0.99

s[1]=cp.PropsSI('S','T',T[1],'P',P[1],gas)
s[2]=cp.PropsSI('S','T',T[2],'P',P[2],gas)
s[3]=cp.PropsSI('S','T',T[3],'P',P[3],gas)
s[4]=cp.PropsSI('S','T',T[4],'P',P[4],gas)
T1=np.linspace(T[2],T[3])
T2=np.linspace(T[1],T[4])

plt.plot(cp.PropsSI('S','T',T1,'P',P[2],gas),T1,'r',linewidth=1.5)
plt.plot(cp.PropsSI('S','T',T2,'P',P[4],gas),T2,'r',linewidth=1.5)
plt.plot([s[1],s[2]],[T[1],T[2]],'--r',linewidth=1.5)
plt.plot([s[3],s[4]],[T[3],T[4]],'--r',linewidth=1.5)


plt.xlabel('Entropy, s (kJ/kg/K)')
plt.ylabel('Temperature, T (K)')
plt.grid('on')
plt.title('Brayton Cycle T-s Real')


#Chart Data
plt.text(s[1]-50,T[1]-50,'(1)\nT={:.1f}\np={:.1f}kPa'.format(T[1],P[1]/1000),
    ha='right',backgroundcolor='white')

plt.text(s[2]+50,T[2]-120,'(2)\nT={:.1f}\np={:.1f}kPa'.format(T[2],P[2]/1000),
    ha='left',backgroundcolor='white')

plt.text(s[3]+150,T[3]-150,'(3)\nT={:.1f}\np={:.1f}kPa'.format(T[3],P[3]/1000),
    ha='left',backgroundcolor='white')

plt.text(s[4]+50,T[4],'(4)\nT={:.1f}\np={:.1f}kPa'.format(T[4],P[4]/1000),
    ha='left',backgroundcolor='white')

plt.text(5000,200,
"""$\dot{{m}}$ = {:.4f}kg/s
$p_r$={:.1f}
$\eta$={:.4f}
$\dot{{W}}_{{net}}$={:1}MW""".format(M_air,P[2]/P[1],eff,Pel),
    backgroundcolor='white')

Ts[2]=cp.PropsSI('T','S',s[1],'P',P[2],gas)
Ts[3]=cp.PropsSI('T','S',s[3],'P',P[3],gas)
ss[2]=cp.PropsSI('S','T',Ts[2],'P',P[2],gas)
ss[3]=cp.PropsSI('S','T',Ts[3],'P',P[3],gas)
ss[4]=cp.PropsSI('S','T',T[4],'P',P[4],gas)
T1s=np.linspace(Ts[2],Ts[3])
T2s=np.linspace(T[1],T[4])

plt.plot(cp.PropsSI('S','T',T1s,'P',P[2],gas),T1s,'g',linewidth=1.5)
plt.plot(cp.PropsSI('S','T',T2s,'P',P[4],gas),T2s,'g',linewidth=1.5)
plt.plot([s[1],ss[2]],[T[1],Ts[2]],'--g',linewidth=1.5)
plt.plot([ss[3],ss[4]],[Ts[3],T[4]],'--g',linewidth=1.5)


plt.xlabel('Entropy, s (kJ/kg/K)')
plt.ylabel('Temperature, T (K)')
plt.grid('on')
plt.title('Brayton Cycle T-s Ideal vs Real')


# Chart Data
#plt.text(s[1]-50,T[1]-50,'(1)\nT={:.1f}\np={:.1f}kPa'.format(T[1],P[1]/1000),
#    ha='right',backgroundcolor='white')
#
plt.text(ss[2]-20,Ts[2]+220,'(2)\nT={:.1f}\np={:.1f}kPa'.format(T[2],P[2]/1000),
    ha='left',backgroundcolor='white')
#
#plt.text(ss[3]+50,Ts[3]+50,'(3)\nT={:.1f}\np={:.1f}kPa'.format(T[3],P[2]/1000),
#    ha='left',backgroundcolor='white')
#
#plt.text(s[4]+50,T[4],'(4)\nT={:.1f}\np={:.1f}kPa'.format(T[4],P[4]/1000),
#    ha='left',backgroundcolor='white')

#plt.text(5200,250,
#"""$\dot{{m}}$ = {:.4f}kg/s
#$p_r$={:.1f}
#$\eta$={:.4f}
#$\dot{{W}}_{{net}}$={:1}MW""".format(M_air,P[2]/P[1],eff,Pel),
#    backgroundcolor='white')
#plt.savefig('Ideal.png')