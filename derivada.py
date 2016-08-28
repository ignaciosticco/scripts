# Sacado de https://scipy.github.io/old-wiki/pages/Cookbook/Matplotlib/LaTeX_Examples.html

import matplotlib.pyplot as plt
import pylab
import numpy as np
import matplotlib.pyplot as plt
import math
from pylab import arange,pi,sin,cos,sqrt,savefig

# a dos columnas: 3+3/8 (ancho, in)
# a una columna : 6.5   (ancho  in)

golden_mean = (math.sqrt(5)-1.0)/2.0        # Aesthetic ratio
fig_width = 3+3/8  			    # width  in inches
fig_height = fig_width*golden_mean          # height in inches
fig_size =  [fig_width,fig_height]

params = {'backend': 'ps',
          'axes.titlesize': 8,
          'axes.labelsize': 9,
          'axes.linewidth': 0.5, 
          'axes.grid': True,
          'axes.labelweight': 'normal',  
          'font.family': 'serif',
          'font.size': 8.0,
          'font.weight': 'normal',
          'text.color': 'black',
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'text.usetex': True,
          'legend.fontsize': 8,
          'figure.dpi': 300,
          'figure.figsize': fig_size,
          'savefig.dpi': 300,
         }

pylab.rcParams.update(params)
#################### DATA ##############################


data = np.genfromtxt('gap_vs_te_40x40.txt', delimiter = '   ')
gap=data[:,0]
te=data[:,1]
yerr=data[:,2]
i=1
derivada=[]
while i<len(te):
	print(i)
	cuenta=(te[i]-te[i-1])/(gap[i]-gap[i-1])
	cuenta=round(cuenta,2)
	derivada.append(cuenta)
	i+=1
print(derivada)


#################### PLOT specification ##############################
pylab.figure(1)
pylab.clf()
plt.plot(gap[1:],derivada,'wo',markersize=3,zorder=3) 
plt.plot(gap[1:],derivada,'k',lw=0.7,zorder=2)     
#plt.errorbar(gap3[::2],te3[::2],yerr3[::2],linestyle='-',marker='.',color='k')    
#pylab.xticks(np.arange(0,24,4))
#pylab.yticks(np.arange(200,410,66))
pylab.xlabel('Gap~(m)')
pylab.ylabel('Derivada')
#pylab.legend()
pylab.ylim(-20,20)
#pylab.xlim(0, 20)
#plt.legend(loc=4)
pylab.grid(True)   
pylab.savefig('derivada_961_v4_zoom.eps', format='eps', dpi=600, bbox_inches='tight')


