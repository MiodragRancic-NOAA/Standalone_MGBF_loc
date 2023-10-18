#
#  plot_z.py
#
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
from numpy import arange
import numpy as np

x = arange(0,65,1)
y = np.loadtxt("v_z.dat")


fig = figure(1)

ax1 = fig.add_subplot(111)

ax1.set_xlim([-0.1,1.8])
ax1.set_ylim([1,65])
ax1.grid(True)
#ax1.grid(axis='z')
ymajor_ticks=np.arange(1,65,5)

ax1.plot(y,x)

#ax1.set_xticks(xmajor_ticks)

#ax1.plot(x, y2,label='1st b')
#ax1.plot(x, y3,label='1st c')
#ax1.plot(x, y4,label='4th b')

#ax1.legend(loc='upper right')



plt.suptitle('3 Generations with half of vertical resolution' ,ha='center',fontsize=14)
plt.title('Line cross-section in z-direction',ha='center',fontsize=11)

#ax1.text(20,10, "$w=(1,1,1,2048)$ \n$a_1=1$ \n$a_2=8$ \n$\psi_m=14.67$ \n$\lambda_h=116.32$", bbox=dict(facecolor='yellow', alpha=0.3))




plt.xlabel('x')
#plt.ylabel('y')
plt.savefig('p3g_z.png',bbox_inches='tight')

plt.show()
