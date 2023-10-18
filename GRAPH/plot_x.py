#
#  plot_x.py
#
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, show
from numpy import arange
import numpy as np

#x = arange(0,41,1)
#x = arange(-2,2.1,0.1)

#x = arange(0,145,1)
#x = arange(96,241,1)
#x = arange(1,241,1)
x = arange(0,792,1)
y = np.loadtxt("v_x.dat")

#y2 = np.loadtxt("y1b.dat")
#y3 = np.loadtxt("y1c.dat")
#y3 = np.loadtxt("y4a.dat")
#y4 = np.loadtxt("y4b.dat")


fig = figure(1)

ax1 = fig.add_subplot(111)

ax1.plot(x,y)

#xmajor_ticks=np.arange(96,240,48)
#xmajor_ticks=np.arange(1,504,60)
xmajor_ticks=np.arange(0,792,72)
ax1.set_xticks(xmajor_ticks)
#ax1.grid(which='x')
ax1.grid(axis='x')

#ax1.plot(x, y2,label='1st b')
#ax1.plot(x, y3,label='1st c')
#ax1.plot(x, y4,label='4th b')

#ax1.legend(loc='upper right')
#ax1.grid(True)

#ax1.set_xlim([0,792])
ax1.set_xlim([300,492])
ax1.set_ylim([-0.1,1.8])  


plt.suptitle('3 Generations with half of vertical resolution' ,ha='center',fontsize=14)
plt.title('Line cross-section in x-direction',ha='center',fontsize=11)

#ax1.text(20,1., "$w=(1,1,1,1,1,401)$ \n$p=1.3$ \n$a_1=1$ \n$a_2=1$ \n$(h_x,h_y,h_z)=(6,6,6)$ \n$\psi_m=8.82$ \n$\lambda_h=83.1$", bbox=dict(facecolor='yellow', alpha=0.3))




plt.xlabel('x')
#plt.ylabel('y')
plt.savefig('p3g_x.png',bbox_inches='tight')

plt.show()
