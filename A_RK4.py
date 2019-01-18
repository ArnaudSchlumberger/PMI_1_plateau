import numpy as np
import matplotlib.pyplot as plt

def equa_plateau_1(t,y1,y0,k,m):
    g = 9.81 #m.s**-1
    return -g-(k/m)*(y1-y0)
    
    
def RK4(equa,t,f,y1,y0,dy1,k,m,h=0.01):
    k1 = equa(t,y1,y0,k,m)
    k2 = equa(t+(h/2),y1+(h/2)*dy1,y0,k,m)
    k3 = equa(t+(h/2),y1+(h/2)*dy1+(h**2/4)*k1,y0,k,m)
    k4 = equa(t+h,y1+h*dy1+(h**(2)/2)*k2,y0,k,m)
    
    new_y1=y1 + h*dy1 + ((h**2)/6)*(k1+k2+k3)
    new_dy1=dy1+(h/6)*(k1+2*k2+2*k3+k4)
    
    return new_y1, new_dy1
    
liste_t = [0]
liste_y0 = [0]
y = [0]
dy = [0]
f = 5 #Hz
tf = 10 #s
i = 0
h = 1/(10*f)
t=0
k=1000
m=1
while t<tf:
    y0_tmp = np.sin(np.pi*2*f*t)
    y1_tmp,dy1_tmp=RK4(equa_plateau_1,t,f,y[i],y0_tmp,dy[i],k,m,h)
    t=t+h
    i=i+1
    liste_y0.append(y0_tmp)
    y.append(y1_tmp)
    dy.append(dy1_tmp)
    liste_t.append(t)
plt.subplot(2,1,1)
plt.plot(liste_t,y)
plt.subplot(2,1,2)
plt.plot(liste_t,liste_y0)
plt.show()