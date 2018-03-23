import ttim
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

sys.path.append(os.path.join('..'))
import subtim.compaction as sub






# from subtim import compaction


z = [200,100,0]
kaq = [1e-6,10] # conductivity
Saq = [1e-10,1e-7] # Storage coefficient
kzoverkh = 1


ten_years = 10 * 365.25


ml = ttim.Model3D(kaq,z,Saq,kzoverkh,tmin=1e-6,tmax=ten_years)


print(ml.tmin)


Qgpm = 100.0
Qcfd = (Qgpm / 7.4801) * (60*24)

tsandQ = [(0,0),(365.25,Qcfd)]

res = .1
x,y = 0.,0.
rw = (12.75/12)*.5

well = ttim.Well(ml,x,y,rw,tsandQ,res,layers=1,label='well')


ml.solve()



time = np.arange(0,365.26*10,365.25)
h = ml.head(x,y,time,[0,1])




# drawdown = h[0][0]-h[0]
drawdown = [h[0][0] - h[0],h[1][1] - h[1]]


LN = [[(180,165),(140,120)],[(90,80),(70,60),(65,50)]]
HC = [200,100]

ske = 5e-6
skv = 3e-4
Sfe = [[ske,ske],[ske,ske,ske]]
Sfv = [[skv,skv],[skv,skv,skv]]

comp = sub.NonDelay(ml,drawdown,z,LN,HC,Sfe,Sfv)#,Com,ComE,ComV)

print(comp)
# print(b)
# print(b.shape)
#
# for i in range(len(b)):
#     print(b[i])


exit()

def compress_coef(dp,dv,V):
    """

    :param dp:
    :param dv:
    :param V:
    :return: a
    """
    a = dv/V * (1/dp)
    return a

def ff(delh,a,b):
    """

    :param delh: change in hydraulic head in the aquifer
    :param a: Aquifer compressibility coef, (calculated from one or more interbeds)
    :param b:
    :return: delb
    """
    delb = delh * a *b
    return delb



def pore_p(h,he,g=32.17405,Pw=62.42796529):
    """

    :param h: head ft
    :param he: elevation head relative to datum
    :param g: gravity in ft/s^2
    :param Pw: lb/ft^3
    :return:
    """
    # 32.17405 feet per second per second.
    # h = (p)/(Pw*g) + he
    p = (h - he)/(Pw*g)
    return p

def v_Estress(pw,g,h):
    """
    theta = vertical effective stress (positive for increase)
    :param pw:
    :param g:
    :param h:
    :return:
    """
    sh= h[0]
    delh = []
    for i in h:
        delh.append(sh-i)
    theta = -pw * g * delh
    return theta





fig, ax = plt.subplots()
ax.plot(time/365.25,h[0])
ax.grid()
ax.set_xlabel('Years')
ax.set_ylabel('Drawdown')



head_c = ml.contour(win=[-10,10,-10,10],t=time.max(),layers=1)








plt.show()

















