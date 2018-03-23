import ttim
from ttim.aquifer_parameters import param_3d, param_maq
import numpy as np
def NonDelay(ml,drawdown,z,LN,HC,Sfe,Sfv):#,Com,ComE,ComV):
    tmin, tmax = ml.tmin, ml.tmax
    print(tmin)
    z = np.array(z)
    nlay = len(z) - 1
    LN_dict = {}
    HC_dict = {}
    for i in range(len(z)-1):
        LN_dict[i] = LN[i]
        HC_dict[i] = HC[i]


    print(LN_dict)

    def del_Estress(drawdown,Pw=62.42796529, g=240177996288.):
        """
        theta = vertical effective stress (positive for increase)
        :param pw: 62.42796529 lb/ft^3
        :param g: 32.17405 ft/sec^2
        :param h:
        :return:
        """
        theta = []
        for lay in range(nlay):
            theta.append(-Pw * g * drawdown[lay])
        return theta

    theta = del_Estress(drawdown,Pw=62.42796529, g=240177996288.)

    pc_dd = []
    for lay in range(nlay):
        pc_dd = z[lay] - HC

    theta_pc = del_Estress(drawdown=pc_dd,Pw=62.42796529, g=240177996288.) # prconsolidation
    b = []
    for lay in range(nlay):
        for Sk in range(len(Sfe)):
            for i in range(len(drawdown)):
                if theta[lay][i] < theta_pc[lay]:
                    S_k = Sfe[lay][Sk]
                elif theta[lay][i] >= theta_pc[lay]:
                    S_k = Sfv[lay][Sk]


                b.append(S_k * drawdown[i])

    comp = np.zeros((nlay))
    print(b[0])
    for lay in range(nlay):
        for i in range(len(b[lay])):
          comp[lay] += b[lay][i]



    return comp








