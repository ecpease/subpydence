import ttim
from ttim.aquifer_parameters import param_3d, param_maq
import numpy as np


def NonDelay(drawdown,z,LN,HC,Sfe,Sfv):
    """

    :param drawdown:
    :param z: elevations of each each geological formation, top, botm1, botm1 etc... array of len nlay +1
    :param LN: Number of delay beds in each geological unit. array of len nlay
    :param HC: preconsolodation head, array of len nlay, contains the initial preconsolodation head.
    :param Sfe: skeletal elastic storage coefficient
    :param Sfv: skeletal inelastic storage coefficient
    :return:
    """
    z = np.array(z)
    nlay = len(z) - 1 # get number of layers

    thick = []
    for t in range(1,len(z)):
        thick.append(z[t-1]-z[t])
    # print(thick)
    # exit()
    i,elements = 0,['LN','HC','Sfe','Sfv']
    for element in [LN,HC,Sfe,Sfv]:
        if nlay != len(element):
            raise ValueError(f'nlay: {nlay} != {elements[i]}: {len(element)}.')
        i+=1
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
        return np.array(theta)

    theta = del_Estress(drawdown,Pw=62.42796529, g=240177996288.)
    pc_dd = []
    for lay in range(nlay):
        pc_dd = z[lay] - HC

    theta_pc = del_Estress(drawdown=pc_dd,Pw=62.42796529, g=240177996288.) # preconsolidation
    # b = []
    # print(theta_pc)
    # for lay in range(nlay):
    #     for i in range(len(drawdown)):
    #         for Sk in range(LN[lay]):
    #             if theta[lay][i] < theta_pc[lay]:
    #                 S_k = Sfe[lay][Sk]
    #             elif theta[lay][i] >= theta_pc[lay]:
    #                 S_k = Sfv[lay][Sk]
    #             b.append(S_k * drawdown[i]*thick[lay])

    b = np.zeros(drawdown.shape)
    # print(b.shape)
    # exit()
    # S_k = np.zeros(drawdown.shape[0])
    # print(Sfe)
    # exit()
    ii = 0
    for lay in range(nlay):
        for i in range(len(drawdown[lay])):
            if theta[lay][i] < theta_pc[lay]:
                # S_k = Sfe[lay]
                # S_k[lay] = Sfe[lay]
                S_k = np.array(Sfe[lay]).mean()  # need to change this to composite mean
            elif theta[lay][i] >= theta_pc[lay]:
                # S_k = Sfv[lay]
                # S_k[lay] = Sfv[lay]
                S_k = np.array(Sfv[lay]).mean()
            # b.append(S_k * drawdown[i])
            # print(S_k)
            # exit()
            thing = S_k * drawdown[lay][i] * thick[lay]
            b[lay][i] = thing

    comp = []
    b = np.array(b)
    print(b.shape)
    # exit()

    for lay in range(nlay):
        comp.append(np.zeros((len(b[0]))))
        comp[lay] += b[lay] * LN[lay]

    return comp

def NonDelayGrid(drawdownGrid,z,LN,HC,Sfe,Sfv): #,Com,ComE,ComV):
    z = np.array(z)
    nlay = len(z) - 1
    def del_Estress(drawdownGrid, Pw=62.42796529, g=240177996288.):
        """

        theta = vertical effective stress (positive for increase)
        :param pw: 62.42796529 lb/ft^3
        :param g: 32.17405 ft/sec^2
        :param h:
        :return:
        """
        theta = [] #np.zeros((nlay))
        for lay in range(nlay):
            theta.append(np.zeros(nlay))
            theta[lay] = -Pw * g * drawdownGrid[lay]
        return theta

    theta = del_Estress(drawdownGrid, Pw=62.42796529, g=240177996288.)
    pc_dd = []
    for lay in range(nlay):
        pc_dd.append(z[lay] - HC)
    fig, ax = plt.subplots()
    plt.imshow()
    theta_pc = del_Estress(drawdownGrid=pc_dd, Pw=62.42796529, g=240177996288.)  # pre consolidation
    b = []
    print(theta_pc)
    for lay in range(nlay):
        for i in range(len(drawdownGrid)): # for each time
            for Sk in range(len(LN)):
                # if theta[lay][i] < theta_pc[lay]:
                #     S_k = Sfe[lay][Sk]
                print(theta[lay][i].shape)
                print(theta_pc[lay].shape)
                Sfe_locs = np.where(theta[lay][i] < theta_pc[lay])
                Sfv_locs = np.where(theta[lay][i] >= theta_pc[lay])
                # elif theta[lay][i] >= theta_pc[lay]:
                #     S_k = Sfv[lay][Sk]
                b.append(S_k * drawdownGrid[i])
    comp = []
    b = np.array(b)
    for lay in range(nlay):
        comp.append(np.zeros((len(b[0]))))
        comp[lay] += b[lay] * LN[lay]

    return compGrid

