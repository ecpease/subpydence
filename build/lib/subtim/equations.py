import numpy

def terzaghi(eff_stress, tot_stress, pfp):
    """
    Terzaghi equation (effective stress = total stress - porefluid pressure; total_stress CAN change through time

    """
    eff_stress = tot_stress - pfp
    return eff_stress
   

def terzaghi_constant(eff_stress, pfp):
    """
    Terzaghi equation (effective stress = - porefluid pressure; total_stress does NOT change through time

    """
    eff_stress = - pfp
    return eff_stress
   

def specific_storage(pw, g, a, Bw, n, ssk, ssw):
    """
    Specific storage equation
    
    """
    ssk = pw * g * a
    ssw = pw * g * n * Bw
    
    return eff_stress

def NonDelay(drawdown, z, LN, HC, Sfe, Sfv):
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

    thickness = []
    for t in range(1, len(z)):
        thickness.append(z[t-1]-z[t])

    i, elements = 0, ['LN','HC','Sfe','Sfv']
    for element in [LN, HC, Sfe, Sfv]:
        if nlay != len(element):
            raise ValueError(f'nlay: {nlay} != {elements[i]}: {len(element)}.')
        i+=1

