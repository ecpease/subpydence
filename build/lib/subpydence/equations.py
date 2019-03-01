def terzaghi(tot_stress, pfp):
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
   

# def specific_storage(pw, g, a, Bw, n, ssk, ssw):
#     """
#     Specific storage equation
    
#     """
#     ssk = pw * g * a
#     ssw = pw * g * n * Bw
    
#     return eff_stress

def specific_storage(p, g, alpha, C, beta, n):
    """
    Specific storage equation
    
    """
    p = 1000 # kg/m3
    g = 9.81 # m/s2
    
    alpha = 0.434 * C * (1 - n) / total_sigma
    
    Ss = p * g * (alpha + (n * beta))
    return Ss
    

    
def Permability(Ko, n, no, m):
    """
    Calculate permeability change as a function of porosity
    Ko = initial permeability
    n = porosity
    m = coefficient (for MX City it is 3)

    """
    m=3
    K = Ko * (n * (1 - no))/(no(1-n)) ^ m
    
    return K
 
    
    
def subsidence(p, g, dZ, alpha, dh):

    """
    Get total subsidence (pg. 48, eq 8 - Rivera et al 1991)
    
    """
    
    sums = []
    
    for i in nlay:
        sums.append(p * g * dZ * alpha * dh)
        
    sums = sum(sums)
    return sums
    
    
# def NonDelay(drawdown, z, LN, HC, Sfe, Sfv):
#     """

#     :param drawdown:
#     :param z: elevations of each each geological formation, top, botm1, botm1 etc... array of len nlay +1
#     :param LN: Number of delay beds in each geological unit. array of len nlay
#     :param HC: preconsolodation head, array of len nlay, contains the initial preconsolodation head.
#     :param Sfe: skeletal elastic storage coefficient
#     :param Sfv: skeletal inelastic storage coefficient
#     :return:
#     """
#     z = np.array(z)
#     nlay = len(z) - 1 # get number of layers

#     thickness = []
#     for t in range(1, len(z)):
#         thickness.append(z[t-1]-z[t])

#     i, elements = 0, ['LN','HC','Sfe','Sfv']
#     for element in [LN, HC, Sfe, Sfv]:
#         if nlay != len(element):
#             raise ValueError(f'nlay: {nlay} != {elements[i]}: {len(element)}.')
#         i+=1

