__author__ = 'ecpease'

import numpy as np
import matplotlib.pyplot as plt
import ttim
class model():

    def __init__(self, nlay, n, Z, alpha, beta, k, dh, sigma_total, pf, p, g, sigma_eff):
        """
        :param n: porosity of the layer
        :param alpha: compressibility of the layer medium
        :param beta: compressibility of water
        :param k: permeability
        :param dh: change in hydraulic head with shape of (nlay, number of times)
        :param Z: is inital formation thickess with shape nlay
        :param sigma_total: total stress (Terzaghi Equation)
        :param pf: pore-fluid pressure
        :param p: density of water 
        :param g: force of gravity
        :param sigma_eff: effective stress
        :param db: change in layer thickness 
        :param b: initial layer thickness
        :param d_sigma_eff = change in the effective stress
        :param Sske: elastic skeletal specific storage
        :param Sskv: inelastic skeletal specific storage
        :param sigma_preconsol: preconsolidation stress
        """
        
        self.p = 1000 # kg/m3
        self.g = 9.8  # m/s2
        self.Z = Z
        self.beta = beta
        self.dh = dh
        self.b = b
        self.d_sigma_eff = d_sigma_eff
        self.pf = pf
        
        

    def formation_compaction(self):
            Zt = self.Z
            dh = np.array(self.dh)
            nlay = self.nlay
            deltaZ = np.zeros((dh.shape))
            ZL = np.zeros((dh.shape)) # thickness of layer with shape nlay, number of times
            for i, deltahead in enumerate(dh):
                Zi = self.p * self.g * deltaZ[i] * deltahead # shape nlay
                ZL[i] = Zi
                
                if i == 0:
                    deltaZ[i] = np.zeros(nlay)
                else:
                    deltaZ[i] = ZL[i] - ZL[i-1]
            return deltaZ
        
    def landSubsidence(self, deltaZ):
            Zt = self.Z
            dh = np.array(self.dh)
            nlay = self.nlay
            
            L = []
            
            for i in range(len(deltaZ[1])):
                L.append(sum(deltaZ[i]))
                
            return L
                
        

def head(pw, p, g, he):
    """
    Get total head
    :param pw: pore-fluid pressure
    :param p: water density
    :param g: gravity
    :param he: elevation head
    """
     
    h = (pw / p * g) + he
    return h
       

def sigma_unconf(sigma_eff, pw, g, dh, n):
    """
    Calculate effective stress for unconfined aquifer
    :param pw: pore-fluid pressure
    :param p: water density
    :param g: gravity
    :param he: elevation head
    """
    
    sigma_eff = -pw * g * (1 - n + nw) * dh
    return sigma_eff


def sigma_conf(sigma_eff, pw, g, dh):
    """
    Calculate effective stress for confined aquifer
    :param pw: pore-fluid pressure
    :param p: water density
    :param g: gravity
    :param he: elevation head
    """
    
    sigma_eff = -pw * g * dh
    return sigma_eff


def alpha_3d(dV, V, d_sigma_eff):
    """
    Calculate sediment compaction
    :param dV: change in volume 
    :param V: initial volume
    :param d_sigma_eff = change in the effective stress
    """
    
    alpha = - (dv / V) / d_sigma_eff
    return alpha


def alpha(db, b, d_sigma_eff):
    """
    Calculate sediment compaction; given no horizonatal displacements
    :param db: change in layer thickness 
    :param b: initial layer thickness
    :param d_sigma_eff = change in the effective stress
    """
    
    alpha = - (db / b) / d_sigma_eff
    return alpha


# def alpha_compressibility(pw, g, alpha, Ssk, b):
#     pw * alpha * g * b = Ssk * b # = Sk = db / dh
    

def specific_storage(p, g, alpha, C, beta, n):
    """
    Specific storage equation
    :param p: water density
    :param g: gravity
    :param alpha: sediment compressibility
    :param C: 
    :param beta: water compressibility
    :param n: porosity
    """
    
    p = 1000 # kg/m3
    g = 9.81 # m/s2

    alpha = 0.434 * C * (1 - n) / total_sigma

    Ss = p * g * (alpha + (n * beta))
    return Ss


def db(Sk, dh):
    db = Sk * dh
    return db


def Ssk(Sske, Sskv, sigma_preconsol, eff_stress):
    """
    Calculate skeletal specific storage when effective stress exceeds preconsolidation stress
    :param Sske: elastic skeletal specific storage
    :param Sskv: inelastic skeletal specific storage
    :param sigma_preconsol: preconsolidation stress
    """
    
    if eff_stress < sigma_preconsol:
        Ssk = Sske
    elif eff_stress >= sigma_preconsol:
        Ssk = Sskv
        
    return Ssk


def Permability(ko, n, no):
    """
    Calculate permeability change as a function of porosity
    Ko = initial permeability
    n = porosity
    m = coefficient (for MX City it is 3)
    """
    
    m = 3
    k = ko * (n * (1 - no))/(no(1-n)) ^ m

    return K



# def subsidence(p, g, dz, alpha, dh):

#     """
#     Get total subsidence (pg. 48, eq 8 - Rivera et al 1991)

#     """

#     sums = []

#     for i in nlay:
#         sums.append(p * g * dz * alpha * dh)
#     sums = sum(sums)
        
        
        
# def terzaghi(sigma_total, pf):
#     """
#     Terzaghi equation; total_stress CAN change through time

#     """
#     eff_stress = tot_stress - pfp
#     return eff_stress


# def terzaghi_constant(eff_stress, pfp):
#     """
#     Terzaghi equation; total_stress does NOT change through time

#     """
#     eff_stress = - pfp
#     return eff_stress