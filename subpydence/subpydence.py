__author__ = 'ecpease'

import numpy as np
import matplotlib.pyplot as plt
import ttim

class model():

    def __init__(self, nlay, h, z, Ssk,pw=1000., g=9.81):
        """
        :param n: porosity of the layer
        :param alpha: compressibility of the layer medium
        :param beta: compressibility of water
        :param k: permeability
        :param dh: change in hydraulic head with shape of (nlay, number of times)
        :param Z: is inital formation thickess with shape nlay
        :param sigma_total: total stress (Terzaghi Equation)
        :param pf: pore-fluid pressure
        :param pw: density of water 
        :param g: force of gravity
        :param sigma_eff: effective stress
        :param db: change in layer thickness 
        :param b: initial layer thickness
        :param d_sigma_eff = change in the effective stress
        :param Sske: elastic skeletal specific storage
        :param Sskv: inelastic skeletal specific storage
        :param sigma_preconsol: preconsolidation stress
        """
        
        self.nlay = nlay
        self.pw = 1000 # kg/m3
        self.g = 9.81  # m/s2
        self.dh = np.diff(h)
        self.h = h
        self.b = -np.diff(z)
        self.Ssk = Ssk
        self.ntimes = len(h[0])
        self.btm = z[-1]
        self.z = z
 
        
        
    """
    Calculate pore-fluid pressure component 
    h = (p/pw*g) + z
    """
    
    def pressure(self):
        """
        Get pore-fluid pressure
        :param p: pore-fluid pressure
        :param pw: water density
        :param g: gravity
        :param btm: elevation head (bottom of the aquifer will be our datum)
        """
        #         p = (h-btm) * self.pw * self.g
        
        p = np.zeros(self.h.shape)
        
        for lay in range(self.nlay):
            for i in range(self.ntimes):
                p[lay][i] = (self.h[lay][i] - self.z[lay+1]) * self.pw * self.g
        return p        

    def d_eff_stress(self):
        """
        Get change in effective stress, given pore-fluid pressure
        deff_stress = dtotal_stress - dp
        dtotal_stress = 0 (no additional loading to the land surface ie construction)
        deff_stress = -dp
        """
        
        d_eff_stress = -self.dh * self.pw * self.g
        eff_stress_list = np.zeros((self.nlay,self.ntimes))
        for lay in range(nlay):
            for i in range(0, len(self.ntimes)):
                d_eff_stress = -self.pw * self.g * self.h[lay][i] -self.h[lay][i-1]
                eff_stress_list[lay][i] = d_eff_stress        
        
        return np.diff(eff_stress_list)
    
    
    def db(self):
#         db = self.Ssk * self.b * self.dh
        
        db = np.zeros((self.nlay, self.ntimes -1))
        nb = np.zeros((self.nlay, self.ntimes)) # new thickness
        for lay in range(self.nlay): # first time
#             print(self.Ssk[lay], self.b[lay], self.dh[lay][0] )
            db[lay][0] = self.Ssk[lay] * self.b[lay] * -self.dh[lay][0] # change in thickness for first time
            nb[lay][0] = self.b[lay] # init thickness
            nb[lay][1] = self.b[lay] - db[lay][0] # thickness after first time
        
        for lay in range(self.nlay):
            for i in range(0,self.ntimes-1):
                db[lay][i] = self.Ssk[lay] * (nb[lay][i] - db[lay][i-1]) * -self.dh[lay][i] # negative since dh is already negative
                nb[lay][i+1] = nb[lay][i] - db[lay][i]
        
        return db,nb
    
    
    def alpha(self, db, b, d_eff_stress):
        """
        Calculate sediment compaction; given no horizonatal displacements
        :param db: change in layer thickness 
        :param b: initial layer thickness
        :param d_sigma_eff = change in the effective stress
        """

        alpha = - (db / self.b) / d_eff_stress
        return alpha
    
    
    

    
    
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
                
        

    def head(p, pw, g, he):
        """
        Get total head
        :param p: pore-fluid pressure
        :param pw: water density
        :param g: gravity
        :param he: elevation head
        """

        h = (p / pw * g) + he
        return h




    def sigma_unconf(pw, g, dh, n):
        """
        Calculate effective stress for unconfined aquifer
        :param pw: pore-fluid pressure
        :param p: water density
        :param g: gravity
        :param he: elevation head
        """

        sigma_eff = -pw * g * (1 - n + nw) * dh
        return sigma_eff


    def sigma_conf(pw, g, dh):
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





    # def alpha_compressibility(pw, g, alpha, Ssk, b):
    #     pw * alpha * g * b = Ssk * b # = Sk = db / dh


    def specific_storage(pw, g, alpha, C, beta, n):
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

        Ss = pw * g * (alpha + (n * beta))
        return Ss


    # def db(Sk, dh):
    #     db = Sk * dh
    #     return db






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