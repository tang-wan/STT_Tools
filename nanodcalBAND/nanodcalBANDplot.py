import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio

from STT_Tool import Tools
color = Tools.ColorList()

# >>>>>>>>>> <<<<<<<<<<
class nanodcalBANDplot():
    def __init__(self, B_Filepath):
        self.BandData = scio.loadmat(B_Filepath)['data'][0][0]

    def Read_AllData(self):
        data = self.BandData
        
        self.a_vec = data[1][0]
        self.b_vec = data[1][1]
        self.c_vec = data[1][2]

        # SymmtryKPoints
        Symm_xlabel = []
        for j in range(len(data[4][0])):
            Symm_xlabel.extend([i for i in data[4][0][j][0]])

        Symm_xtick  = np.array(data[5][0])
        self.Symm_xtick  = Symm_xtick
        
        Symm_xlabel      = np.transpose(np.array(Symm_xlabel)[:,0])
        self.Symm_xlabel = Symm_xlabel

        # fermiEnergy
        Fermi_E = data[6][0]
        self.Fermi_E = Fermi_E

        # BandGap
        BandGap = data[10][0][0]

        # BandEnergy
        BandEnergy = (data[7][0][0]-Fermi_E+BandGap/2)*27.211386245981
        Kpath = np.linspace(1, len(BandEnergy[0]), len(BandEnergy[0]))
        self.BandEnergy = BandEnergy
        self.Kpath = Kpath
        
        LABELoutput = (Symm_xlabel, Symm_xtick)
        BANDoutput  = (Kpath,       BandEnergy)
        return BANDoutput, LABELoutput
        
    
    def Print_info(self):
        print()
        print('Lattice vector is: ', '\n', 
              self.a_vec, '\n', 
              self.b_vec, '\n', 
              self.c_vec, '\n')
        print()
        print(self.Symm_xlabel)
        print(self.Symm_xtick)
        print()
        print('FermiEnergy = ', self.Fermi_E)
        print()
        print('Band Number = ', np.shape(self.BandEnergy))
