import numpy as np
import matplotlib.pyplot as plt

class BandPlot():
    def __init__(self, atom, PLBasis, SOC, ndBasis, r):
        self.atom    = atom
        self.PLBasis = PLBasis
        self.SOC     = SOC
        self.ndBasis = ndBasis
        self.r       = r

        self.nd_path = f'{atom}/{PLBasis}_r={r}/Band_single/{ndBasis}_{PLBasis}_{SOC.lower()}_band.txt'
        self.vp_path = f'{atom}/VASP_r={r}/{PLBasis}_{SOC.lower()}_zx0_mz0_band'

        if SOC == 'soc':
            re_SOC = 'ncl'
        elif SOC == 'ncl':
            re_SOC = 'soc'

        self.re_SOC = re_SOC
        self.re_vp_path = f'{atom}/VASP_r={r}/{PLBasis}_{re_SOC.lower()}_zx0_mz0_band'
        self.re_nd_path = f'{atom}/{PLBasis}_r={r}/Band_single/{ndBasis}_{PLBasis}_{re_SOC}_band.txt'
        
    def NdPlot(self, tick, label):
        nd_path = self.nd_path
        
        BandData = np.loadtxt(f'{nd_path}')
        
        plt.figure(figsize=(8, 6))
        for data in BandData:
            py = data
            px = np.linspace(0, 1, len(py))

            plt.plot(px, py, color='#3976af')
            plt.xticks(ticks=tick, labels=label)
            plt.ylim(ymin=-10, ymax=10)
        
        plt.grid('--')
        plt.show()
    
    def VpPlot(self, tick, label):
        vp_path = self.vp_path

        BandData = np.loadtxt(f'{vp_path}')
        
        plt.figure(figsize=(8, 6))
        for i in range(1, len(BandData[0])):
            px = BandData[:,0]
            px = px/np.max(px)
            py = BandData[:,i]
            
            plt.plot(px, py, color='#3976af')
            plt.xticks(ticks=tick, labels=label)
            plt.ylim(ymin=-10, ymax=10)
            
        plt.grid('--')
        plt.show()

    def NdVpPlot(self, tick, label):

        PLBasis = self.PLBasis
        SOC = self.SOC
        ndBasis = self.ndBasis

        Nd_Data = np.loadtxt(f'{self.nd_path}')
        Vp_Data = np.loadtxt(f'{self.vp_path}')

        index = 0
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_title(f'{ndBasis}_{PLBasis}_{SOC.lower()}_band\nNanodcal vs. Vasp', fontsize=14)
        for data in Nd_Data:
            py = data
            px = np.linspace(0, 1, len(py))
            if index == 0:
                ax.plot(px, py, 
                        color='#738a79', 
                        linestyle='--',
                        label='nd')
            else:
                ax.plot(px, py, 
                        color='#738a79', 
                        linestyle='--')
            index = index + 1
        
        for i in range(1, len(Vp_Data[0])):
            px = Vp_Data[:,0]
            px = px/np.max(px)
            py = Vp_Data[:,i]
            if i == 1:
                ax.plot(px, py, 
                         color='#d28063', 
                         linestyle=':',
                         label='vp')
            else:
                ax.plot(px, py, 
                         color='#d28063', 
                         linestyle=':')
        ax.legend(loc='best')
        ax.set_xticks(ticks=tick, labels=label)
        ax.set_ylim(ymin=-10, ymax=10)
        ax.grid('--')

        plt.tight_layout()

    def NdNdPlot_SocNcl(self, tick, label):
        
        atom    = self.atom
        PLBasis = self.PLBasis
        SOC     = self.SOC
        re_SOC  = self.re_SOC
        ndBasis = self.ndBasis
        r       = self.r

        Nd_Data = np.loadtxt(f'{self.nd_path}')
        re_Nd_Data = np.loadtxt(f'{self.re_nd_path}')

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_title(f'{ndBasis}_{PLBasis}_band\nsoc vs. ncl', fontsize=14)
        # for data in Nd_Data:
        for i in range(len(Nd_Data)):
            py = Nd_Data[i]
            px = np.linspace(0, 1, len(py))
            repy = re_Nd_Data[i]
            if i == 0:
                ax.plot(px, py, 
                        color='#738a79', 
                        linestyle='--',
                        label=f'{SOC}')
                
                ax.plot(px, repy,
                        color='#b57fb3', 
                        linestyle=':',
                        label=f'{re_SOC}'
                        )
            else:
                ax.plot(px, py, 
                        color='#738a79', 
                        linestyle='--')
                ax.plot(px, repy,
                        color='#b57fb3', 
                        linestyle=':'
                        )
        ax.legend(loc='best')
        ax.set_xticks(ticks=tick, labels=label)
        ax.set_ylim(ymin=-10, ymax=10)
        ax.grid('--')

        plt.tight_layout()
    
    def NdNdPlot_r(self, tick, label):
        
        atom    = self.atom
        PLBasis = self.PLBasis
        SOC     = self.SOC
        ndBasis = self.ndBasis
        r       = self.r

        Nd_Data = np.loadtxt(f'{self.nd_path}')
        if r == 0.5:
            re_r = 0
        elif r == 0:
            re_r = 0.5
        re_Nd_Data = np.loadtxt(f'{atom}/{PLBasis}_r={re_r}/Band_single/{ndBasis}_{PLBasis}_{SOC}_band.txt')

        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_title(f'{ndBasis}_{PLBasis}_{SOC.lower()}_band', fontsize=14)
        # for data in Nd_Data:
        for i in range(len(Nd_Data)):
            py = Nd_Data[i]
            px = np.linspace(0, 1, len(py))
            repy = re_Nd_Data[i]
            if i == 0:
                ax.plot(px, py, 
                        color='#738a79', 
                        linestyle='--',
                        label=f'{r}')
                
                ax.plot(px, repy,
                        color='#b57fb3', 
                        linestyle=':',
                        label=f'{re_r}'
                        )
            else:
                ax.plot(px, py, 
                        color='#738a79', 
                        linestyle='--')
                ax.plot(px, repy,
                        color='#b57fb3', 
                        linestyle=':'
                        )
        ax.legend(loc='best')
        
        ax.set_xticks(ticks=tick, labels=label)
        # ax.set_ylim(ymin=-2, ymax=2)
        ax.set_ylim(ymin=-10, ymax=10)
        ax.grid('--')

        plt.tight_layout()
        
    def VpVpPlot(self, tick, label):

        PLBasis = self.PLBasis
        SOC     = self.SOC
        re_SOC  = self.re_SOC
        ndBasis = self.ndBasis
        atom    = self.atom

        re_Vp_Data = np.loadtxt(f'{self.re_vp_path}')
        Vp_Data = np.loadtxt(f'{self.vp_path}')

        index = 0
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.set_title(f'Vasp Band\nsoc vs. ncl', fontsize=14)
        
        for i in range(1, len(Vp_Data[0])):
            px = Vp_Data[:,0]
            px = px/np.max(px)
            py = Vp_Data[:,i]
            if i == 1:
                ax.plot(px, py, 
                        color='#738a79', 
                        linestyle='--',
                        label=f'{SOC}')
            else:
                ax.plot(px, py, 
                        color='#738a79', 
                        linestyle='--',)
                
        for i in range(1, len(re_Vp_Data[0])):
            px = re_Vp_Data[:,0]
            px = px/np.max(px)
            repy = re_Vp_Data[:,i]
            if i == 1:
                ax.plot(px, repy, 
                        color='#b57fb3', 
                        linestyle=':',
                        label=f'{re_SOC}')
            else:
                ax.plot(px, repy, 
                        color='#b57fb3', 
                        linestyle=':',)
                
        ax.legend(loc='best')
        ax.set_xticks(ticks=tick, labels=label)
        ax.set_ylim(ymin=-10, ymax=10)
        ax.grid('--')

        plt.tight_layout()

    def VpNdNdPlot(self, tick, label):
        atom    = self.atom
        PLBasis = self.PLBasis
        SOC     = self.SOC
        re_SOC  = self.re_SOC
        ndBasis = self.ndBasis

        Nd_Data = np.loadtxt(f'{self.nd_path}')
        re_Nd_Data = np.loadtxt(f'{self.re_nd_path}')
        Vp_Data = np.loadtxt(f'{self.vp_path}')
        # Vp_Data = np.loadtxt(f'{self.re_vp_path}')

        fig = plt.figure()
        ax = plt.subplot(111)
        # ax.set_title(f'Nanodcal ncl, soc and Vasp {SOC}', fontsize=14)
        # for data in Nd_Data:
        for i in range(len(Nd_Data)):
            py = Nd_Data[i]
            px = np.linspace(0, 1, len(py))
            repy = re_Nd_Data[i]
            
            if i == 0:
                ax.plot(px, py, 
                        color='#738a79', 
                        linestyle='-',
                        linewidth=2,
                        label=f'NanoDCAL_TM-{PLBasis}_{SOC}')
                
                ax.plot(px, repy,
                        color='#b57fb3', 
                        linestyle='--',
                        linewidth=1.5,
                        label=f'NanoDCAL_TM-{PLBasis}_{re_SOC}'
                        )
                
                
            else:
                ax.plot(px, py, 
                        color='#738a79', 
                        linestyle='-',
                        linewidth=2,)
                ax.plot(px, repy,
                        color='#b57fb3', 
                        linestyle='--',
                        linewidth=1.5,
                        )

        for i in range(1, len(Vp_Data[0])):
            px = Vp_Data[:,0]
            px = px/np.max(px)
            py = Vp_Data[:,i]
            if i == 1:
                ax.plot(px, py,
                    color='#d28063', 
                    linestyle=':',
                    label=f'Vasp_{PLBasis}_{SOC}',
                    linewidth=1.5
                    )
            else:
                ax.plot(px, py,
                    color='#d28063', 
                    linestyle=':',
                    linewidth=1.5
                    )
                    
        ax.legend(loc='upper right')
        ax.set_xticks(ticks=tick, labels=label)
        ax.set_ylim(ymin=-10, ymax=10)
        ax.grid('--')
        # ax.hlines(y=-9, xmin=0, xmax=1)
        plt.tight_layout()

