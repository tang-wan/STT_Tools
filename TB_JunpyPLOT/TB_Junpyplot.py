import numpy as np
import matplotlib.pyplot as plt
import os

from STT_Tool import Tools
color = Tools.ColorList()

class SingleR_TB_Plot():
    def __init__(self, DataPath:str, Rtot:int):
        
        with open(f"{DataPath}") as R_File:
            data = R_File.readlines()

        data_caseList = [[] for _ in range(Rtot)]
        data_RList = [[] for _ in range(Rtot)]
        for case in data:
            case = case.split()
            del case[1]
            Rgroup = int(case[0][1:])-1
            
            for i in case[1:]:
                data_RList[Rgroup].append(float(i))

            data_caseList[Rgroup].append(np.float64(case[1:]))
        
        self.data_caseList = data_caseList
        self.data_RList    = data_RList
        self.Rtot          = Rtot

    def _BoundarySet(self, data, Bound:tuple):
        maximum = Bound[0]
        minimum = Bound[1]
        
        if np.max(data)>maximum:
            maximum = np.max(data)
        else:
            pass
        if np.min(data)<minimum:
            minimum = np.min(data)
        else:
            pass

        return maximum, minimum    
        
    def casePlot(self, rescale=1e6, maxminBound=(-1e10, 1e10), p=False):
        
        maximum = maxminBound[0]
        minimum = maxminBound[1]
        scale   = rescale

        data_caseList = self.data_caseList
        Rtot          = self.Rtot
        
        # ==========

        plt.figure(figsize=(8, 4))
        plt.rcParams['ytick.direction']='in'

        xpos = 0
        x_xpos = -0.5
        xpos_list = []
        for case in data_caseList[::-1]:
            for subcase in case:
                x = [xpos for _ in range(len(subcase))]
                plt.scatter(x, subcase/scale, marker='*',
                            c='k',
                            alpha=0.5
                            )
                plt.errorbar(xpos, np.mean(subcase)/scale,
                            fmt='o', color="blue",
                            yerr=np.std(subcase)/scale,
                            ecolor='red',
                            capsize=4, elinewidth=2
                            )
                xpos = xpos+1

                maximum, minimum = self._BoundarySet(data=subcase, Bound=(maximum, minimum))
            
            xpos_list.append((x_xpos+(xpos-0.5))/2)
            plt.vlines(x=xpos-0.5, 
                       ymin=(minimum-scale*10)/scale, ymax=(maximum+scale*10)/scale, 
                       color=color[-1]
                       )
            x_xpos = xpos-0.5
    
        xlabel = [f"R{i}" for i in range(Rtot, 0, -1)]
        plt.xlim(-0.5, xpos-0.5)
        plt.xticks(ticks=xpos_list,
                   labels=xlabel
                   )
        plt.ylim(
            (minimum-0.5*scale)/scale, (maximum+0.5*scale)/scale
            )
        plt.ylabel("Resistance (MΩ)")
        plt.grid("--", axis='y')

        if p:
            plt.savefig("Rcase_Result.png")
        else:
            plt.show()
        # =====
        LABELoutput  = (xlabel, xpos_list)
        RESULToutput = data_caseList
        # =====
        self.rescale     = rescale
        self.maxminBound = maxminBound
        # =====
        return LABELoutput, RESULToutput

    def RPlot(self, p=False):
        
        maximum = self.maxminBound[0]
        minimum = self.maxminBound[1]
        scale   = self.rescale

        data_RList = self.data_RList
        Rtot       = self.Rtot
        
        # ==========

        plt.figure()
        plt.rcParams['ytick.direction']='in'
        xpos=0
        xpos_list = []
        for subR in data_RList[::-1]:
            x = [xpos for _ in range(len(subR))]
            plt.scatter(x, np.array(subR)/scale, marker='*', 
                        c='k', 
                        alpha=0.5
                        )
            plt.errorbar(xpos, np.mean(subR)/scale,
                            fmt='o', color="blue",
                            yerr=np.std(subR)/scale,
                            ecolor='red',
                            capsize=4, elinewidth=2
                            )
            xpos_list.append(xpos)
            xpos = xpos+1
            maximum, minimum = self._BoundarySet(data=subR, Bound=(maximum, minimum))
            
            plt.vlines(x=xpos-0.5, 
                       ymin=(minimum-scale*10)/scale, ymax=(maximum+scale*10)/scale, 
                       color=color[-1]
                       )


        plt.ylim(
            (minimum-0.5*10**(6))/scale, (maximum+0.5*10**(6))/scale
            )

        xlabel = [f"R{i}" for i in range(Rtot, 0, -1)]
        plt.xlim(-0.5, xpos-0.5)
        plt.xticks(ticks=xpos_list,
                labels=xlabel)

        plt.ylabel("Resistance (MΩ)")

        plt.grid("--", axis='y')
        plt.show()

        if p:
            plt.savefig("Rcase_Result.png")
        else:
            plt.show()
        # =====
        LABELoutput  = (xlabel, xpos_list)
        RESULToutput = data_RList
        # =====

class MultiR_TB_Plot():
    def __init__(self, Rtot, DataPath="Temperature_Rparse"):
        
        path = os.listdir(f"{DataPath}")

        targetPath = []
        for p in path:
            if p[-4:] == ".dat":
                targetPath.append(p)
        
        self.targetPath = targetPath
        self.Rtot       = Rtot
        print(targetPath)
    
    def Tool_BoundarySet(self, data, Bound:tuple):
        maximum = Bound[0]
        minimum = Bound[1]
        
        if np.max(data)>maximum:
            maximum = np.max(data)
        else:
            pass
        if np.min(data)<minimum:
            minimum = np.min(data)
        else:
            pass

        return maximum, minimum       

    def Tool_LoadData(self, DataPath):
        
        Rtot = self.Rtot

        with open(f"{DataPath}") as R_File:
            data = R_File.readlines()

        data_caseList = [[] for _ in range(Rtot)]
        data_RList = [[] for _ in range(Rtot)]
        for case in data:
            case = case.split()
            del case[1]
            Rgroup = int(case[0][1:])-1
            
            for i in case[1:]:
                data_RList[Rgroup].append(float(i))

            data_caseList[Rgroup].append(np.float64(case[1:]))
        
        return data_caseList, data_RList

    def Tool_casePlot(self, ax, case_data, 
                      color, label=False,
                      scale=1e6, maxminBound=(-1e10, 1e10)):
        
        data_caseList = case_data
        targetPath = self.targetPath
        maximum = maxminBound[0]
        minimum = maxminBound[1]

        xpos = 0
        x_xpos = -0.5
        xpos_list = []
        for k, case in enumerate(data_caseList[::-1]):
            for subcase in case:
                x = [xpos for _ in range(len(subcase))]
                if label and k==0:
                    ax.errorbar(xpos, np.mean(subcase)/scale,
                                fmt='o', color=color[0],
                                yerr=np.std(subcase)/scale,
                                ecolor=color[1],
                                capsize=4, elinewidth=2,
                                label=label
                                )
                else:
                    ax.errorbar(xpos, np.mean(subcase)/scale,
                                fmt='o', color=color[0],
                                yerr=np.std(subcase)/scale,
                                ecolor=color[1],
                                capsize=4, elinewidth=2,
                                )
                xpos = xpos+1

                maximum, minimum = self.Tool_BoundarySet(data=subcase, Bound=(maximum, minimum))
            
            xpos_list.append((x_xpos+(xpos-0.5))/2)
            ax.vlines(x=xpos-0.5, 
                       ymin=(minimum-scale*10)/scale, ymax=(maximum+scale*10)/scale, 
                       color='k', alpha=0.5
                       )
            x_xpos = xpos-0.5

        return maximum, minimum, xpos, xpos_list

    def casePlot(self, rescale=1e6, maxminBound=(-1e10, 1e10), p=False):
        maximum = maxminBound[0]
        minimum = maxminBound[1]
        scale   = rescale

        Rtot          = self.Rtot
        targetPath    = self.targetPath
        
        # ==========

        plt.figure(figsize=(8, 6))
        plt.rcParams['ytick.direction']='in'
        ax = plt.subplot(111)

        for i, datapath in enumerate(targetPath):
            data_caseList, data_RList = self.Tool_LoadData(DataPath=f"Temperature_Rparse/{datapath}")

            maximum, minimum, xpos, xpos_list = self.Tool_casePlot(ax, data_caseList, 
                                                                   color=(color[i], color[i]), 
                                                                   label=f"{targetPath[i][:-4]}",
                                                                   scale=1e6, 
                                                                   maxminBound=(maximum, minimum)
                                                                   )
    
        xlabel = [f"R{i}" for i in range(Rtot, 0, -1)]
        plt.xlim(-0.5, xpos-0.5)
        plt.xticks(ticks  = xpos_list,
                   labels = xlabel
                   )
        plt.ylim(
            (minimum-0.5*scale)/scale, (maximum+0.5*scale)/scale
            )
        plt.ylabel("Resistance (MΩ)")
        plt.grid("--", axis='y')
        plt.legend(loc='best')

        if p:
            plt.savefig("Rcase_Result.png")
        else:
            plt.show()



        
        
    