import numpy as np
import matplotlib.pyplot as plt
import itertools
import os
from scipy.optimize import curve_fit

from STT_Tool import Tools as stt
from STT_Tool import Tools
color = Tools.ColorList()

class SingleR_TB_Plot():
    def __init__(self, DataPath:str, Rtot:int):

        self.DataPath = DataPath
        self.Rtot = Rtot

    def _LoadData(self, DataPath):
        
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

        self.data_caseList = data_caseList
        self.data_RList    = data_RList

        return data_caseList, data_RList

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

    def _NormalDisMask(self, stdNum=3):
        def GaussionFit(x0, y0, fitt):
            def fit_function(x, A, mu, sigma):
                factor = 1/np.sqrt(2*np.pi*sigma**2)
                main   = np.e**(-((x-mu)**2)/2*sigma**2)
                return A*factor*main
            popt, pcov = curve_fit(fit_function, x0, y0, p0=fitt)
            new_x = np.linspace(np.min(x0), np.max(x0), 1000)
            new_y = fit_function(new_x, *popt)
            # print("A={}, mu={}, sigma={}".format(*popt))
            return new_x, new_y, popt
 
        Ntot = self.Rtot
        Ntot_Num_array = np.array([])
        for NLL in range(1, Ntot+1, 1):
            NL_set  = stt.Sum_Combination(NLL, Ntot)
            totNum = 0
            for NL in NL_set:
                sub_NL_set = list(set(itertools.permutations(NL)))
                totNum = totNum+len(sub_NL_set)
            Ntot_Num = np.array([NLL, totNum])
            Ntot_Num_array = np.append(Ntot_Num_array, Ntot_Num)
        Ntot_Num_array = Ntot_Num_array.reshape(Ntot, 2)
        # print(Ntot_Num_array)
        #  ----------
        x0 = Ntot_Num_array[:,0]
        y0 = Ntot_Num_array[:,1]/np.max(Ntot_Num_array[:,1])
        new_x, new_y, popt = GaussionFit(x0, y0, fitt=(2, 7, 0.5))
        INstdpos = np.where(np.logical_and(
                            x0<=popt[1]+stdNum*popt[2], x0>=popt[1]-stdNum*popt[2])
                            )
        return INstdpos[0]

    def _LabelSet(self, rescale):
        match rescale:
            case 1e6:
                plt.ylabel(r"Resistance (M$\Omega$)", fontsize=14)
            
            case 1e3:
                plt.ylabel(r"Resistance (k$\Omega$)", fontsize=14)
            case 1:
                plt.ylabel(r"Resistance ($\Omega$)", fontsize=14)
            case _:
                plt.ylabel("NO Support This Rescale", fontsize=18)

    def casePlot(self, rescale=1e6, maxminBound=(-1e10, 1e10), Gaussian=False, p=False):
        
        maximum = maxminBound[0]
        minimum = maxminBound[1]
        scale   = rescale

        data_caseList, data_RList = self._LoadData(self.DataPath)
        Rtot          = self.Rtot
        
        INstdpos = self._NormalDisMask()
        print("These number of R are in the Gaussian Distribution: " + str((INstdpos+1)[::-1]))
        # ==========

        plt.figure(figsize=(8, 4))
        plt.rcParams['ytick.direction']='in'

        xpos = 0
        x_xpos = -0.5
        xpos_list = []
        for i, case in enumerate(data_caseList[::-1]):
            for subcase in case:
                x = [xpos for _ in range(len(subcase))]
                if Gaussian:
                    if i in INstdpos:
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
                    else:
                        pass
                else:
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
        self._LabelSet(rescale=rescale)
        plt.grid("--", axis='y')

        if p:
            plt.savefig("case_Result.png",
                        transparent=True
                        )
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

    def RPlot(self, Gaussian=False, p=False):
        
        maximum = self.maxminBound[0]
        minimum = self.maxminBound[1]
        scale   = self.rescale

        data_caseList, data_RList = self._LoadData(self.DataPath)
        Rtot       = self.Rtot
        
        INstdpos = self._NormalDisMask()
        # ==========

        plt.figure()
        plt.rcParams['ytick.direction']='in'
        xpos=0
        xpos_list = []
        ax = plt.subplot(111)
        for i, subR in enumerate(data_RList[::-1]):
            x = [xpos for _ in range(len(subR))]
            if Gaussian:
                if i in INstdpos:
                    ax.scatter(x, np.array(subR)/scale, marker='*', 
                                c='k', 
                                alpha=0.5
                                )
                    ax.errorbar(xpos, np.mean(subR)/scale,
                                    fmt='o', color="blue",
                                    yerr=np.std(subR)/scale,
                                    ecolor='red',
                                    capsize=4, elinewidth=2
                                    )
                else:
                    pass
            else:
                if p=="More":
                    pass
                else:
                    ax.errorbar(xpos, np.mean(subR)/scale,
                                fmt='o', color="blue",
                                yerr=np.std(subR)/scale,
                                ecolor='red',
                                capsize=4, elinewidth=2
                                )
                ax.scatter(x, np.array(subR)/scale, marker='*', 
                            c='k', 
                            alpha=0.5
                            )
                

            xpos_list.append(xpos)
            xpos = xpos+1
            maximum, minimum = self._BoundarySet(data=subR, Bound=(maximum, minimum))
            
            ax.vlines(x=xpos-0.5, 
                       ymin=0, ymax=(maximum+scale*10)/scale, 
                       color=color[-1]
                       )


        plt.ylim(
            (minimum-0.5*10**(6))/scale, (maximum+0.5*10**(6))/scale
            )

        xlabel = [f"R{i}" for i in range(Rtot, 0, -1)]
        plt.xticks(ticks=xpos_list,
                labels=xlabel)
        # plt.xlim(xpos_list[np.min(INstdpos)]-0.5, xpos_list[np.max(INstdpos)]-0.5)
        plt.xlim(-0.5, xpos-0.5)

        self._LabelSet(rescale=scale)

        plt.grid("--", axis='y')

        if p=="True":
            plt.savefig("case_R_Result.png",
                        transparent=True
                        )
        elif p=="More":
            pass

        else:
            plt.show()
        # =====
        LABELoutput  = (xlabel, xpos_list)
        RESULToutput = data_RList
        # =====
        return ax

class MultiR_TB_Plot_Angle():
    def __init__(self, Rtot:int, DataPath="Temperature_Rparse"):
        
        path = os.listdir(f"{DataPath}")

        targetPath = []
        for p in path:
            if p[-4:] == ".dat":
                targetPath.append(p)
        
        self.DataPath   = DataPath
        self.targetPath = targetPath
        self.Rtot       = Rtot
        print(targetPath)
    
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

    def _NormalDisMask(self, stdNum=3):
        def GaussionFit(x0, y0, fitt):
            def fit_function(x, A, mu, sigma):
                factor = 1/np.sqrt(2*np.pi*sigma**2)
                main   = np.e**(-((x-mu)**2)/2*sigma**2)
                return A*factor*main
            popt, pcov = curve_fit(fit_function, x0, y0, p0=fitt)
            new_x = np.linspace(np.min(x0), np.max(x0), 1000)
            new_y = fit_function(new_x, *popt)
            print("A={}, mu={}, sigma={}".format(*popt))
            return new_x, new_y, popt
        
        Ntot = self.Rtot
        Ntot_Num_array = np.array([])
        for NLL in range(1, Ntot+1, 1):
            NL_set  = stt.Sum_Combination(NLL, Ntot)
            totNum = 0
            for NL in NL_set:
                sub_NL_set = list(set(itertools.permutations(NL)))
                totNum = totNum+len(sub_NL_set)
            Ntot_Num = np.array([NLL, totNum])
            Ntot_Num_array = np.append(Ntot_Num_array, Ntot_Num)
        Ntot_Num_array = Ntot_Num_array.reshape(Ntot, 2)
        # print(Ntot_Num_array)
        #  ----------
        x0 = Ntot_Num_array[:,0]
        y0 = Ntot_Num_array[:,1]/np.max(Ntot_Num_array[:,1])
        new_x, new_y, popt = GaussionFit(x0, y0, fitt=(2, 7, 0.5))
        INstdpos = np.where(np.logical_and(
                            x0<=popt[1]+stdNum*popt[2], x0>=popt[1]-stdNum*popt[2])
                            )
        return INstdpos

    def _LoadData(self, DataPath):
        
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

    def _casePlot(self, ax, case_data, 
                      color, label=False,
                      scale=1e6, maxminBound=(-1e10, 1e10)):
        
        data_caseList = case_data
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

                maximum, minimum = self._BoundarySet(data=subcase, Bound=(maximum, minimum))
            
            xpos_list.append((x_xpos+(xpos-0.5))/2)
            ax.vlines(x=xpos-0.5, 
                       ymin=(minimum-scale*10)/scale, ymax=(maximum+scale*10)/scale, 
                       color='k', alpha=0.5
                       )
            x_xpos = xpos-0.5

        return maximum, minimum, xpos, xpos_list

    def _LabelSet(self, ax, rescale):
        match rescale:
            case 1e6:
                ax.set_ylabel(r"Resistance (M$\Omega$)", fontsize=14)
            
            case 1e3:
                ax.set_ylabel(r"Resistance (k$\Omega$)", fontsize=14)
            case 1:
                ax.set_ylabel(r"Resistance ($\Omega$)", fontsize=14)
            case _:
                ax.set_ylabel("NO Support This Rescale", fontsize=18)

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
            data_caseList, data_RList = self._LoadData(DataPath=f"{self.DataPath}/{datapath}")

            maximum, minimum, xpos, xpos_list = self._casePlot(ax, data_caseList, 
                                                            color=(color[i], color[i]), 
                                                            label=f"{targetPath[i][:-4]}",
                                                            scale=rescale, 
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
        self._LabelSet(ax, rescale=rescale)
        plt.grid("--", axis='y')
        plt.legend(loc='best')

        if p:
            plt.savefig("diffAngle_Result.png", 
                        transparent=True
                        )
        else:
            plt.show()

class MultiR_TB_Plot_Rtot():
    def __init__(self, RtotList:tuple, AngleSet=(90, 180, 90), DataPath="Temperature_Rparse"):
        
        self.DataPath = DataPath
        self.RtotList = RtotList
        self.AngleSet = AngleSet

        AngleDict = {}
        for N in RtotList:
            targetPathList = []
            path = os.listdir(f"{DataPath}/{N}_Ntot")
            for p in path:
                if p[-4:] == ".dat":
                    targetPathList.append(p)
            AngleDict[f"{N}_Ntot"]=targetPathList
        self.AngleDict = AngleDict

    def _LoadData(self, Rtot, DataPath):

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

    def _RPlot(self, ax, R_data, start,
                    color, label,
                    scale=1e6, maxminBound=(-1e10, 1e10)):
    
        data_RList = R_data
        maximum  = maxminBound[0]
        minimum  = maxminBound[1]
        scale    = scale
        RtotList = self.RtotList

        # =====

        xpos = (0+start)
        xpos_list = []
        for i, subR in enumerate(data_RList[::-1]):
            x = [(xpos+start) for _ in range(len(subR))]
            if i==0:
                plt.errorbar(xpos, np.mean(subR)/scale,
                                fmt='o', color=color,
                                yerr=np.std(subR)/scale,
                                ecolor=color,
                                capsize=4, elinewidth=2,
                                label=label
                                )
            else:
                plt.errorbar(xpos, np.mean(subR)/scale,
                                fmt='o', color=color,
                                yerr=np.std(subR)/scale,
                                ecolor=color,
                                capsize=4, elinewidth=2,
                                )


            xpos_list.append(xpos)
            xpos = xpos+1
        plt.legend(loc='best')

        return maximum, minimum, xpos, xpos_list
    
    def _LabelSet(self, ax, rescale):
        match rescale:
            case 1e6:
                ax.set_ylabel(r"Resistance (M$\Omega$)", fontsize=14)
            
            case 1e3:
                ax.set_ylabel(r"Resistance (k$\Omega$)", fontsize=14)
            case 1:
                ax.set_ylabel(r"Resistance ($\Omega$)", fontsize=14)
            case _:
                ax.set_ylabel("NO Support This Rescale", fontsize=18)

    def _NormalDisMask(self, Ntot, stdNum=3):
        def GaussionFit(x0, y0, fitt):
            def fit_function(x, A, mu, sigma):
                factor = 1/np.sqrt(2*np.pi*sigma**2)
                main   = np.e**(-((x-mu)**2)/2*sigma**2)
                return A*factor*main
            popt, pcov = curve_fit(fit_function, x0, y0, p0=fitt)
            new_x = np.linspace(np.min(x0), np.max(x0), 1000)
            new_y = fit_function(new_x, *popt)
            # print("A={}, mu={}, sigma={}".format(*popt))
            return new_x, new_y, popt

        Ntot_Num_array = np.array([])
        for NLL in range(1, Ntot+1, 1):
            NL_set  = stt.Sum_Combination(NLL, Ntot)
            totNum = 0
            for NL in NL_set:
                sub_NL_set = list(set(itertools.permutations(NL)))
                totNum = totNum+len(sub_NL_set)
            Ntot_Num = np.array([NLL, totNum])
            Ntot_Num_array = np.append(Ntot_Num_array, Ntot_Num)
        Ntot_Num_array = Ntot_Num_array.reshape(Ntot, 2)
        # print(Ntot_Num_array)
        #  ----------
        x0 = Ntot_Num_array[:,0]
        y0 = Ntot_Num_array[:,1]/np.max(Ntot_Num_array[:,1])
        new_x, new_y, popt = GaussionFit(x0, y0, fitt=(2, 7, 0.5))
        INstdpos = np.where(np.logical_and(
                            x0<=popt[1]+stdNum*popt[2], x0>=popt[1]-stdNum*popt[2])
                            )
        return INstdpos[0]

    def RPlot(self, rescale=1e6, maxminBound=(-1e10, 1e10), p=False):

        maximum = maxminBound[0]
        minimum = maxminBound[1]
        scale   = rescale

        RtotList = self.RtotList
        AngleSet = self.AngleSet
        DataPath = self.DataPath
        
        # ==========

        plt.figure(figsize=(8, 6))
        plt.rcParams['ytick.direction']='in'
        ax = plt.subplot(111)
        
        for i, Rtot in enumerate(RtotList):
            if len(AngleSet)==3:
                data_caseList, data_RList = self._LoadData(Rtot=Rtot, DataPath=f"{DataPath}/{Rtot}_Ntot/{AngleSet[0]}_{AngleSet[1]}_{AngleSet[2]}.dat")
            else:
                data_caseList, data_RList = self._LoadData(Rtot=Rtot, DataPath=f"{DataPath}/{Rtot}_Ntot/{AngleSet[i][0]}_{AngleSet[i][1]}_{AngleSet[i][2]}.dat")

            maximum, minimum, xposR, xpos_listR = self._RPlot(ax, data_RList, start=i,
                                                              color=color[i], label=f"Rtot={Rtot}",
                                                              scale=rescale, 
                                                              maxminBound=(maximum, minimum)
                                                              )
            
            if i==0:
                xpos, xpos_list = xposR, xpos_listR
            else:
                pass
            

        xlabel = [f"R{i}" for i in range(np.max(RtotList), 0, -1)]
        plt.xticks(ticks=xpos_list,
                    labels=xlabel)

        plt.xlim(-0.5, xpos-0.5)

        self._LabelSet(ax, rescale=rescale)

        plt.grid("--", axis='y')
        plt.show()

    def RtotPlot(self, rescale=1e6, maxminBound=(-1e10, 1e10), p=False):
        maximum = maxminBound[0]
        minimum = maxminBound[1]
        scale   = rescale

        RtotList = self.RtotList
        AngleSet = self.AngleSet
        DataPath = self.DataPath
        
        # ==========

        plt.figure(figsize=(8, 6))
        plt.rcParams['ytick.direction']='in'
        ax = plt.subplot(111)
        
        xlabel = []
        for i, Rtot in enumerate(RtotList):
            if len(AngleSet)==3:
                data_caseList, data_RList_sep = self._LoadData(Rtot=Rtot, DataPath=f"{DataPath}/{Rtot}_Ntot/{AngleSet[0]}_{AngleSet[1]}_{AngleSet[2]}.dat")
            else:
                data_caseList, data_RList_sep = self._LoadData(Rtot=Rtot, DataPath=f"{DataPath}/{Rtot}_Ntot/{AngleSet[i][0]}_{AngleSet[i][1]}_{AngleSet[i][2]}.dat")
            
            data_RList = []
            for data in data_RList_sep:
                data_RList = data_RList+data
            data_RList = [data_RList]
            if len(AngleSet)==3:
                maximum, minimum, xposR, xpos_listR = self._RPlot(ax, data_RList, start=i,
                                                                color=color[i], label=f"Rtot={Rtot}: {AngleSet[0]}_{AngleSet[1]}_{AngleSet[2]}",
                                                                scale=rescale, 
                                                                maxminBound=(maximum, minimum)
                                                                )
            else:
                maximum, minimum, xposR, xpos_listR = self._RPlot(ax, data_RList, start=i,
                                                                color=color[i], label=f"Rtot={Rtot}: {AngleSet[i][0]}_{AngleSet[i][1]}_{AngleSet[i][2]}",
                                                                scale=rescale, 
                                                                maxminBound=(maximum, minimum)
                                                                )

            xlabel.append(f"{Rtot}Rtot")
        xpos_list = np.linspace(0, len(xlabel)-1, len(xlabel))

        plt.xticks(ticks=xpos_list,
                    labels=xlabel)
        self._LabelSet(ax, rescale=rescale)
        plt.grid("--", axis='y')
        plt.show()

    def Rtot_AnglePlot(self, rescale=1e6, maxminBound=(-1e10, 1e10), p=False, colorB="Case"):
        maximum = maxminBound[0]
        minimum = maxminBound[1]
        scale   = rescale

        RtotList = self.RtotList
        DataPath = self.DataPath
        AngleSet = self.AngleDict
        
        # ==========

        plt.figure(figsize=(8, 6))
        plt.rcParams['ytick.direction']='in'
        ax = plt.subplot(111)
        
        xlabel = []
        xpos_list = []
        pos = 0
        for j, Rtot in enumerate(RtotList):
            AngleDat = AngleSet[f"{Rtot}_Ntot"]
            for i, angle in enumerate(AngleDat):
                data_caseList, data_RList_sep = self._LoadData(Rtot=Rtot, DataPath=f"{DataPath}/{Rtot}_Ntot/{angle}")
            
                data_RList = []
                for data in data_RList_sep:
                    data_RList = data_RList+data
                data_RList = [data_RList]
                if colorB == "Case":
                    c = color[j]
                    l = None
                elif colorB == "Angle":
                    c = color[i]
                    l = fr"$\gamma$L={angle[0:2]}"
            
                if j == 0:
                    maximum, minimum, xposR, xpos_listR = self._RPlot(ax, data_RList, start=pos,
                                                                    color=c, 
                                                                    label=l,
                                                                    scale=rescale, 
                                                                    maxminBound=(maximum, minimum)
                                                                    )
                else:
                    maximum, minimum, xposR, xpos_listR = self._RPlot(ax, data_RList, start=pos,
                                                                    color=c, 
                                                                    label=None,
                                                                    scale=rescale, 
                                                                    maxminBound=(maximum, minimum)
                                                                    )
                if i==len(AngleDat)-1:
                    pass
                else:
                    pos+=1

            xlabel.append(f"{Rtot}Rtot")
            xpos_list.append(pos)

        plt.xticks(ticks=xpos_list,
                   labels=xlabel)              
        self._LabelSet(ax, rescale=rescale)
        plt.grid("--")
        # plt.show()