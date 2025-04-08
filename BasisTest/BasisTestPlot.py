import numpy as np
import matplotlib.pyplot as plt

from STT_Tool.nanodcalBAND import nanodcalBANDplot as nBp
from STT_Tool.vaspBAND import vaspBANDplot_v2 as vBp
from STT_Tool import Tools
color = Tools.ColorList()

class BasisTest():
    def __init__(self, ndPath_soc, 
                       ndPath_ncl, 
                       vpPath_soc, vpFermi_soc, 
                       vpPath_ncl, vpFermi_ncl,
                    #    newPath_soc,
                    #    newPath_ncl,
                       ):

        def ReadNdData(path):
            ReadData = nBp.nanodcalBANDplot(path)
            BANDoutput, LABELoutput = ReadData.Read_AllData()
            return BANDoutput, LABELoutput
        
        def ReadVpData(path, EF):
            ReadData = vBp.procarBNADplot(FilePath=path, fermiEnergy=EF, PROCARtype='vasp')
            BANDoutput, LABELoutput, kwargs_plot1   = ReadData.Read_AllData_Band()
            return BANDoutput, LABELoutput

        self.ndBandData_soc, self.ndLabelData_soc = ReadNdData(ndPath_soc)
        self.ndBandData_ncl, self.ndLabelData_ncl = ReadNdData(ndPath_ncl)

        self.vpBandData_soc, self.vpLabelData_soc = ReadVpData(vpPath_soc, vpFermi_soc)
        self.vpBandData_ncl, self.vpLabelData_ncl = ReadVpData(vpPath_ncl, vpFermi_ncl)

        # if newPath_soc:
        #     self.ndBandData_soc_new, self.ndLabelData_soc_new = ReadNdData(newPath_soc)
        
        # if newPath_ncl:
        #     self.ndBandData_ncl_new, self.ndLabelData_ncl_new = ReadNdData(newPath_ncl)
    
    def _NormalizeTool(self, data):
        data = data-np.min(data)
        data = data/np.max(data)
        return data

    def _PlotTool(self, ax, kpath, band, lineinfo:tuple):
        # ax = plt.subplot(111)
        for i, b in enumerate(range(len(band))):
            if i == 0:
                ax.plot(kpath, band[b], 
                        color=lineinfo[2], 
                        linestyle=lineinfo[0], linewidth=1,
                        label=lineinfo[1]
                        )
            else:
                ax.plot(kpath, band[b], 
                        color=lineinfo[2], 
                        linestyle=lineinfo[0], linewidth=1
                        )

    def SingleBand(self, CalType, socType, boundary=(-10, 15), s=False):
        match CalType:
            case "vasp":
                match socType:
                    case "soc":
                        data1_Kpath, data1_Band = self.vpBandData_soc
                        data1_Label, data1_Tick = self.vpLabelData_soc
                    case "ncl":
                        data1_Kpath, data1_Band = self.vpBandData_ncl
                        data1_Label, data1_Tick = self.vpLabelData_ncl
                    case _:
                        print("No this kind of basis type")
            # ===============
            case "nanodcal":
                match socType:
                    case "soc":
                        data1_Kpath, data1_Band = self.ndBandData_soc
                        data1_Label, data1_Tick = self.ndLabelData_soc
                    case "ncl":
                        data1_Kpath, data1_Band = self.ndBandData_ncl
                        data1_Label, data1_Tick = self.ndLabelData_ncl
                    case _:
                        print("No this kind of basis type")
            case _:
                print("No this kind of calculation type")
        
        lineStyle1 = '-'
        caseType1  = f"{CalType}_{socType}"

        data_Kpath = self._NormalizeTool(data1_Kpath)
        data_Tick  = self._NormalizeTool(data1_Tick)
        # data_Kpath = data1_Kpath
        # data_Tick  = data1_Tick

        title = f"{CalType} with {socType}"
        # >>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ax = plt.subplot(111)
        self._PlotTool(ax, kpath=data_Kpath, band=data1_Band, 
                       lineinfo=(lineStyle1, caseType1, color[-1]))
        ax.set_ybound(lower=boundary[0], upper=boundary[1])
        ax.set_xticks(ticks=data_Tick,
                    labels=data1_Label)
        ax.set_title(title)
        ax.grid('--')
        if s:
            plt.savefig(f"SingleBand_{CalType}_with_{socType}.png", transparent=True, bbox_inches='tight')
        else:
            plt.show()

    def CompBand(self, CompType, title, boundary=(-10, 15), s=False):
        match CompType:
            case "nanodcal":
                caseType1 = "soc"
                data1_Kpath, data1_Band = self.ndBandData_soc
                caseType2 = "ncl"
                data2_Kpath, data2_Band = self.ndBandData_ncl
            case "vasp":
                caseType1 = "soc"
                data1_Kpath, data1_Band = self.vpBandData_soc
                caseType2 = "ncl"
                data2_Kpath, data2_Band = self.vpBandData_ncl
            case "soc":
                caseType1 = "vasp"
                data1_Kpath, data1_Band = self.vpBandData_soc
                caseType2 = "nanodcal"
                data2_Kpath, data2_Band = self.ndBandData_soc
            case "ncl":
                caseType1 = "vasp"
                data1_Kpath, data1_Band = self.vpBandData_ncl
                caseType2 = "nanodcal"
                data2_Kpath, data2_Band = self.ndBandData_ncl
        
        data1_Kpath = self._NormalizeTool(data1_Kpath)
        lineStyle1  = '-'

        data2_Kpath = self._NormalizeTool(data2_Kpath)
        lineStyle2  = '--'

        xlabel, xtick = self.vpLabelData_soc
        xtick = self._NormalizeTool(xtick)

        # >>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ax = plt.subplot(111)
        self._PlotTool(ax, kpath=data1_Kpath, band=data1_Band, 
                       lineinfo=(lineStyle1, caseType1, color[-1]))
        self._PlotTool(ax, kpath=data2_Kpath, band=data2_Band, 
                       lineinfo=(lineStyle2, caseType2, color[1]))
        
        ax.set_ybound(lower=boundary[0], upper=boundary[1])
        ax.set_xticks(ticks=xtick,
                      labels=xlabel)
        if title:
            ax.set_title(title)
        ax.grid('--')
        plt.legend(loc='upper right')
        if s:
            plt.savefig(f"CompareBand_{CompType}.png")
        else:
            plt.show()

    def CompBand_selfLabel(self, CompType, title:str, label:tuple, boundary=(-10, 15), s=False):
        match CompType:
            case "nanodcal":
                data1_Kpath, data1_Band = self.ndBandData_soc
                data2_Kpath, data2_Band = self.ndBandData_ncl
            case "vasp":
                data1_Kpath, data1_Band = self.vpBandData_soc
                data2_Kpath, data2_Band = self.vpBandData_ncl
            case "soc":
                data1_Kpath, data1_Band = self.vpBandData_soc
                data2_Kpath, data2_Band = self.ndBandData_soc
            case "ncl":
                data1_Kpath, data1_Band = self.vpBandData_ncl
                data2_Kpath, data2_Band = self.ndBandData_ncl
        
        caseType1=label[0]
        caseType2=label[1]

        data1_Kpath = self._NormalizeTool(data1_Kpath)
        lineStyle1  = '-'

        data2_Kpath = self._NormalizeTool(data2_Kpath)
        lineStyle2  = '--'

        xlabel, xtick = self.vpLabelData_soc
        xtick = self._NormalizeTool(xtick)

        # >>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ax = plt.subplot(111)
        self._PlotTool(ax, kpath=data1_Kpath, band=data1_Band, 
                       lineinfo=(lineStyle1, caseType1, color[-1]))
        self._PlotTool(ax, kpath=data2_Kpath, band=data2_Band, 
                       lineinfo=(lineStyle2, caseType2, color[1]))
        
        ax.set_ybound(lower=boundary[0], upper=boundary[1])
        ax.set_xticks(ticks=xtick,
                      labels=xlabel)
        if title:
            ax.set_title(title)
        ax.grid('--')
        plt.legend(loc='upper right')
        if s:
            plt.savefig(f"CompareBand_{CompType}.png")
        else:
            plt.show()

    def CompBand_3kinds(self, wo_CalType, wo_socType, title, boundary=(-10, 15), s=False):
        match wo_CalType:
            case "vasp":
                match wo_socType:
                    case "soc":
                        caseType1 = "nd_soc"
                        data1_Kpath, data1_Band = self.ndBandData_soc

                        caseType2 = "nd_ncl"
                        data2_Kpath, data2_Band = self.ndBandData_ncl

                        caseType3 = "vp_ncl"
                        data3_Kpath, data3_Band = self.vpBandData_ncl
                        
                    case "ncl":
                        caseType1 = "nd_soc"
                        data1_Kpath, data1_Band = self.ndBandData_soc

                        caseType2 = "nd_ncl"
                        data2_Kpath, data2_Band = self.ndBandData_ncl

                        caseType3 = "vp_soc"
                        data3_Kpath, data3_Band = self.vpBandData_soc

                    case _:
                        print("No this kind of basis type")
            # ===============
            case "nanodcal":
                match wo_socType:
                    case "soc":
                        caseType1 = "vp_soc"
                        data1_Kpath, data1_Band = self.vpBandData_soc

                        caseType2 = "vp_ncl"
                        data2_Kpath, data2_Band = self.vpBandData_ncl

                        caseType3 = "nd_ncl"
                        data3_Kpath, data3_Band = self.ndBandData_ncl

                    case "ncl":
                        caseType1 = "vp_soc"
                        data1_Kpath, data1_Band = self.vpBandData_soc

                        caseType2 = "vp_ncl"
                        data2_Kpath, data2_Band = self.vpBandData_ncl

                        caseType3 = "nd_soc"
                        data3_Kpath, data3_Band = self.ndBandData_soc

                    case _:
                        print("No this kind of basis type")
            case _:
                print("No this kind of calculation type")
        
        data1_Kpath = self._NormalizeTool(data1_Kpath)
        lineStyle1  = '-'

        data2_Kpath = self._NormalizeTool(data2_Kpath)
        lineStyle2  = '--'

        data3_Kpath = self._NormalizeTool(data3_Kpath)
        lineStyle3  = ':'

        xlabel, xtick = self.vpLabelData_soc
        xtick = self._NormalizeTool(xtick)

        # >>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ax = plt.subplot(111)
        self._PlotTool(ax, kpath=data1_Kpath, band=data1_Band, 
                       lineinfo=(lineStyle1, caseType1, color[-1]))
        self._PlotTool(ax, kpath=data2_Kpath, band=data2_Band, 
                       lineinfo=(lineStyle2, caseType2, color[1]))
        self._PlotTool(ax, kpath=data3_Kpath, band=data3_Band, 
                       lineinfo=(lineStyle3, caseType3, color[4]))
        
        ax.set_ybound(lower=boundary[0], upper=boundary[1])
        ax.set_xticks(ticks=xtick,
                      labels=xlabel)
        if title:
            ax.set_title(title)
        ax.grid('--')
        plt.legend(loc='upper right')
        if s:
            plt.savefig(f"CompareBand-3kinds_wo-{wo_CalType}_{wo_socType}.png")
        else:
            plt.show()

    def CompBand_3kinds_selfLabel(self, wo_CalType, wo_socType, 
                                  title:str, label:tuple, 
                                  boundary=(-10, 15), s=False):
        match wo_CalType:
            case "vasp":
                match wo_socType:
                    case "soc":
                        data1_Kpath, data1_Band = self.ndBandData_soc
                        data2_Kpath, data2_Band = self.ndBandData_ncl
                        data3_Kpath, data3_Band = self.vpBandData_ncl
                        
                    case "ncl":
                        data1_Kpath, data1_Band = self.ndBandData_soc
                        data2_Kpath, data2_Band = self.ndBandData_ncl
                        data3_Kpath, data3_Band = self.vpBandData_soc

                    case _:
                        print("No this kind of basis type")
            # ===============
            case "nanodcal":
                match wo_socType:
                    case "soc":
                        data1_Kpath, data1_Band = self.vpBandData_soc
                        data2_Kpath, data2_Band = self.vpBandData_ncl
                        data3_Kpath, data3_Band = self.ndBandData_ncl

                    case "ncl":
                        data1_Kpath, data1_Band = self.vpBandData_soc
                        data2_Kpath, data2_Band = self.vpBandData_ncl
                        data3_Kpath, data3_Band = self.ndBandData_soc

                    case _:
                        print("No this kind of basis type")
            case _:
                print("No this kind of calculation type")
        
        caseType1=label[0]
        caseType2=label[1]
        caseType3=label[2]

        data1_Kpath = self._NormalizeTool(data1_Kpath)
        lineStyle1  = '-'

        data2_Kpath = self._NormalizeTool(data2_Kpath)
        lineStyle2  = '--'

        data3_Kpath = self._NormalizeTool(data3_Kpath)
        lineStyle3  = ':'

        xlabel, xtick = self.vpLabelData_soc
        xtick = self._NormalizeTool(xtick)

        # >>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ax = plt.subplot(111)
        self._PlotTool(ax, kpath=data1_Kpath, band=data1_Band, 
                       lineinfo=(lineStyle1, caseType1, color[-1]))
        self._PlotTool(ax, kpath=data2_Kpath, band=data2_Band, 
                       lineinfo=(lineStyle2, caseType2, color[1]))
        self._PlotTool(ax, kpath=data3_Kpath, band=data3_Band, 
                       lineinfo=(lineStyle3, caseType3, color[4]))
        
        ax.set_ybound(lower=boundary[0], upper=boundary[1])
        ax.set_xticks(ticks=xtick,
                      labels=xlabel)
        if title:
            ax.set_title(title)
        ax.grid('--')
        plt.legend(loc='upper right')
        if s:
            plt.savefig(f"CompareBand-3kinds_wo-{wo_CalType}_{wo_socType}.png")
        else:
            plt.show()

    def AllBand(self, title, boundary=(-10, 15), s=False):
        caseType1 = "vp_soc"
        data1_Kpath, data1_Band = self.vpBandData_soc

        caseType2 = "vp_ncl"
        data2_Kpath, data2_Band = self.vpBandData_ncl

        caseType3 = "nd_soc"
        data3_Kpath, data3_Band = self.ndBandData_soc

        caseType4 = "nd_ncl"
        data4_Kpath, data4_Band = self.ndBandData_ncl
        
        data1_Kpath = self._NormalizeTool(data1_Kpath)
        lineStyle1  = '-'

        data2_Kpath = self._NormalizeTool(data2_Kpath)
        lineStyle2  = '--'

        data3_Kpath = self._NormalizeTool(data3_Kpath)
        lineStyle3  = ':'

        data4_Kpath = self._NormalizeTool(data4_Kpath)
        lineStyle4  = '-.'

        xlabel, xtick = self.vpLabelData_soc
        xtick = self._NormalizeTool(xtick)

        # >>>>>>>>>>>>>>>>>>>>>>>  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        ax = plt.subplot(111)
        self._PlotTool(ax, kpath=data1_Kpath, band=data1_Band, 
                       lineinfo=(lineStyle1, caseType1, color[-1]))
        self._PlotTool(ax, kpath=data2_Kpath, band=data2_Band, 
                       lineinfo=(lineStyle2, caseType2, color[1]))
        self._PlotTool(ax, kpath=data3_Kpath, band=data3_Band, 
                       lineinfo=(lineStyle3, caseType3, color[4]))
        self._PlotTool(ax, kpath=data4_Kpath, band=data4_Band, 
                       lineinfo=(lineStyle4, caseType4, color[3]))
        
        ax.set_ybound(lower=boundary[0], upper=boundary[1])
        ax.set_xticks(ticks=xtick,
                      labels=xlabel)
        if title:
            ax.set_title(title)
        ax.grid('--')
        plt.legend(loc='upper right')
        if s:
            plt.savefig(f"CompareBand-ALL.png")
        else:
            plt.show()