import numpy as np
import matplotlib.pyplot as plt

from STT_Tool import Tools
color = Tools.ColorList()

class DataLoad():
    def __init__(self, dataPath, Border=2, mType='m', 
                 BasisParameter=[('t', 1), ('m', 3), ('Bext', 3)],
                 FurthParameter=[('Region-m1', 3), ('Region-m2', 3)]
                 ):
        
        data = np.loadtxt(dataPath)
        if mType=="m":
            data = data
        elif mType=="M":
            data = np.delete(data, [1, 2, 3], axis=1)

        dataDict = {}
        
        StartNum = 0
        EndNum   = StartNum + BasisParameter[0][1]
        dataDict[BasisParameter[0][0]] = data[:,StartNum:EndNum]

        StartNum = EndNum
        EndNum   = StartNum + BasisParameter[1][1]
        dataDict[BasisParameter[1][0]] = data[:,StartNum:EndNum]

        StartNum = EndNum
        EndNum   = StartNum + BasisParameter[2][1]
        dataDict[BasisParameter[2][0]] = np.round(data[:,StartNum:EndNum], Border)

        StartNum = EndNum
        for name, Dnum in FurthParameter:
            dataDict[name] = data[:,(StartNum):(StartNum+Dnum)]
            StartNum = StartNum+Dnum
        
        self.BasisParameter = BasisParameter
        self.FurthParameter = FurthParameter
        self.dataDict = dataDict

    def ALLDataLoad(self):
        dataDict = self.dataDict
        return dataDict
    
    def AngleLoad(self, BextAxis='z', zeroAxis='z', LoadDir='F',
                  CalParameter=['m', 'Region-m1', 'Region-m2']):
        
        dataDict = self.dataDict
        BextDict = {}
        BextArray = dataDict['Bext']
        match BextAxis:
            case 'x':
                BextInd = 0
            case 'y':
                BextInd = 1
            case 'z':
                BextInd = 2
            case _:
                print("Wrong Input, plz check the BextAxis")
        
        if LoadDir.upper()=='F':
            start = 0
            end   = (len(BextArray)//2)+1
        elif LoadDir.upper()=='B':
            start = (len(BextArray)//2)
            end   = (len(BextArray))+1
        else:
            print("!! No This direction !!")
        

        for i, Barray in enumerate(BextArray[start:end]):
            Car_r_array  = np.array([])
            theta_array  = np.array([])
            phi_array    = np.array([])
            for Cal in CalParameter:
                data = dataDict[Cal][start+i]
                mx, my, mz = data

                theta  = np.arctan((np.sqrt(mx**2+my**2))/mz)*180/np.pi
                if mx>0 and mz>0:
                    theta = theta
                elif mx>0 and mz<0:
                    theta = 180 + theta
                elif mx<0 and mz>0:
                    theta = -theta
                elif mx<0 and mz<0:
                    theta = theta
                elif mx==0 and mz<0:
                    theta = 180
                elif mx==0 and mz>0:
                    theta = 0
                    
                if mx == 0:
                    phi =0
                else:
                    phi = np.arctan(my/mx)*180/np.pi
                
                Car_r_array = np.append(Car_r_array, data)
                theta_array = np.append(theta_array, theta)
                phi_array   = np.append(phi_array, phi)

            Car_r_array = Car_r_array.reshape(len(Car_r_array)//3, 3)
            BextDict[str(Barray[BextInd])] = {'Car_r': Car_r_array,
                                              'theta': theta_array, 
                                              'phi': phi_array}
            # print(theta_array)
        return BextDict
            
            


    
    