import numpy as np
import pandas as pd
import re

class data(object):
    def __init__(self,filelocation,log):
        ###get stress strain data
        headers = ['col1','col2','col3','col4','col5','col6','col7','col8','col9','col10','col11','col12']
        data = pd.read_csv(filelocation, sep="\s+",names = headers )
        df = pd.DataFrame(data)
        maxRow = df['col1'].size
        for i in range(8):
            if i in [1,2,6]:
                n = 0
                while(8*n<maxRow):
                    for j in range(1,11):
                        col = 'col' + str(12-j)
                        colShifted = 'col' + str(12-j+1)
                        every = i + 8*n
                        df.at[every,colShifted] = df[col][every]
                    n+=1
            else:
                n = 0
                while(8*n<maxRow):
                    every = i + 8*n
                    df.at[every,'col1'] = df['col1'][every]+'_'+df['col2'][every]
                    n+=1

        df = df.drop(['col2','col12'],axis=1)

        ColNames = {'col1':'Variable','col3':'e11','col4':'e12','col5':'e13'
                    ,'col6':'e21','col7':'e22','col8':'e23','col9':'e31'
                    ,'col10':'e32','col11':'e33'}
        self.__df = df.rename(index=str, columns = ColNames)

        ###get the log data
        self.__IterationCount=[]
        with open(log) as f:
            for line in f:
                pattern = re.compile(r'(?<=iteration: )\d+')
                self.__IterationCount += pattern.findall(line)
        self.__IterationCount = list(map(int,self.__IterationCount))
        
    def getDataFrame(self):
        return self.__df
    
    def getIterCount(self):
        return self.__IterationCount
    
    def getValue(self,**kwargs):
        #Value = {}
        #Value = []
        maxRow = len(self.__df.iloc[:,0])
        for key, arg in kwargs.items():
            temp = []
            for i in range(maxRow):
                if self.__df['Variable'][i]==key:
                    temp.append(self.__df[arg][i])
            #Value.append(temp)
        return temp
    
    def getPureElasticStress(self,eij='e12'):
        #Value = self.getValue(Stress=eij,Strain_mean=eij)
        stress = self.getValue(Stress=eij)
        strain = self.getValue(Strain_mean=eij)
        tangent = (stress[2]-stress[0])/2/(strain[1]-strain[0])
        Elastic = [stress[0]]
        updateStress = stress[0]
        for i in range(1,len(strain)):
            stepsize = strain[i]-strain[i-1]
            updateStress += tangent*stepsize
            Elastic.append(updateStress)
        return Elastic
        
    def getCpToLinearErr(self):
        elastic = self.getPureElasticStress()
        plastic = self.getValue(Stress='e12')
        err = [0]
        for i in range(1,len(plastic)):
            err.append(abs(elastic[i]-plastic[i])/elastic[i]*100)
            
        return err
            
    def getPointCoord(self,percentage):
        stress = self.getValue(Stress='e12')
        strain = self.getValue(Strain='e12')
        pointStress = np.percentile([np.min(stress),np.max(stress)],percentage)

        for i in range(len(stress)):
            if stress[i]>pointStress:
                x0 = stress[i-1]
                x1 = stress[i]
                count = i-1
                break
        pointStrain = strain[count] + (pointStress-x0)/(x1-x0)*(strain[count+1]-strain[count])
        return (pointStrain,pointStress,count)

regularisation = ["bardella_lin","tanh_lin","powerlaw_sech"]
gamma0dot = ['0.01','0.0025','0.0006','0.00024','0.00006']
p = {"bardella_lin":['0.4','0.16','0.064','0.0256','0.01'],
     "tanh_lin":['0.16','0.064'],
     "powerlaw_sech":['0.16']}


with open("data.txt","w")as file:
    for r in regularisation:
        for i in p[r]:
            IterAtPoint80 = []
            ErrAtPoint80 = []
            IterAtPoint90 = []
            ErrAtPoint90 = []
            for g in gamma0dot:
                text = "stress_strain_"+r+"_"+g+"_"+str(i)+".txt"
                log = "deallog_"+r+"_"+g+"_"+str(i)
                filename = '/home/conner/Downloads/Results/' + text
                logname = '/home/conner/Downloads/Results/' + log
                dataset = data(filename,logname)

                name = r+"_"+g+"_"+i
               
                file.write(name) #stress-strain data

                file.write("\r\n")
                file.write("stress-strain:")
                file.write("\r\n")
                stress = dataset.getValue(Stress = "e12")
                strain = dataset.getValue(Strain_mean="e12")
                for x in zip(strain,stress):
                    file.write(str(x))
                    file.write("\r\n")
                    
                file.write("\r\n")
                file.write("plastic strain - strain:")
                file.write("\r\n")
                plasticStrain = dataset.getValue(Plastic_strain_mean="e12")
                for x in zip(strain,plasticStrain):
                    file.write(str(x))
                    file.write("\r\n")
                    
                file.write("\r\n")
                file.write("IterationCount-strain:")
                file.write("\r\n")
                IterCount = dataset.getIterCount()
                for x in zip(strain,IterCount):
                    file.write(str(x))
                    file.write("\r\n")
                    
                file.write("\r\n")
                file.write("error - strain:")
                file.write("\r\n")
                err = dataset.getCpToLinearErr()
                for x in zip(strain,err):
                    file.write(str(x))
                    file.write("\r\n")

                file.write("\r\n")
                file.write("Attentioned Point:80% of Stress")
                file.write("\r\n")
                point80 = dataset.getPointCoord(80)
                pointCoord80 = (point80[0],point80[1])
                neighborLoc80 = point80[2]
                file.write(str(pointCoord80))
                file.write("\r\n")
                file.write("------------------------\n")

                ErrAtPoint80.append(err[neighborLoc80])
                IterAtPoint80.append(IterCount[neighborLoc80])

                file.write("\r\n")
                file.write("Attentioned Point:90% of Stress")
                file.write("\r\n")
                point90 = dataset.getPointCoord(90)
                pointCoord90 = (point90[0],point90[1])
                neighborLoc90 = point90[2]
                file.write(str(pointCoord90))
                file.write("\r\n")
                file.write("------------------------\n")

                ErrAtPoint90.append(err[neighborLoc90])
                IterAtPoint90.append(IterCount[neighborLoc90])

            file.write(">>>>>>>>>>>>>>>>>\n")
            file.write("Error at attentioned point(80%):")
            file.write("\r\n")
            for x in zip(list(map(float,gamma0dot)),ErrAtPoint80):
                file.write(str(x))
                file.write("\r\n")

            file.write("\r\n")
            file.write("Iteraion Count at attentioned point(80):")
            file.write("\r\n")
            for x in zip(list(map(float,gamma0dot)),IterAtPoint80):
                file.write(str(x))
                file.write("\r\n")
            file.write("<<<<<<<<<<<<<<<<<\n")
            file.write("-----------------\n")

            file.write(">>>>>>>>>>>>>>>>>\n")
            file.write("Error at attentioned point(90%):")
            file.write("\r\n")
            for x in zip(list(map(float,gamma0dot)),ErrAtPoint90):
                file.write(str(x))
                file.write("\r\n")
            file.write("\r\n")

            file.write("\r\n")
            file.write("Iteraion Count at attentioned point(90):")
            file.write("\r\n")
            for x in zip(list(map(float,gamma0dot)),IterAtPoint90):
                file.write(str(x))
                file.write("\r\n")
            file.write("<<<<<<<<<<<<<<<<<\n")
            


