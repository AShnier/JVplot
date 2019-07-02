# JVplot v0.6 by Adam Shnier 20/02/2018
# python 3.4.3
# Adam is not responsible for any calculation errors or local minima this script produces.
# In using this program you agree that you are responsible for checking that your data is 
# consistant with any output this program produces.

import matplotlib              #
matplotlib.use('agg')          # must be before importing pyplot
import matplotlib.pyplot as pl #
import sys, ctypes, xlwt
import configparser, os.path
from tkinter import filedialog, Tk, Label, Button, Entry, Frame, Scrollbar, Text, Menu
from tkinter import BooleanVar, DoubleVar
from scipy.stats import linregress
import numpy as np

global files

## Config
#*****************************************************
#***********************Config************************
#*****************************************************

def makeConfig(path):
    global opendir, pRsh, pRs, Vmin, Vmax, Imin, Imax, bLim, bSave2, bJPG, bxlraw
    global bArea, rArea, bIncid, rIncid, rIin, sIout, rIout, sWout, rWout, bPWR
    config = configparser.ConfigParser()
    config["Settings"] = {"opendir": opendir, "raw_in_excel":bxlraw, "bpwr":bPWR} #comma between values
    config["Processing"] = {"pRsh":pRsh, "pRs":pRs, "s_current_out":sIout, "r_current_out":rIout, "Convert_to_Amps":rIin, "s_power_out":sWout, "r_power_out":rWout,}
    config["Graphs"] = {"Vmin":Vmin, "Vmax":Vmax, "Imin":Imin, "Imax":Imax, "bLim":bLim, "bSave2":bSave2, "bJPG":bJPG}
    config["Area"] = {"bArea":bArea, "rArea":rArea}
    config["PCE"] = {"useIncidence":bIncid, "Incidence_W_per_cm2":rIncid}
    with open(path, "w") as configfile:
        config.write(configfile)

def doConfig(path):
    global opendir, pRsh, pRs, Vmin, Vmax, Imin, Imax, bLim, bSave2, bJPG, bxlraw
    global bArea, rArea, bIncid, rIncid, rIin, sIout, rIout, sWout, rWout, bPWR
    if os.path.isfile(path):
        config = configparser.ConfigParser()
        config.read(path)
        try:
            Vmin = float(config["Graphs"]["Vmin"]) #Vmin and Vmax are voltage values
            Vmax = float(config["Graphs"]["Vmax"])
            Imin = float(config["Graphs"]["Imin"]) #Imin and Imax are scaling factors for Isc
            Imax = float(config["Graphs"]["Imax"])
            #unfortunatly "False" wasnt recognised as bool, this is my work around
            if config["Graphs"]["bLim"] == "True" or config["Graphs"]["bLim"] == "1": #To implement Imin and Imax
                bLim = True
            else: bLim = False
            if config["Graphs"]["bSave2"] == "True" or config["Graphs"]["bSave2"] == "1": #To implement Imin and Imax
                bSave2 = True
            else: bSave2 = False
            if config["Graphs"]["bJPG"] == "True" or config["Graphs"]["bJPG"] == "1":
                bJPG = True
            else: bJPG = False
        except:
            print("Graph settings set as default")
        try:
            opendir = config["Settings"]["opendir"]
            if config["Settings"]["raw_in_excel"] == "True" or config["Settings"]["raw_in_excel"] == "1":
                bxlraw = True
            else: bxlraw = False
            if config["Settings"]["bPWR"] == "True" or config["Settings"]["bPWR"] == "1":
                bPWR = True
            else: bPWR = False
        except:
            print("Primary settings set as default")
        try:
            pRsh = int(config["Processing"]["pRsh"])
            pRs = int(config["Processing"]["pRs"])
            rIin = float(config["Processing"]["Convert_to_Amps"])
            sIout = str(config["Processing"]["s_current_out"])
            rIout = float(config["Processing"]["r_current_out"])
            sWout = str(config["Processing"]["s_power_out"])
            rWout = float(config["Processing"]["r_power_out"])
        except:
            print("Processing settings set as default")
        try:
            rArea = float(config["Area"]["rArea"])
            if config["Area"]["bArea"] == "True" or config["Area"]["bArea"] == "1":
                bArea = True
            else: bArea = False
        except:
            print("Area settings set as default")
        try:
            rIncid = float(config["PCE"]["Incidence_W_per_cm2"])
            if config["PCE"]["useIncidence"] == "True" or config["PCE"]["useIncidence"] == "1":
                bIncid = True
            else: bIncid = False
        except:
            print("PCE settings set as default")
    else:
        makeConfig(path)
        doConfig(path)

## Data handling
#*****************************************************
#*******************Data handling*********************
#*****************************************************

# CAUTION: this binary search may result in local minima
# a quick iterative search for accuracy
def absmin(array): # for array with a 1:1 mapping that passes through 0
    n1=0
    hi=len(array)-1
    lo=0
    if array[hi]*array[lo] < 0: #checks if this algorithm can be used
        for i in range(0, len(array)-1):
            mid=(int((hi+lo)/2))
            if array[mid]*array[lo]>0: #checking of same sign
                lo = mid
            else: #elif array[mid]*array[hi]>0:
                hi = mid
            if abs(lo-hi) <= 1: # truncation of int() would cause inf loop with 0
                break  
        n1 = mid
        n2 = 0
        m1=abs(array[mid])
        for i in range(mid-5,mid+5):
            try: #to handle possible out of range error
                if abs(array[i]) < m1:
                    m1=abs(array[i])
                    n2=i
            except:
                pass
        if abs(array[n2]) < abs(array[n1]):
            return n2
        else:
            return n1

    else:
        print("absmin exception")
        return None                                 

# averaged over 2 values for a more reliable answer
# It can be made more efficent by looking at every 3rd ([i]+[i+1])/2 before looking at closely at (i,i+1,i-1) 
def absmaxpwr(arrA,arrB,ci,cf):
    try:
        if ci > cf:
            cm = ci
            ci = cf
            cf = cm
        ## Can optimise calc for array
        v1=(arrA[ci]+arrA[ci+1])/2
        v2=(arrB[ci]+arrB[ci+1])/2
        m1=abs(v1*v2)
        for i in range(ci,cf-1):
            v1 = (arrA[i]+arrA[i+1])/2
            v2 = (arrB[i]+arrB[i+1])/2
            mloop = abs(v1*v2)
            if mloop > m1:
                m1 = mloop
        return m1
    except:
        ctypes.windll.user32.MessageBoxW(0, "max power error", "JVplot Error", 1)
        return -9876543210

def slopefit(arrI,arrV,iIsc,qRsh,iVoc,qRs): # Current on x-axis
    if qRsh > 0:
        x1=arrI[iIsc-qRsh:iIsc+qRsh+1]
        y1=arrV[iIsc-qRsh:iIsc+qRsh+1]
        slope1, intercept1, r_value1, p_value1, std_err1 = linregress(x1,y1)
    else:
        slope1, r_value1, std_err1 = -1,-1,-1
    if qRs > 0:
        x2=arrI[iVoc-qRs:iVoc+qRs+1]
        y2=arrV[iVoc-qRs:iVoc+qRs+1]
        slope2, intercept2, r_value2, p_value2, std_err2 = linregress(x2,y2)
    else:
        slope2, r_value2, std_err2 = -1,-1,-1
        
    if (r_value1 < 0.90) and (qRsh > 1):
        qRsh = qRsh-1
        return slopefit(arrI,arrV,iIsc,qRsh,iVoc,qRs) 
    if (r_value2 < 0.90) and (qRs > 1):
        qRs = qRs-1
        return slopefit(arrI,arrV,iIsc,qRsh,iVoc,qRs)
    ret = [slope1, (2*qRsh+1), r_value1, std_err1, slope2,(2*qRs+1), r_value2, std_err2]
    return ret
                    
## Read and process data to a set
def data2set(file):
    global arr1I, arr1V, array2, pRs, pRsh, bArea, rArea, bIncid, rIncid, lfiles
    global rIin, rIout, rWout
    ## Read in files
    print("reading files to sets")
    arr1I=[]
    arr1V=[]
    for j in range(0, lfiles):
        with open(file[j], "r") as ins:
            # X and y sets for arr1I,arr1V 
            arr1I.append([])
            arr1V.append([])
            for line in ins:
                try:
                    string=line.split("\t")
                    float(string[0])
                    float(string[1])
                    float(string[2])
                    float(string[3])
                    arr1I[j].append(float(string[1])) #current [2*j]
                    arr1V[j].append(float(string[3])) #voltage [2*j+1]
                    # at j=0 indexes 0 and 1 at j=1 indexes 2 and 3
                except:
                    pass
        
        # a fix for the truncation in our .lvm files
        try:
            v1trunc2 = arr1V[j][len(arr1V[j])-2] - arr1V[j][len(arr1V[j])-3]
            v1trunc = round(v1trunc2, 3) 
            if arr1V[j][len(arr1V[j])-1] == 0:
                if round(arr1V[j][len(arr1V[j])-1] - arr1V[j][len(arr1V[j])-2],3) != v1trunc:
                    arr1V[j][len(arr1V[j])-1] = (arr1V[j][len(arr1V[j])-2] + v1trunc2)
##                if (arr1V[j][len(arr1V[j])-2] == (1-v1trunc)):
##                    arr1V[j][len(arr1V[j])-1]=1
##                else:
##                    ctypes.windll.user32.MessageBoxW(0, "Truncation Error", "JVplot Error", 1)
            else:
                print("JVplot Error: fix truncation thing, 1")
        except:
            print("JVplot Error: fix truncation thing, 2")

    ## Convert input current
    if bArea:
        rMPin = rIin/rArea
    else:
        rMPin = rIin
    arr1V = np.array(arr1V)
    arr1I = np.array(arr1I)
    arr1I = np.multiply(arr1I,rMPin)
    
    ## Generate IV curve parameters
    print("processing data")
    array2 = []
    for j in range(0, lfiles):
        i=len(file[j])-3
        fname=""
        while i > 1:
            if file[j][i-1]=="/":
                fname=file[j][i:len(file[j])-4]
                break
            i=i-1
        Bpara = False
##        try:
        iIsc = absmin(arr1V[j]) # returns None type if array1 is unsuitable
        Isc = arr1I[j][iIsc]
        iVoc = absmin(arr1I[j])
        #Voc
        if arr1I[j][iVoc]*arr1I[j][iVoc-1] < 0:
            # voltage + volt_incriment(current_on_side/current_step)
            Voc=arr1V[j][iVoc-1]+(abs(arr1I[j][iVoc-1])*((arr1V[j][iVoc]-arr1V[j][iVoc-1])/(abs(arr1I[j][iVoc])+abs(arr1I[j][iVoc-1]))))
        else:
            Voc=arr1V[j][iVoc]+(abs(arr1I[j][iVoc])*((arr1V[j][iVoc+1]-arr1V[j][iVoc])/(abs(arr1I[j][iVoc+1])+abs(arr1I[j][iVoc]))))
        if Isc < -0.0000001 and Voc > 0 and iIsc != iVoc: #any current < 100 nA is assumed negligable
            pwr = absmaxpwr(arr1I[j],arr1V[j],iIsc,iVoc)
            ff = round(pwr/abs(Isc*Voc),3)
            setR = slopefit(arr1I[j],arr1V[j],iIsc,pRsh,iVoc,pRs)
            if (bArea and bIncid):
                rPCE = round((pwr/rIncid)*100, 3)
            else:
                rPCE = -1
            array2.append([file[j],Isc,Voc,pwr,ff]+[rPCE]+setR)
            Bpara = True
            
        else:
            print(1)
            array2.append([file[j],0,0,-1,-1]+[-1,-1,-1,-1,-1,-1,-1,-1,-1])   
##        except:
##            print(2)
##            array2.append([file[j],0,0,-1,-1]+[-1,-1,-1,-1,-1,-1,-1,-1,-1])

    ## Convert to output units (current)
    arr1I = np.multiply(arr1I,rIout)
    for j in range(0, lfiles):
        array2[j][1] = array2[j][1]*rIout
        array2[j][3] = array2[j][3]*rWout

    print("Generating images")
##    from multiprocessing.dummy import Pool # to speed up generating plots
##    pl = Pool(2)
##    arrTh = []
##    for i in range(0, lfiles):
##        arrTh.append(i)
##    print(arrTh)
##    pl.map(plot, arrTh)
##    pl.close()
    
##    pl.join()

    for j in range(0, lfiles):
        plot(j)

## Save plots by index
def plot(j):
    global bArea, arr1I, arr1V, array2, bLim, bSave2, Vmin, Vmax, Imin, Imax, files, sIout
    Bpara = False
    if array2[j][2] > 0:
        Bpara = True
        Isc = array2[j][1]
        Voc = array2[j][2]
        ff = array2[j][4]
    i=len(files[j])-3
    fname=""
    while i > 1:
        if files[j][i-1]=="/":
            fname=files[j][i:len(files[j])-4]
            break
        i=i-1
        
##    pl.close()
##    pl.cla()
##    pl.clf()
    pl.figure(fname)
    pl.plot(arr1V[j],arr1I[j])
    pl.grid(which='major')
    pl.axhline(linewidth=1, color="black")
    pl.axvline(linewidth=1, color="black")
    
    if Bpara:
        if bArea:
            pl.title(fname + "\n Voc:" + str(round(Voc,3)) + "V Jsc:" + str(round(Isc,4)) + sIout+"/cm^2 FF:" + str(ff))
        else:
            pl.title(fname + "\n Voc:" + str(round(Voc,3)) + "V Isc:" + str(round(Isc,4)) + sIout+" FF:" + str(ff))
    else:
        pl.title(fname)
    pl.xlabel("Voltage(V)")
    if bArea:
        pl.ylabel("Current density("+sIout+"/cm^2)")
    else:
        pl.ylabel("Current("+sIout+")")
        ## Full graph
    if bJPG: # If running full program
        Psaved = False
        Pcount = 0
        if not os.path.exists(files[j]+"Full.png"):
            pl.savefig(files[j]+"Full.png")
            Psaved = True
        elif not(bSave2):
            while not Psaved:
                Pcount += 1
                if not os.path.exists(files[j]+"Full("+str(Pcount)+").png"):
                    pl.savefig(files[j]+"Full("+str(Pcount)+").png")
                    Psaved = True
                else:
                    if Pcount > 20:
                        print("Picture could not be saved")
                        Psaved = True # ends loop
                        
        ## Focused graph
        if bLim and Bpara and Isc < -0.0000001 and Voc > 0:
            pl.xlim(-1*abs(Vmin*Voc), abs(Vmax*Voc))
            pl.ylim(-1*abs(Imin*Isc), abs(Imax*Isc))
            Psaved = False
            Pcount = 0
            if not os.path.exists(files[j]+"Zoom.png"):
                pl.savefig(files[j]+"Zoom.png")
                Psaved = True
            elif not(bSave2):
                while not Psaved:
                    Pcount += 1
                    if not os.path.exists(files[j]+"Zoom("+str(Pcount)+").png"):
                        pl.savefig(files[j]+"Zoom("+str(Pcount)+").png")
                        Psaved = True
                    else:
                        if Pcount > 20:
                            print("Picture could not be saved")
                            Psaved = True # ends loop
    pl.close()
    
#*****************************************************
#***********************Excel*************************
#*****************************************************
## Excel
def excel():
    global arr1I, arr1V, array2, pRs, pRsh, opendir, bxlraw, bArea, rArea, bIncid, rIncid
    global sIout, rIout, sIin, rIin, sWout, rWout, bPWR 
    print("Writing to excel")
    book = xlwt.Workbook()
    sheets = []
    sheets.append(book.add_sheet("Summary"))
        #array2.append([file[j],Isc,Voc,pwr,ff])
        # Voc and Isc order swopped when written to excel
    sheets[0].write(0, 1, "file") # row, column, value
    sheets[0].write(0, 2, "Voc (V)")

    r = 0
    if bArea:            
        sheets[0].write(0, 3, "Jsc ("+sIout+"/cm^2)")
        sheets[0].write(0, 4, "pwr ("+sWout+"/cm^2)")
        sheets[0].write(0, 5, "FF")
        if bIncid:
            r = 1
            sheets[0].write(0, 6, "PCE %")
            sheets[0].write(4, 15+r, "Incidence/W cm^-2")
            sheets[0].write(5, 16+r, rIncid)
        #1/slope1, slope1, r_value1, std_err1, 1/slope2, slope2, r_value2, std_err2
        sheets[0].write(0, 6+r, "Approx Rsh (omega cm^2)")
        sheets[0].write(0, 10+r, "Approx Rs (omega cm^2)")
        sheets[0].write(3, 15+r, "Area/cm^2")
        sheets[0].write(3, 16+r, rArea)
    else:
        sheets[0].write(0, 3, "Isc ("+sIout+")")
        sheets[0].write(0, 4, "pwr ("+sWout+")")
        sheets[0].write(0, 5, "FF")
        #1/slope1, slope1, r_value1, std_err1, 1/slope2, slope2, r_value2, std_err2
        sheets[0].write(0, 6, "Approx Rsh")
        sheets[0].write(0, 10, "Approx Rs")
    sheets[0].write(0, 7+r, "Points used")
    sheets[0].write(0, 8+r, "r_value")
    sheets[0].write(0, 9+r, "std_error")
    sheets[0].write(0, 11+r, "Points used")
    sheets[0].write(0, 12+r, "r_value")
    sheets[0].write(0, 13+r, "std_error")
    sheets[0].write(0, 15+r, "Points used for gradients")
    sheets[0].write(1, 15+r, "Shunt (omega)")
    sheets[0].write(1, 16+r, 2*pRsh+1)
    sheets[0].write(2, 15+r, "Series (omega)")
    sheets[0].write(2, 16+r, 2*pRs+1)
    sheets[0].write(6, 15+r, "Unit conversions")
    sheets[0].write(7, 15+r, "Current input conversion")
    sheets[0].write(8, 15+r, "Multiplier")
    sheets[0].write(8, 16+r, rIin)
    sheets[0].write(9, 15+r, "-For the conversion to amps")
    sheets[0].write(10, 15+r, "Current output conversion")
    sheets[0].write(11, 15+r, "Multiplier")
    sheets[0].write(11, 16+r, rIout)
    sheets[0].write(12, 15+r, "Unit name")
    sheets[0].write(12, 16+r, sIout)
    sheets[0].write(13, 15+r, "Power output conversion")
    sheets[0].write(14, 15+r, "Multiplier")
    sheets[0].write(14, 16+r, rWout)
    sheets[0].write(15, 15+r, "Unit name")
    sheets[0].write(15, 16+r, sWout)
    

    q = 0
    for n in range(0,len(array2)):
        k = len(array2[n][0])-4
        Rfile="file "+str(n) #for the unlikley possibility k <= 0
        while (k > 0):
            if array2[n][0][k] == "/":
                Rfile = (array2[n][0][k+1:len(array2[n][0])-4])
                k = 0 # ends loop
            else:
                k -= 1
        # Summary values
        sheets[0].write(n+1, 0, Rfile)
        sheets[0].write(n+1, 1, array2[n][0]) # row, column, value
        sheets[0].write(n+1, 2, array2[n][2])
        sheets[0].write(n+1, 3, array2[n][1])
        sheets[0].write(n+1, 4, array2[n][3])
        sheets[0].write(n+1, 5, array2[n][4])
        if bArea and bIncid:
            sheets[0].write(n+1, 6, array2[n][5])
        try:
            for i in range (6,len(array2[n])): # the last entry of array2 for PCE
                sheets[0].write(n+1,i+r, array2[n][i])
        except:
            pass
        t=0
        if bxlraw:
            if n == 20*q:
                sheets.append(book.add_sheet("Raw_Data"+str(q)))
                q += 1
            t=(q-1)*20
            sheets[q].write(2, 3*(n-t), "V") # row, column, value
            if bArea:
                sCur = "J ("+sIout+"/cm^2)"
            else:
                sCur = "A ("+sIout+")"
            sheets[q].write(2, 3*(n-t)+1, sCur)
            sheets[q].write(0, 3*(n-t), Rfile)
            sheets[q].write(1, 3*(n-t), array2[n][0])
            if bPWR:
                if bArea:
                    sheets[q].write(2, 3*(n-t)+2, "pwr ("+sWout+"/cm^2)")
                else:
                    sheets[q].write(2, 3*(n-t)+2, "pwr ("+sWout+")")
                
            for p in range(0,len(arr1V[n])):
                sheets[q].write(p+3, 3*(n-t), arr1V[n][p]) # row, column, value
                sheets[q].write(p+3, 3*(n-t)+1, arr1I[n][p])
                if bPWR:
                    sheets[q].write(p+3, 3*(n-t)+2, (arr1V[n][p]*arr1I[n][p]))
                    

    Esaved = False
    nameadd = ""
    lopendir = len(opendir)
    for s in range(0, lopendir):
        if opendir[lopendir-s-1] == "/":
            nameadd = opendir[lopendir-s-1:lopendir]+"_JVplot" # includes "/"
            break
    if not os.path.exists(opendir+nameadd+".xls"):
        book.save(opendir+nameadd+".xls")
        print(opendir+nameadd+".xls")
        Esaved = True
    else:
        Ecount = 0
        while not Esaved:
            Ecount += 1
            if not os.path.exists(opendir+nameadd+"("+str(Ecount)+").xls"):
                book.save(opendir+nameadd+"("+str(Ecount)+").xls")
                Esaved = True
            else:
                if Ecount > 20:
                    print("Excel file could not be saved")
                    Esaved = True # ends loop

#*****************************************************
#**************Interface and direct run***************
#*****************************************************
        
## interface
class gui_tk(Tk):    
    def __init__(self,parent):
        Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()

    def initialize(self):
        global files, opendir, bSave2, bJPG, bxlraw, rArea, bArea, rIncid, bIncid, bPWR
        def setlist():
            t1wid = 60
            #text1.delete(1.0,"end")
            for m in range(0, len(files)):
                text1.insert("end", files[m] + "\n")
                if len(files[m]) > t1wid:
                    t1wid = len(files[m])
            lfiles=len(files)
            if lfiles > 10:
                t1hi = 20
            else:
                t1hi = 10
            text1.config(width=t1wid, height=t1hi)

        ## set opendir # not used for folder2set
        def setopendir():
            global opendir, files
            lfiles=len(files)
            ilfile = []
            lfile0 = files[0]
            lfileL = files[lfiles-1]
            iFind = 0
            while True :
                iFind += 1
                try:
                    if lfile0[iFind] == "/":
                        if lfileL[iFind] == "/":
                            ilfile.append(iFind)
                        else:
                            break
                except:
                    break
            c1 = len(ilfile)
            for k in range (0,c1):
                print(lfile0[0:ilfile[c1-k-1]])
                if lfile0[0:ilfile[c1-k-1]] == lfileL[0:ilfile[c1-k-1]]:
                    opendir = lfile0[0:ilfile[c1-k-1]]
                    break
                
        def clearfiles():
            global files
            files = ()
            text1.delete(1.0,"end")
            text1.config(width=60, height=15)
            menu1.entryconfig("Run", state="disabled")

        def clearlist():
            text1.delete(1.0,"end")
            text1.config(width=60, height=15)
            menu1.entryconfig("Run", state="disabled")
            
        ## Read in files
        def folder2set():
            global files, opendir
            opendir = filedialog.askdirectory(initialdir = opendir)
            filetemp = list(files)
            for root, dirs, wfiles in os.walk(opendir):
                for name in wfiles:
                    lname = len(name)
                    if name[lname-3:lname].lower() == "lvm":
                        filetemp.append(os.path.join(root, name).replace("\\","/"))
            files = tuple(filetemp)
            clearlist()
            setlist()
            menu1.entryconfig("Run", state="normal")

        def files2set():
            global files
            files2 = (filedialog.askopenfilenames(initialdir = opendir, filetypes = (("Labview","lvm"),("Text files","txt"),("all files","*.*"))))
            files = tuple(files+files2)
            clearlist()
            setlist()
            setopendir()
            menu1.entryconfig("Run", state="normal")

        

        def OpenConfig():
            os.startfile("C:\\ProgramData\JVplot\JVconfig.ini")

        def areaon():
            global bArea, bIncid
            bArea = bvarArea.get()
            if bArea:
                tArea.grid(column=1, row=2)
                lArea.grid(column=0, row=2)
                if bIncid:
                    tIncid.grid(column=3, row=2)
                    lIncid.grid(column=2, row=2)
            else:
                tArea.grid_remove()
                lArea.grid_remove()
                tIncid.grid_remove()
                lIncid.grid_remove()
        def incidon():
            global barea, bIncid
            bIncid = bvarIncid.get()
            if bIncid and bArea:
                tIncid.grid(column=3, row=2)
                lIncid.grid(column=2, row=2)
            else:
                tIncid.grid_remove()
                lIncid.grid_remove()

        def run():
            global files, bSave2, opendir, bJPG, bxlraw, rArea, bArea, lfiles, bPWR
            lfiles = len(files)
            if lfiles > 0:
                bSave2 = varSave2.get()
                bJPG = varJPG.get()
                bxlraw = varxlraw.get()
                rArea = rvarArea.get()
                bArea = bvarArea.get()
                rIncid = rvarIncid.get()
                bIncid = bvarIncid.get()
                bPWR = bvarPWR.get()
                data2set(files)
                excel()
                makeConfig("C:\\ProgramData\JVplot\JVconfig.ini")
                ctypes.windll.user32.MessageBoxW(0, "Finished, check the relevant directory to see output", "JVplot Completed", 1)
                os.startfile(opendir)

        def close():
            global arr1I, arr1V, array2
            arr1I, arr1V, array2 = [], [], []
            app.destroy()
                
        ## GUI objects
        self.grid()
        self.minsize(width=400, height=200)
            
        menu1 = Menu(self)
        filemenu = Menu(menu1, tearoff=0)
        filemenu.add_command(label="Open Files", command=files2set)
        filemenu.add_command(label="Read Folder", command=folder2set)
        filemenu.add_separator()
        filemenu.add_command(label="Clear Files", command=clearfiles)
        menu1.add_cascade(label="Files", menu=filemenu)
        optmenu = Menu(menu1, tearoff=0)

        varxlraw=BooleanVar()
        if bxlraw:
            varxlraw.set(True)
        else:
            varxlraw.set(False)
        optmenu.add_checkbutton(label="Raw data in Excel", variable=varxlraw ,onvalue = True,offvalue = False)
        bvarPWR=BooleanVar()
        if bPWR:
            bvarPWR.set(True)
        else:
            bvarPWR.set(False)
        optmenu.add_checkbutton(label="power in raw data", variable=bvarPWR ,onvalue = True,offvalue = False)
        optmenu.add_separator()
        
        varSave2=BooleanVar()
        if bSave2:
            varSave2.set(True)
        else:
            varSave2.set(False)
        optmenu.add_checkbutton(label="New jpg's only", variable=varSave2 ,onvalue = True,offvalue = False)
        varJPG = BooleanVar()
        if bJPG:
            varJPG.set(True)
        else:
            varJPG.set(False)
        optmenu.add_checkbutton(label="Save images", variable=varJPG ,onvalue = True,offvalue = False)
        optmenu.add_separator()
        
        bvarArea=BooleanVar()
        if bArea:
            bvarArea.set(True)
        else:
            bvarArea.set(False)
        optmenu.add_checkbutton(label="Use cell area", variable=bvarArea ,onvalue = True,offvalue = False, command=areaon)
        bvarIncid=BooleanVar()
        if bIncid:
            bvarIncid.set(True)
        else:
            bvarIncid.set(False)
        optmenu.add_checkbutton(label="Calculate PCE", variable=bvarIncid ,onvalue = True,offvalue = False, command=incidon)
        optmenu.add_separator()
        
        optmenu.add_command(label="Open config", command=OpenConfig)
        menu1.add_cascade(label="Options", menu=optmenu)
        menu1.add_command(label="Run", command=run, state="disabled")
        
        helpmenu = Menu(menu1, tearoff=0)
        helpmenu.add_command(label="Email Adam")
        menu1.add_cascade(label="Help", menu=helpmenu)
        menu1.add_command(label="Exit", command=close)
        
        self.config(menu=menu1)

        label1 = Label(self, text="Be patient after you click run, The file directories will show progress")
        label1.grid(column=0, row=1, columnspan=4)

        lArea = Label(self, text="Area /cm2")
        rvarArea = DoubleVar()
        rvarArea.set(rArea)
        tArea = Entry(self, textvariable=rvarArea)
        lIncid = Label(self, text="Incidence /W per cm2")
        rvarIncid = DoubleVar()
        rvarIncid.set(rIncid)
        tIncid = Entry(self, textvariable=rvarIncid)
        if bArea:
            tArea.grid(column=1, row=2)
            lArea.grid(column=0, row=2)
            if bIncid:
                tIncid.grid(column=3, row=2)
                lIncid.grid(column=2, row=2) 

        frame1 = Frame(self,width=60, height=15)
        frame1.grid(column=0, row=3, columnspan=5)
        scrollbar = Scrollbar(frame1) 
        scrollbar.pack(side="right", fill="y")
        text1 = Text(frame1, width=60, height=15, yscrollcommand=scrollbar.set)
        text1.pack()
        scrollbar.config(command=text1.yview)
        
## Read or make config
global opendir, pRsh, pRs, Vmin, Vmax, Imin, Imax, bLim, files, bSave2, bJPG
global bArea, rArea, bIncid, rIncid, rIin, rIout, rWout, sIout, sWout, bPWR
opendir = "C:\\users"
pRsh, pRs = 2, 2
Vmin, Vmax, Imin, Imax = 0.2, 1.2, 1.2, 0.2 #Imin and Imax are scaling factors for Isc
bLim, bSave2, bJPG, bxlraw = True, True, True, True
bIncid, bArea = False, False
files = ()
rArea, rIncid = 0.07, 0.1
rIin, rIout, rWout = 1, 1000, 1000
sIout, sWout = "mA", "mW"
bPWR = False
if not os.path.exists("C:\\ProgramData\JVplot"):
    os.makedirs("C:\\ProgramData\JVplot")
doConfig("C:\\ProgramData\JVplot\JVconfig.ini")

## calls interface or plots
##if __name__ == "__main__" and len(sys.argv) != 3: # determine run type
app = gui_tk(None)
app.title('JVplot v0.5.6 by Adam Shnier (Sept2017)')
app.lift()
app.mainloop()
##else:
##    files = [sys.argv[2]]
##    data2set(files)



#***************
#future features
#***************
#redundancy between Rfile and fname by finding name seperatly
# In Excel show current in mA and power, mW

#Seperate current and voltage to seperate arrays
    #use initial, final-1, step size to compare if voltage ranges are the same
#My use of numpy is pathetic, Fix
#seperate sections read process plot
