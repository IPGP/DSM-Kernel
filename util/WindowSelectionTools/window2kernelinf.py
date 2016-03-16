'''
Created on 11 aout 2013

@author: hugojaegler
'''
import numpy as np
import matplotlib.pyplot as plt
from pylab  import vlines
from obspy.core import read
import subprocess
from Tkinter import Tk, Canvas, Frame, Scrollbar, Label, Entry, Button, Menu, Menubutton, Listbox, StringVar, IntVar, RAISED, Checkbutton
from shutil import Error
             
class Path:
    def __init__(self, root):
        
        self.canvas = Canvas(root, borderwidth=1, background="#ffffff")
        self.frame = Frame(self.canvas, background="#ffffff")
        self.vsb = Scrollbar(root, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)
        self.vsb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)
        
        self.canvas.create_window((4,4), window=self.frame, anchor="nw", tags="self.frame")

        self.frame.bind("<Configure>", self.OnFrameConfigure)

        self.data()
        
    def data(self): 
        
        global textPath
        textPath = StringVar()
        global text0a
        text0a = StringVar()
        global text0b 
        text0b = StringVar()
        global text2a 
        text2a = StringVar()
        global text3 
        text3 = StringVar()
        global alphaVar
        alphaVar = IntVar()
        global betaVar 
        betaVar = IntVar()
        global allVar
        allVar = IntVar()
        global text6a
        text6a = "0"
        global filterVar
        filterVar = IntVar()
        global text6b
        text6b = StringVar()
        global t1x
        t1x = ""
        global t2x
        t2x = ""
        global t3x
        t3x = ""
        global t4x
        t4x = ""
        global text8_0
        text8_0 = StringVar()
        global text8_1
        text8_1 = StringVar()
        
        Label(self.frame,text="Path ? ").grid(row=0, column=0)
        Entry(self.frame,textvariable=textPath).grid(row=1, column=0)
        Button(self.frame, text="Valider et afficher", command = affiche_recap).grid(row=1, column=1)
    
        Label(self.frame, text="Green function database information file\n (for a certain depth only for the instance) ?").grid(row=3)
        Entry(self.frame, textvariable=text0a).grid(row=4)
        
        Label(self.frame, text="Output directory (parentdir) ?").grid(row=5)
        Entry(self.frame, textvariable=text0b).grid(row=6)
            
        Label(self.frame, text="Phase name ?").grid(row=9)
        Entry(self.frame, textvariable=text3).grid(row=10)
        
        def afficheAlpha():
            seismicPara["text"]="alpha"
            betaVar.set(0)
            allVar.set(0)
        def afficheBeta():
            seismicPara["text"]="beta"
            alphaVar.set(0)
            allVar.set(0)
        def afficheAll():
            seismicPara["text"]="all"
            alphaVar.set(0)
            betaVar.set(0)
        
        seismicPara = Menubutton(self.frame, text="Seismic Parameter", relief=RAISED)
        seismicPara.grid(row=0)
        seismicPara.menu = Menu(seismicPara, tearoff = 0)
        seismicPara["menu"] = seismicPara.menu

        
        seismicPara.menu.add_checkbutton(label="alpha", variable = alphaVar, command = afficheAlpha)
        seismicPara.menu.add_checkbutton(label="beta", variable = betaVar, command = afficheBeta)
        seismicPara.menu.add_checkbutton(label="all", variable = allVar, command = afficheAll)
        seismicPara.grid(row=11)
        
        
        
        Label(self.frame, text="Filter name ?").grid(row=12)
        Entry(self.frame, textvariable=text6b).grid(row=13)
        
        
        
        Label(self.frame, text="time window t1 ?").grid(row=14)
        Labelt1 = Label(self.frame, text="-->").grid(row=15)
        Button(self.frame, text="time 1", command=self.time1).grid(row=15, column=1)
        
        Label(self.frame, text="time window t2 ?").grid(row=16)
        Labelt1 = Label(self.frame, text="-->").grid(row=17)
        Button(self.frame, text="time 2", command=self.time2).grid(row=17, column=1)
        '''
        Label(self.frame, text="time window t3 ?").grid(row=18)
        Labelt1 = Label(self.frame, text="-->").grid(row=19)        
        Button(self.frame, text="time 3", command=self.time3).grid(row=19, column=1)
        
        Label(self.frame, text="time window t4 ?").grid(row=20)
        Labelt1 = Label(self.frame, text="-->").grid(row=21)
        Button(self.frame, text="time 4", command=self.time4).grid(row=21, column=1)
        '''
        def affiche0():
            convertPara["text"]="No conversion"
            text8_1.set(0)
            
        def affiche1():
            convertPara["text"]="Conversion"
            text8_0.set(0)
    
        convertPara = Menubutton(self.frame, text="Geodetic latitude to geocentric latitude conversion", relief=RAISED)
        convertPara.grid(row=0)
        convertPara.menu = Menu(convertPara, tearoff = 0)
        convertPara["menu"] = convertPara.menu

        convertPara.menu.add_checkbutton(label="No conversion", variable = text8_0, command = affiche0)
        convertPara.menu.add_checkbutton(label="Conversion", variable = text8_1, command = affiche1)
        
        convertPara.grid(row=22)
        b = Checkbutton(self.frame, text = "apply filter", variable = filterVar)
        b.grid(row=23, column = 0)
        Button(self.frame, text="continue", command=self.quitter).grid(row=23, column=1)
        
    def time1(self):
        global t1x
        global t1y
        t1x, t1y = Pointage()
        print type(t1x)
        print t1y

    def time2(self):
        global t2x
        global t2y
        t2x, t2y = Pointage()
        print t2x
        print t2y
            
    def time3(self):
        t3x, t3y = Pointage()
        print t3x
        print t3y
            
    def time4(self):
            t4x, t4y = Pointage()
            print t4x
            print t4y
            
            
    def quitter(self):
        root.destroy()
    
    
    def OnFrameConfigure(self, event):
        '''Reset the scroll region to encompass the inner frame'''
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        
class RecapCalculs:
    '''
    Interface graphique recapitulant les caracteristique du sismogramme
    presentant les options de filtrage et de calculs du noyau de sensibilite
    '''
    def __init__(self,root):
        self.canvas = Canvas(root, borderwidth=1, background="#ffffff")
        self.frame = Frame(self.canvas, background="#ffffff")
        self.vsb = Scrollbar(root, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)
        self.vsb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)
        self.canvas.create_window((4,4), window=self.frame, anchor="nw", tags="self.frame")

        self.frame.bind("<Configure>", self.OnFrameConfigure)

        self.data()
    def data(self):
        
        self.message = Label(self.frame, text="Recapitulatif du sismogramme").grid(row=0)
       
        self.recap = Listbox(self.frame, height = 15, width = 50)
           
        self.recap.insert(1, "network: {}\n".format(X[0].stats.network))
        self.recap.insert(2, "station: {}\n".format(X[0].stats.station))
        self.recap.insert(3, "location: {}\n".format(X[0].stats.location))
        self.recap.insert(4, "channel: {}\n".format(X[0].stats.channel))
        self.recap.insert(5, "start time: {}\n".format(X[0].stats.starttime))
        self.recap.insert(6, "end time: {}\n".format(X[0].stats.endtime))
        self.recap.insert(7, "sampling rate: {}\n".format(X[0].stats.sampling_rate))
        self.recap.insert(8, "delta: {}\n".format(X[0].stats.delta))
        self.recap.insert(9, "number points: {}\n".format(X[0].stats.npts))
        self.recap.insert(10, "calibration: {}\n".format(X[0].stats.calib))
        self.recap.insert(11, "event latitude: {}\n".format(X[0].stats.sac.evla))
        self.recap.insert(12, "event longitude: {}\n".format(X[0].stats.sac.evlo))
        self.recap.insert(13, "event depth: {}\n".format(X[0].stats.sac.evdp))
        self.recap.insert(14, "station latitude: {}\n".format(X[0].stats.sac.stla))
        self.recap.insert(15, "station longitude: {}\n".format(X[0].stats.sac.stlo))
        self.recap.grid(row=0)       
        
    def OnFrameConfigure(self, event):
        '''Reset the scroll region to encompass the inner frame'''
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))
  
        
def affiche_recap():
    global X
    X=read(textPath.get())
    plt.figure(1)
    t = np.arange(0, X[0].stats.npts / X[0].stats.sampling_rate, X[0].stats.delta)
    plt.subplot(111)
    plt.plot(t, X[0].data, 'k')
    plt.ylabel('Raw Data')
    plt.xlabel('Time [s]')
    plt.suptitle(textPath.get())
    plt.show()

    root = Tk()   
    fenetre = RecapCalculs(root)
    root.geometry("500x300+200+0")
    root.mainloop()
    
    
class filterswindow:
    '''
    Interface graphique recapitulant les caracteristique du sismogramme
    presentant les options de filtrage et de calculs du noyau de sensibilite
    '''
    def __init__(self,racine):
        self.canvas = Canvas(racine, borderwidth=1, background="#ffffff")
        self.frame = Frame(self.canvas, background="#ffffff")
        self.vsb = Scrollbar(racine, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)
        self.vsb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)
        self.canvas.create_window((4,4), window=self.frame, anchor="nw", tags="self.frame")

        self.frame.bind("<Configure>", self.OnFrameConfigure)

        self.data()
            
    def data(self):
        global filterVar
        filterVar = 1
        global text6a
        text6a = "1"
        global text6c1
        text6c1 = StringVar()
        global text6c2
        text6c2 = StringVar()
        global text6c3
        text6c3 = StringVar()
           
        Label(self.frame, text="Option Filter").grid(row=0)
        Label(self.frame, text="\n").grid(row=1)
        
        Label(self.frame, text="lowest frequency ?").grid(row=4)
        e1 = Entry(self.frame, textvariable=text6c1)
        e1.grid(row=5)
        
        Label(self.frame, text="highest frequency ?").grid(row=20)
        e2 = Entry(self.frame, textvariable=text6c2)
        e2.grid(row=21)
        
        Label(self.frame, text="number of poles ?").grid(row=22)
        e3 = Entry(self.frame, textvariable=text6c3)
        e3.grid(row=23)
                    
        Button(self.frame, text="continue", command=self.quitter).grid(row=24)
                  
    def quitter(self):
        global racine
        racine.destroy()
        afficheSismoFiltre(textPath.get(), float(text6c1.get()), float(text6c2.get()), float(text6c3.get()))
            
    def OnFrameConfigure(self, event):
        '''Reset the scroll region to encompass the inner frame'''
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        
def optionFilter():
    global racine
    
    racine = Tk()   
    fenetre = filterswindow(racine)
    racine.geometry("500x300+200+0")
    racine.mainloop()    
    
    
def afficheSismoFiltre(Path, frequenceMin, frequenceMax, nbPoles):
    
    # Read the seismogram
    st = read(Path)

    # There is only one trace in the Stream object, let's work on that trace...
    tr = st[0]

    # Filtering with a lowpass on a copy of the original Trace
    tr_filt = tr.copy()
    tr_filt.filter('bandpass', freqmin = frequenceMin, freqmax = frequenceMax, corners = nbPoles, zerophase=True)

    # Now let's plot the raw and filtered data...
    t = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)
    plt.subplot(211)
    plt.plot(t, tr.data, 'k')
    plt.ylabel('Raw Data')
    plt.subplot(212)
    plt.plot(t, tr_filt.data, 'k')
    plt.ylabel('filtered Data')
    plt.xlabel('Time [s]')
    plt.suptitle(tr.stats.starttime)
    plt.show()



class LineBuilder :
    '''
    Fonction permettant de creer des lignes verticales sur un sismo et
    d'enregistrer les coordonees des endroits pointes sur la courbe
    ''' 
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        line.figure.canvas.mpl_connect('button_press_event', self)
        
    def __call__(self, event):
        print 'click', event
        if event.inaxes!=self.line.axes: return
                
        self.xs = [event.xdata, event.xdata]
        self.ys = [ymin, ymax]
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()
        return True
        
def onclick(event): 
    if enregx == []:
        enregx.append(event.xdata)
        enregy.append(X[0].data[event.xdata])
        vlines(event.xdata, ymin, ymax, color='k', linestyles='dashed')
      
def Pointage():
    global enregx
    global enregy
    enregx=[]
    enregy=[]
    
    
        
    global X
    X=read(textPath.get())   
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('click on graphic to store data')
    line, = ax.plot(np.arange(0, len(X[0])) , X[0].data, 'k')
    plt.xlim(0, len(X[0]))
    global ymin
    global ymax
    ymin, ymax = plt.ylim(X[0].data.min()*1.1, X[0].data.max()*1.1)
    
    fig.canvas.mpl_connect('button_press_event', onclick)
        
    line, = ax.plot([0], [0])
    linebuilder = LineBuilder(line)
    plt.show()
    print enregx    # Abscisse des points selectionnes 
    print enregy    # Ordonnees lues sur le sismo
    
    return enregx[0], enregy[0]
       
def main():
    
    '''
    Premiere partie : l'utilisateur entre le chemin du sismo. On affiche
    alors les caracteristiques du sismo
    '''
   
    global root
    root = Tk()
    fenetrePath = Path(root)
    root.geometry("800x600")
    root.mainloop() 
    
    print filterVar.get()
    if filterVar.get() == 1 :
        global racine
    
        racine = Tk()   
        fenetre = filterswindow(racine)
        racine.geometry("500x300+200+0")
        racine.mainloop()  
        text6c = text6c1.get() + ", " + text6c2.get() + ", " + text6c3.get() + "\n"
        
    global X
    X=read(textPath.get())
    splited_path = textPath.get().split('.')
    
    seismicParaVar = [alphaVar.get(),betaVar.get(),allVar.get()] 
    print seismicParaVar   
    
    conversionPara = "0"
    if text8_1 == 1: conversionPara = "1"
     
    
    with open('donnees.inf', 'w') as fichier:
        fichier.write('# 0a. Green function database information file (for a certain depth only for the instance)\n')  
        fichier.write(text0a.get() + '\n')
        fichier.write('# 0b. output directory (parentdir)\n')
        fichier.write(text0b.get() + '\n')
        fichier.write('# 1a. event name\n')
        fichier.write(splited_path[1] + '\n')
        fichier.write('# 1b. event latitude, longitude, depth (however, interpolation for depths won\'t be performed)\n')
        fichier.write(str(X[0].stats.sac.evla) + ", ")
        fichier.write(str(X[0].stats.sac.evlo) + ", ")
        fichier.write(str(X[0].stats.sac.evdp) + '\n')
        fichier.write('# 1c. Mrr, Mtt, Mpp, Mrt, Mrp, Mtp\n')
        fichier.write('1.0 1.0 0.0 1.0 0.0 1.0\n')
        fichier.write('# 2a. station name\n')
        fichier.write(splited_path[0] + '\n')
        fichier.write('# 2b. station latitude, longitude\n')
        fichier.write(str(X[0].stats.sac.stla) + ", ")
        fichier.write(str(X[0].stats.sac.stlo) + "\n")
        fichier.write('# 3. phase name\n')
        fichier.write(text3.get() + '\n')
        fichier.write('# 4. component (Z,R,T)\n')
        fichier.write(splited_path[2] + '\n')
        fichier.write('# 5. seismic parameter (alpha, beta, or all for this version)\n#    if you choose "test" the program will only give you the synthetic\n#          (fort.13 in your directory too)\n')
        if seismicParaVar[0]==1: fichier.write('alpha\n')
        elif seismicParaVar[1]==1 : fichier.write('beta\n')
        else : fichier.write('all\n')
        fichier.write('# 6a. Butterworth filter (if 1 on; if 0 off)\n')
        fichier.write(text6a +"\n")
        fichier.write('# 6b. filter name (mandatory even if 6a=0 )\n')
        fichier.write(text6b.get() + "\n")
        fichier.write("# 6c. if butterworth = 1; lowest freq., highest freq., number of poles\n#     if butterworth = 0; just comment out those parameters (subroutine won\'t read them)\n")
        if filterVar == 1: 
            fichier.write(text6c)
        else : fichier.write('#\n')
        fichier.write('# 7. time window t1, t2, t3, t4 \n#  (if t1=t2 and t3=t4, fwin(:) will be rectangular)\n#  (normally taper functions are sine functions)\n')
        fichier.write(str(t1x) + ", ")
        fichier.write(str(t2x) + "\n")
        #fichier.write(t3x + ", ")
        #fichier.write(t4x + "\n")
        fichier.write('# 8. itranslat (1 if you convert geodetic latitude to geocentric latitude)\n')
        fichier.write(conversionPara + '\n')
        fichier.write('#\n#\n#  Below are minor parameters for kernel calculations\n#                         (i.e. you can leave them as they are to start with)\n#\n# Aa. SINC interpolation window (ipdistance deg) (it works well with 10-20 degrees)\n')
        fichier.write('10.d0\n')
        fichier.write('# Ab. reducing slowness for interplation (c_red_reci s/deg) (if 0.d0 we do not perform slowness reduction)\n')
        fichier.write('0.d0\n')
        fichier.write('# Ba. fast FFT (if 1 on; if 0 off)\n#   you can re-define imin and imax for FFT of Green functions\n#   thence you can avoid reading frequencies for which you don\'t have to account.\n#\n')
        fichier.write('0\n')
        fichier.write('# Bb. if fast FFT = 1; lowest i(freq), highest i(freq) (note that freq = i(freq)/tlen)\n#     if fast FFT = 0; just comment out those parameters (subroutine won\'t read them)\n')
        fichier.write('#0  256 \n')
        fichier.write('# Ca. gridding and extent in R(longitudinal) direction (dph, ph1)\n')
        fichier.write('2.5d-1 5.d0\n')
        fichier.write('# Cb. gridding and extent in T(transverse) direction (dth, thw)\n')
        fichier.write('2.5d-1 5.d0\n')
        fichier.write('# Cd. gridding in radius (rmin, rmax, deltar : should correspond to grids in catalogue)\n# if you put 0.d0 0.d0 0.d0 then the program will take the grids in catalogue\n')
        fichier.write('0.d0 0.d0 0.d0\n')
        fichier.write('# Da. time window (start, end in sec)\n')
        fichier.write('0.d0 1.3d3\n')
        fichier.write('# Db. sampling Hz\n')
        fichier.write('2.d0\n')
        fichier.write('# Ea. ignoring criteria (calculrapide: we ignore the values below; if 0.d0 we don\'t use the algo)\n#         (in Fuji et al. 2012b, we chose 5.d-3)\n')
        fichier.write('0.d0\n')
        fichier.write('#\n# Eb. number of kernel types for the ignoring scheme (if Ea. = 0.d0, just comment out all)\n')
        fichier.write('1\n')
        fichier.write('# Ec. kernel type for ignoring scheme (if 0 we calculate for the envelop) note them vertically\n')
        fichier.write('0\n')
        fichier.write('# F. PSV/SH (PSV only = 2; SH only = 1; PSV+SH =3)\n')
        fichier.write('3\n')
        fichier.write('# don\'t forget write \'end\' at the end\n')
        fichier.write('end\n')
        
    
    newPath = textPath.get() + ".inf"
    print("enregistre sous :%s" % newPath)
        
    try:
        subprocess.call('cp donnees.inf %s' % newPath, shell = True)
    except (Error, OSError), e:
        print "Attempt to copy failed: %s" % e
        
    subprocess.call('rm donnees.inf', shell = True)
    '''
    Deuxieme partie : on affiche le sismo. L'utilisateur peut pointer
    sur le sismo et on enregistrer ces coordonnees sur la courbe.
    '''

if __name__ == "__main__":
    main()
