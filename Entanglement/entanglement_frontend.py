import tkinter as tk
from tkinter import ttk
import matplotlib
matplotlib.use("TkAgg") # Use Tkinter as teh matplotlib back-end
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg #, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import numpy as np
import matplotlib.animation as animation
from matplotlib import style

LARGE_FONT= ("Verdana", 12)
style.use('dark_background')

fig = Figure(figsize=(15,4), dpi=100)
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

def animate(i):
    pullData = open('spdcdata.txt','r').read()
    dataArray = pullData.split('\n')
    det1=[]
    det2=[]
    coinc=[]

    for eachLine in dataArray:
        if len(eachLine)>1:
            d1,d2,c2x = eachLine.split(',')
            det1.append(float(d1))
            det2.append(float(d2))
            coinc.append(float(c2x))    
    
    ax1.clear()
    ax1.plot(det1,'-o',c='crimson',markersize=10)
    ax1.grid(True)    
    ax1.set_ylim(ymin=0)
    
    ax2.clear()
    ax2.plot(det2,'-o',c='cyan',markersize=10)
    ax2.grid(True)    
    ax2.set_ylim(ymin=0)
    
    ax3.clear()
    ax3.plot(coinc,'-o',c='magenta',markersize=10)
    ax3.grid(True)    
    ax3.set_ylim(ymin=0)
    
    if len(det1) is not 0:
        ax1.annotate('Left: '+str(round(det1[-1])),xy=(280, 360), xycoords='figure pixels',size=30)
        ax2.annotate('Right: '+str(round(det2[-1])),xy=(660, 360), xycoords='figure pixels',size=30)
        ax3.annotate('Dual: '+str(round(coinc[-1])),xy=(1090, 360), xycoords='figure pixels',size=30)
    else:
        print('Empty list!!!!')
        
# Inheritence is given in parenthesis
class mainApp(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        
        tk.Tk.wm_title(self, "SPDC Coincidence Count")
        
        container = tk.Frame(self)
        # Nicely pack up components
        container.pack(side = "top",fill = "both",expand = True) 
        container.grid_rowconfigure(0,weight=1)
        container.grid_columnconfigure(0,weight=1)
    
        self.frames = {}
        
        for fr in (StartPage,PlotPage):
            frame = fr(container, self)
            self.frames[fr] = frame
            frame.grid(row=0, column = 0, sticky="nsew") #North South East West
    
        self.show_frame(StartPage)
    
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()
        
class StartPage(tk.Frame):

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        
        label = ttk.Label(self, text="Press \"Start\" to display analysis", font=LARGE_FONT)
        label.pack(pady=10,padx=10)
        
        button1 = ttk.Button(self, text = "Start",command=lambda: controller.show_frame(PlotPage))
        button1.pack()
        
class PlotPage(tk.Frame):
    
    def __init__(self, parent, controller):
        
        ttk.Frame.__init__(self, parent)
        label = tk.Label(self, text="For graphing", font=LARGE_FONT)
        button1 = tk.Button(self, text = "Go Back",command=lambda: controller.show_frame(StartPage))
        button1.pack()        
        
        canvas = FigureCanvasTkAgg(fig, self)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        
        # toolbar = NavigationToolbar2TkAgg(canvas, self)
        # toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
app = mainApp()
ani = animation.FuncAnimation(fig,animate, interval=200)
app.mainloop()