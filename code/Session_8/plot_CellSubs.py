from pyMCDS import pyMCDS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys

  
def plot_1subs(mcds,n,figure,gs,filenameOut,fileFig,SavePNG):
    cx = mcds.data['discrete_cells']['position_x'];
    cy = mcds.data['discrete_cells']['position_y'];
    
    # First figure
    cellFig   = plt.imread(fileFig)
    ax0 = figure.add_subplot(gs[:,0])
    ax0.axis('off')
    ax0.imshow(cellFig)
    
    # Second figure
    ax1 = figure.add_subplot(gs[0,1])
    subs1 = mcds.get_concentrations( mcds.get_substrate_names()[0] );
    X1,Y1 = mcds.get_2D_mesh();
    CS0 = ax1.contourf(X1,Y1,subs1[:,:,0],cmap='viridis');
    ax1.set_aspect('equal')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = figure.colorbar(CS0,cax=cax)
    ax1.set_title(mcds.get_substrate_names()[0])
    
  
    if (SavePNG):
      figure.savefig(filenameOut)
      #plt.savefig(filenameOut)
    else:
      plt.draw()
      plt.waitforbuttonpress(0) # this will wait for indefinite time
      #plt.pause(0.2)
    plt.clf()

def plot_2subs(mcds,n,figure,gs,filenameOut,fileFig,SavePNG):
    cx = mcds.data['discrete_cells']['position_x'];
    cy = mcds.data['discrete_cells']['position_y'];
    
    # First figure
    cellFig   = plt.imread(fileFig)
    ax0 = figure.add_subplot(gs[:,0])
    ax0.axis('off')
    ax0.imshow(cellFig)
    
    # Second figure
    ax1 = figure.add_subplot(gs[0,1])
    subs1 = mcds.get_concentrations( mcds.get_substrate_names()[0] );
    X1,Y1 = mcds.get_2D_mesh();
    CS0 = ax1.contourf(X1,Y1,subs1[:,:,0],cmap='viridis');
    ax1.set_aspect('equal')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = figure.colorbar(CS0,cax=cax)
    ax1.set_title(mcds.get_substrate_names()[0])
    
    # Third figure
    ax2 = figure.add_subplot(gs[1,1])
    subs2 = mcds.get_concentrations( mcds.get_substrate_names()[1] );
    CS1 = ax2.contourf(X1,Y1,subs2[:,:,0],cmap='viridis');
    ax2.set_aspect('equal')
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = figure.colorbar(CS1,cax=cax)
    ax2.set_title(mcds.get_substrate_names()[1])
    

    if (SavePNG):
      figure.savefig(filenameOut)
      #plt.savefig(filenameOut)
    else:
      plt.draw()
      plt.waitforbuttonpress(0) # this will wait for indefinite time
      #plt.pause(0.2)
    plt.clf()
    
def Example_1subs(initial_index,last_index,folder,SavePNG,func=plot_1subs):
    figure = plt.figure(figsize=(10,8))
    gs = gridspec.GridSpec(nrows=1, ncols=2)
    for n in range( initial_index,last_index+1 ):
        filenameOut=folder+'/output'+"%08i"%n+'.png'
        filename= "output%08i"%n+'.xml'
        fileFig = folder+'/snapshot'+"%08i"%n+'.jpg'
        mcds=pyMCDS(filename,folder)
        func(mcds,n,figure,gs,filenameOut,fileFig,SavePNG)

def Example_2subs(initial_index,last_index,folder,SavePNG,func=plot_2subs):
    figure = plt.figure(figsize=(10,8))
    gs = gridspec.GridSpec(nrows=2, ncols=2, width_ratios=[2,1], height_ratios=[1,1])
    for n in range( initial_index,last_index+1 ):
        filenameOut=folder+'/output'+"%08i"%n+'.png'
        filename= "output%08i"%n+'.xml'
        fileFig = folder+'/snapshot'+"%08i"%n+'.jpg'
        mcds=pyMCDS(filename,folder)
        func(mcds,n,figure,gs,filenameOut,fileFig,SavePNG)        


if __name__ == '__main__':
    if (len(sys.argv) != 5):
      print("Please provide 4 args: InitialTime LastTime, SavePNG (bool), and folder")
      sys.exit(1)
    initial_index = int(sys.argv[1]);
    last_index = int(sys.argv[2]);
    SavePNG = int(sys.argv[3])
    folder = sys.argv[4]
    
    mcds=pyMCDS("initial.xml",folder)
    if (len(mcds.get_substrate_names()) == 1):
        Example_1subs(initial_index,last_index,folder,SavePNG)
    if (len(mcds.get_substrate_names()) == 2):
        Example_2subs(initial_index,last_index,folder,SavePNG)
    

