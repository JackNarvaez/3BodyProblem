import matplotlib.pyplot as plt
plt.style.use('dark_background')

from sys import argv
from numpy import loadtxt, linspace

def Plot_Energy(E, Names='test particle', Color='cyan'):
    '''
    ----------------------------------------------------------------
    Plot Energy of the system.
    ----------------------------------------------------------------
    Arguments:
    E: Energy
    ----------------------------------------------------------------
    '''

    plt.rcParams['axes.grid'] = True
    plt.rcParams['axes.formatter.limits']= -3, 3
    fig, axs = plt.subplots(figsize=(16,8))

    fig.suptitle('Energy ' + Names[0] + "-" + Names[1] + "-" + Names[2],fontsize=20)

    #Plot E
    axs.plot(E, color=Color)
    axs.set_ylabel('E')
    axs.set_xlabel('it')
    axs.set_xlim(0, len(E))

    plt.show()

# Body's name
skip2 = int(argv[1])-1
skip3 = int(argv[2])-1
Second_name = str(loadtxt('Data.asc', skiprows=14+skip2, max_rows=1,  usecols=[0], dtype='str')) # Name of third body
Third_name = str(loadtxt('Data.asc', skiprows=14+skip3, max_rows=1,  usecols=[0], dtype='str')) # Name of third body
Names = ['Sun', Second_name, Third_name] # Remember to update the 3rd body!!!

# Path to the files
path_E = './Files/Energy_' + Third_name + '.txt'

E =  loadtxt(path_E, skiprows=3, unpack=True)

# Plots of Energy.
Plot_Energy(E, Names)