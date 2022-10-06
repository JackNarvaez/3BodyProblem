import matplotlib.pyplot as plt
plt.style.use('dark_background')

from sys import argv
from numpy import loadtxt

def Plot_PS(w, e, Name='test particle', Color='cyan'):
    '''
    ----------------------------------------------------------------
    Plot phase space: e vs w.
    ----------------------------------------------------------------
    Arguments:
    e: Eccentricity
    w: The argument of perihelion
    Name: Body's name.
    Color: Plot's color. cyan by default.
    ----------------------------------------------------------------
    '''
        
    plt.rcParams['axes.grid'] = True
    plt.rcParams['axes.formatter.limits']= -3, 3
    fig = plt.figure(figsize=(16,8))

    fig.suptitle('Phase Space of ' + Name,fontsize=20)

    #Plot in Cartesian Coordinates
    axs0 = fig.add_subplot(121)
    axs0.plot(w, e, color=Color)
    axs0.set_title('Cartesian Coordinates')
    axs0.set_ylabel('e')
    axs0.set_xlabel('$\omega$ [rads]')

    #Plot in Polar Coordinates
    axs1 = fig.add_subplot(122, projection='polar')
    axs1.plot(w, e, color=Color)
    axs1.set_title('Polar Coordinates', va='bottom')

    plt.subplots_adjust(wspace=0.4,hspace=0.4)
                        
    plt.show()

# Body's name
skip = int(argv[1])-1
Name = str(loadtxt('Data.asc', skiprows=14+skip, max_rows=1,  usecols=[0], dtype='str')) # Name of third body

# Path to the files
path_OP_B2 = './Files/OP_' + Name + '.txt'

skip_rows = 2

# Plots of orbital elements.
ecc, omega = loadtxt(path_OP_B2, skiprows=skip_rows, usecols=[1, 3], unpack=True, 
                    dtype={'names': ('ecc','omega'), 'formats': (float, float)})
Plot_PS(omega, ecc, Name)