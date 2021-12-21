import matplotlib.pyplot as plt
plt.style.use('dark_background')

from numpy import loadtxt, linspace

def Plot_OrbElem(t, Orbital_elements, Name='test particle', Color='cyan'):
    '''
    ----------------------------------------------------------------
    Plot_OrbElem(t, Orbital_elements):
    Draw the orbital elements over a time interval.
    ----------------------------------------------------------------
    Arguments:
    t:  Time interval.
    Orbital_elements: Numpy array with the orbital elements of each body
        over a time interval, following the next format:
        Orbital_elements=[a,e,i,omega, Omega]
    Name: Body's name.
    Color: Plot's color. Crimson by default.
    ----------------------------------------------------------------
    '''

    plt.rcParams['axes.grid'] = True
    plt.rcParams['axes.formatter.limits']= -3, 3
    fig, axs = plt.subplots(2, 3, figsize=(16,8))

    fig.suptitle('COEs of ' + Name,fontsize=20)

    #Plot Semi-major axis
    axs[0, 0].plot(t, Orbital_elements[:,0], color=Color)
    axs[0, 0].set_title('$a$')
    axs[0, 0].set_ylabel('[au]')
    axs[0, 0].set_xlabel('$t$ [yr]')

    #Plot Eccentricity
    axs[0, 1].plot(t, Orbital_elements[:,1], color=Color)
    axs[0, 1].set_title('$e$')
    axs[0, 1].set_xlabel('$t$ [yr]')

    #Plot Inclination
    axs[0, 2].plot(t, Orbital_elements[:,2], color=Color)
    axs[0, 2].set_title('$\iota$')
    axs[0, 2].set_ylabel('[rads]')
    axs[0, 2].set_xlabel('$t$ [yr]')

    #Plot argument of periapse
    axs[1, 0].plot(t, Orbital_elements[:,3], color=Color)
    axs[1, 0].set_title('$\omega$')
    axs[1, 0].set_ylabel('[rads]')
    axs[1, 0].set_xlabel('$t$ [yr]')

    #Plot Longitude of ascendind node
    axs[1, 1].plot(t, Orbital_elements[:, 4], color=Color)
    axs[1, 1].set_title('$\Omega$')
    axs[1, 1].set_ylabel('[rads]')
    axs[1, 1].set_xlabel('$t$ [yr]')

    plt.subplots_adjust(wspace=0.4,hspace=0.4)
    plt.setp(axs, xlim=[t[0],t[-1]])
                        
    axs[1, 0].set_position([0.24,0.125,0.228,0.343])
    axs[1, 1].set_position([0.55,0.125,0.228,0.343])
    axs[1, 2].set_visible(False)
    plt.show()

# Body's name
Names = ['Sun', 'Jupiter','(1373)Cincinnati'] # Remember to update the 3rd body!!!
#Body's masses [M_sun]
Masses = [1., 9.54792e-04, 1.0e-16]

# Path to the files
path_OP_B1 = './Files/OP_' + Names[1] + '.txt'
path_OP_B2 = './Files/OP_' + Names[2] + '.txt'
path_CC_B0 = './Files/CC_' + Names[0] + '.txt'
path_CC_B1 = './Files/CC_' + Names[1] + '.txt'
path_CC_B2 = './Files/CC_' + Names[2] + '.txt'

# dt: Discrete time step.
# n: Number of time-iterations in one outer period.
# k: Number of outer periods.
# jump: Jump size to store data in files
dt, n, k, jump =  loadtxt(path_OP_B1, skiprows=1, max_rows=1, unpack=True,
                          dtype={'names': ('dt','n','k','jump'), 'formats': (float, int, int, int)})

# Array to store time information
it = int(n/ jump)
t = linspace(0, dt*n*k, k*it)
skip_rows = 2
# Plots of orbital elements.
Plot_OrbElem(t, loadtxt(path_OP_B1, skiprows=skip_rows), Names[1])
Plot_OrbElem(t, loadtxt(path_OP_B2, skiprows=skip_rows), Names[2])