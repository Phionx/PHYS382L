#################################################################################
#
# This file contains the functions most likely to be modified by all students.
#
# See the lab handout for more details.
#
# First method:
#  set_input:          Define allowed input parameters and defaults
#
# There are five generators, each named 'gen_<name>' that provide 
#   five separate inputs for the main loop that calculated values at each 
#   final temperature. They are:
#  
#  gen_T: Tempature at each step
#  gen_B: Magnetic field at each step
#  gen_collect_EM:  When to sample the average E and M values of each site
#  gen_collect_SC:  When to sample the spin correlation
#  gen_prog_update: This is purely cosmetic and constrolls update outputs to screen
#                   If in doubt, it's ok to just always yield False
#
# The remaining functions are:
#  run_ising_lattice:  Run the simulation for each step
#  make_T_array:       Make the default array of temperatures on which
#                      to run the simulation. Default is evenly spaced from
#                      T_min to T_max with T_spacing
#
#################################################################################

import sys
import time
import numpy as np
import math
#------------ new code ---------
import random
#-----------end new code -------
import matplotlib.pyplot as plt
from sys import exit, argv
from library_interface import IsingLattice
from fn_ising import *
import curses

def set_input(cmd_line_args, parse_cmd_line_args):
    """Parse command-line parameters into an input dictionary.
    Provide default values.
    Add new parameters as you wish.

    note: This fills a Python dictionary called 'inp' which is passed
          around throughout the Python part of the program.

          Feel free to add your own input parameters.
    
    input:   sys.argv
             use syntax of keyword:value on command line
    return:  dict[key] = value
    note:    Any value which can be turned into a number will be a float 
             if it has a '.', otherwise it will be an int.

    """

    inp = dict()

    inp['T_min']      = 1.2    # minimum temperature
    inp['T_max']      = 3.0    # maximum temperature
    inp['T_spacing']  = 0.2    # step size from min to max temperature
    inp['N']          = 100     # sqrt(lattice size) (i.e. lattice = N^2 points

    inp['T0_anneal']    = 4.0    # start temperature (arbitrary; feel free to change)
    inp['steps_anneal']    = 10000  # number of lattice steps in simulation
    inp['steps_burnin']    = 2000   # optional parameter, used as naive default
    inp['EM_samples']      = 5000
    inp['EM_sample_spacing'] = 1

    inp['SC_samples']      = 100
    inp['SC_algorithm']    = 0 # 0 for legacy version, 1 for <ab>-<a><b>
    # default to spreading out SC samples throughout the EM_sample steps

    inp['B']                 = 0.0    # magnetic field strength
    inp['r_flip']            = 0.1    # ratio of sites examined to flip in each step
    inp['dir_out']           = 'data' # output directory for file output
    inp['plots']             = False  # whether or not plots are generated
                                 
    inp['print_inp']         = False  # optional flag
    inp['use_cpp']           = True   # use 1 for True and 0 for False

    inp['date_output']       = False
    inp['file_prefix']       = ''
    inp['multiprocess']      = False
    # inp['skip_prog_print']   = False
    inp['curses']            = True # use the terminal cursor output
    inp['draw_lattice']     = False # is using ncurses, print the lattice at the
                                     # end of each time step
    inp['num_bins']         = 10
    inp['k_b']              = 1

    #function call from support_methods
    parse_cmd_line_args(inp, cmd_line_args)

    if inp['print_inp']:
        print('Printed list of input keys:')
        for key in sorted(inp.keys()):
            print('%-20s'%key,' ',inp[key])
    return inp

def gen_T(inp, T_final):
    '''Yield values of temperature.
    Default implementation: 
    (1) Start at T0_anneal and linearly drop to T_final in steps anneal
    (2) Yield T_final for all remaining steps.'''

    for T in np.linspace(inp['T0_anneal'],T_final,inp['steps_anneal']):
        yield T

    for T in range(inp['steps_burnin'] +
                   inp['EM_samples']*inp['EM_sample_spacing']):
        yield T_final

def gen_B(inp, T_final=None):
    '''Yield magnetic field values.
    Default implementation: Always yield inp['B'].'''
    while True:
        yield inp['B']

def gen_collect_EM(inp, T_final=None):
    '''Yield when to sample E and M values.
    Default implementation: Don't sample during steps_anneal
      and steps_burning. Then sample EM_samples separated by
      EM_sample_spacing-1 between each sample.'''
    for _ in range(inp['steps_anneal']+inp['steps_burnin']):
        yield False
    
    for _ in range(inp['EM_samples']):
        for __ in range(inp['EM_sample_spacing']-1):
            yield False
        yield True

def gen_collect_SC(inp, T_final=None):
    '''Yield when to collect Spin Correlation values.
    Default implementation: Collect EM_samples evenly spaced 
             throughout where E and M values are sampled.'''
    for _ in range(inp['steps_anneal']+inp['steps_burnin']):
        yield False
    
    space_corr = int( (inp['EM_samples']*inp['EM_sample_spacing'])
                      /inp['SC_samples'])
    if space_corr == 0:
        print('Fatal error.',
         'Collecting more SC_samples than EM_samples*EM_sample_spacing.')
        exit(1)
    if space_corr == 1:
        space_corr = 2

    for _ in range(inp['SC_samples']):
        for __ in range(space_corr-1):
            yield False
        yield True

    while True:
        yield False

def gen_prog_update(inp):
    '''Yield when to update the screen. 
    Do not use with multiprocess.'''

    #if using multiprocess, always return False
    if inp['multiprocess']:
        while True:
            yield False

    n_total = (inp['steps_anneal']+inp['steps_burnin']
              +inp['EM_samples']*inp['EM_sample_spacing'])
    n_spacing = int( n_total/100);
    # print (f'\n\n {n_spacing} \n')
    if n_spacing == 0:
        n_spacing = 1
    n_complete = 0

    while True:
        for _ in range(n_spacing-1):
            n_complete += 1
            yield False
            
        n_complete += 1
        yield (n_complete, n_total)

def run_ising_lattice(inp, T_final, updates=True):
    '''Run a 2-D Ising model on a lattice. Use the above
    five generators for loop values.
    Return mean and std values for Energy, Magnetization, and Spin Correlations
    '''
    if updates:
        inp['print'].start_Tprogress(T_final)

    lattice = IsingLattice(inp)
    if inp['SC_algorithm'] == 0:
        get_SC = lattice.get_SC_v0;
    else:
        get_SC = lattice.get_SC_v1;

    try:
        E_collect = []
        M_collect = []
        SC_collect = []
        
#---------- new code ----------------------- 
        #initialize bins
        E_std_devs = []
        M_std_devs = []
        k = math.floor(inp['EM_samples']/inp['num_bins']) #number of analysis steps per bin
        k_b = inp['k_b'] #pull botzmann's constant from input
#------- end new code ------------------------
        
        for T, B, sample_EM, sample_SC, prog_update in zip(
            gen_T(inp, T_final),
            gen_B(inp),
            gen_collect_EM(inp, T_final),
            gen_collect_SC(inp, T_final),
            gen_prog_update(inp)
        ):
            lattice.step(T,B)

            if sample_EM:
                E_collect.append(lattice.get_E())
                M_collect.append(lattice.get_M())
            if sample_SC:
                SC_collect.append( get_SC() )
            if updates and prog_update:
                inp['print'].update_Tprogress(prog_update)

        if updates:
            inp['print'].update_Tprogress(prog_update,finished=True)
        
        if inp['draw_lattice'] and not inp['multiprocess']:
            inp['print'].draw_lattice(lattice, T=T_final)

        lattice.free_memory()
        SC_data = np.array(SC_collect)
        
        shuffled_E = E_collect
        shuffled_M = M_collect
        random.shuffle(shuffled_E)
        random.shuffle(shuffled_M)
        
#---------- new code -----------------------
        while len(shuffled_E) > k:
            E_std_devs.append(np.std(shuffled_E[:k]))
            shuffled_E = shuffled_E[k:]
            M_std_devs.append(np.std(shuffled_M[:k]))
            shuffled_M = shuffled_M[k:]
#-------- end new code ----------------------
            
                
#---------- new code -----------------------
        #calculate specific heat and susceptibility and std devs
        spec_heat = (1/(k_b*(T_final**2)))*(np.mean(E_std_devs)**2)
        std_spec_heat = (((1/(k_b*(T_final**2)))*2*np.mean(E_std_devs))**2 * (np.std(E_std_devs)**2))**(.5)
        chi = (1/(k_b * T_final))*(np.mean(M_std_devs)**2)
        std_chi = (((1/(k_b*T_final))*2*np.mean(M_std_devs))**2 * np.std(M_std_devs)**2)**(.5)
#-------- end new code ----------------------
        
        return (
            (T_final, len(E_collect), np.mean(E_collect), np.std(E_collect), np.mean(M_collect), np.std (M_collect), spec_heat, std_spec_heat, chi, std_chi),
            (T_final, np.shape(SC_data)[0], np.array(SC_data.mean(0)), np.array(SC_data.std (0)))
        )

    except KeyboardInterrupt:
        try:
            lattice.free_memory()
        except:
            pass
        print("\n\nProgram terminated by keyboard. Good Bye!")
        sys.exit()

def make_T_array(inp):
    if inp['T_max'] <= inp['T_min']:
        return [inp['T_min'],]
    else:
        # Add a slight offset to the upper bound so T_max won't be excluded
        # when T_max-T_min is evenly divisible by T_spacing. 
        return np.arange(inp['T_min'], inp['T_max']+1E-5, inp['T_spacing'])
