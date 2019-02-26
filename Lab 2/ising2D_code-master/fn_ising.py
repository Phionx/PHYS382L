import os
import time
import csv
import numpy as np
import logging
import curses
import matplotlib.pyplot as plt
'''This is a file with methods (and classes) called by mod_ising.py.
You are welcome to modify these. However, it is not anticipated that you
should need to do so.
'''

# colors for the ncurses screen
'''Up ()  down []'''
UP_FG = 16
UP_BG = 205
DN_FG = 16
DN_BG = 63
i_UP  = 10
i_DN  = 11


def get_filenames(inp): #make data folder if doesn't exist, then specify filename
    '''Generate the output file names for the EM (energy and megnetism) and SC (spin correlation) files'''
    try:
        dir_out = inp['dir_out']
        prefix  = inp['file_prefix']
        if inp['date_output']:
            dir_out += str(time.strftime("_%Y%m%d-%H%M%S"))

        if not os.path.isdir(dir_out):
            os.makedirs(dir_out)

        # file name = [file_prefix]##_EM_v#.csv if only one temperature (example: runA_4.20_EM_v0.csv)
        #             [file_prefix]##T##_EM_v#.csv if there are two temperatures (example: runA_4.2T5.3_EM_v0.csv) 
        # the other file name is identical, but with "SC" (for spin correlation)) instead of EM
        if inp['T_max'] <= inp['T_min']:
            t_name = '%.2f'%inp['T_min']
        else:
            t_name = '%.2fT%.2f'%(inp['T_min'],inp['T_max'])

        v = 0
        while True:
            EM_file = f'{dir_out}/{prefix}{t_name}_EM_v{v}.csv'
            SC_file = f'{dir_out}/{prefix}{t_name}_SC_v{v}.csv'
            if not (os.path.isfile(EM_file) or os.path.isfile(SC_file)):
                break
            v += 1
        return EM_file, SC_file

    except:
        print ('fatal: Failed to make output file names')
        sys.exit()

def parse_cmd_line_args(inp, cmd_line_args):
    '''Parse arguments from the command line into input dictionary inp'''
    for x in cmd_line_args[1:]:
        if ':' in x:
            try:
                key, val = x.split(':')
                if key not in inp:
                    print(f'fatal error: invalid command-line input "{key}"')
                    print('  If you want to add a new input, add it in\n'
                          '  the "set_input" function in mod_ising.py')
                    raise KeyError
                try:
                    if '.' in val:
                        inp[key] = float(val)
                        print('%-20s'%('inp["%s"]'%key),'set to float  ',inp[key])
                    elif val.lower() == 'false' or val.lower() == 'f':
                        inp[key] = False
                    elif val.lower() == 'true' or val.lower() == 't':
                        inp[key] = True
                    else:
                        inp[key] = int(val)
                        print('%-20s'%('inp["%s"]'%key),'set to int    ',inp[key])
                except:
                    inp[key] = val
                    print('%-20s'%('inp["%s"]'%key),'set to string ',inp[key])
            except KeyError:
                exit(2)
            except:
                print('warning: input "%s" not added to arguments'%x)
        else:
            print(f'fatal error: invalid command-line format "{x}"')
            print( '             proper format is key:value key:value etc...')
            exit(2)

class checkered_nums:
    '''Class to write increasing digits (always two columns limit) to column and row.
       Print in checkered colors: black-on-white then white-on-black'''
    def __init__(self, stdscr_in):
        self.stdscr = stdscr_in
        self.i = 1
    def print(self, row, col):
        if col+2 > curses.COLS or row >= curses.LINES:
            return
        num = self.i
        if num > 100:
            num = num % 100
        num = '%2i'%num
        # self.stdscr.addstr(col,row,'%i'%num, curses.color_pair(5))
        self.stdscr.addstr(row,col,num,curses.color_pair(5+self.i%2))
        self.i += 1

class printer:
    '''class to print to the screen'''
    def __init__(self, inp, stdscr=None):
        ''' Use stdscr (curses terminal screen) if present.  Otherwise use stdout ''' 
        self.start_time = time.time()
        self.inp = inp
        self.stdscr = stdscr

        if self.stdscr:
            self.col = 0
            self.row = 1
            curses.start_color()
            curses.use_default_colors()
            curses.init_pair(1, curses.COLOR_RED,   254) #curses.COLOR_WHITE)
            curses.init_pair(2, curses.COLOR_WHITE, curses.COLOR_MAGENTA)
            curses.init_pair(3, curses.COLOR_GREEN, 16) #curses.COLOR_BLACK)
            curses.init_pair(4, curses.COLOR_WHITE, curses.COLOR_CYAN)
            curses.init_pair(5, curses.COLOR_WHITE, 245) #16) #curses.COLOR_BLACK)
            curses.init_pair(6, curses.COLOR_BLACK, curses.COLOR_WHITE)
            curses.init_pair(7, curses.COLOR_RED, 16)
            curses.init_pair(8, curses.COLOR_WHITE, 16)
            curses.init_pair(i_UP, UP_FG, UP_BG)
            curses.init_pair(i_DN, DN_FG, DN_BG)

    def start_singlecorr(self):
        '''Write start of program message to screen'''
        if self.stdscr:
            self.stdscr.addstr(0,0, f'{self.inp["N"]}', curses.color_pair(7))
            self.stdscr.addstr('x')
            self.stdscr.addstr(f'{self.inp["N"]}', curses.color_pair(7))
            self.stdscr.addstr(' 2D Ising Lattice: Single Core')
                    # 2D Ising Lattice: Single Core',
                   # curses.color_pair(7) )
            self.stdscr.refresh()
        else:
            print(f'{self.inp["N"]}x{self.inp["N"]} 2D Ising Lattice: Single Core')

    def start_Tprogress(self, Tfinal):
        # get number of finished steps
        '''Write initial progress line to screen for a given temperature'''
        self.Tfinal = Tfinal;
        self.Tstart_time = time.time()
        self.fmt_prog = ('T:{Tfinal:>.3} '
                'steps: {prog[0]:>7}/{prog[1]:>7}, {ratio:3}% '
                   'time: {time_pass} '
                   'est.time-to-go: {est_time}'
                   )
        self.fmt_finish = ('T:{Tfinal:>.2} '
                   'steps: {prog[1]:>7}/{prog[1]:>7}, 100% '
                   'time: {time_pass} '
                   'est.time-to-go: done!'
                   )
        if self.stdscr:
            pass
        else:
            print(f'T:{self.Tfinal:>.2}',end='\r')

    def update_Tprogress(self, prog,finished=False):
        time_pass = time.time() - self.Tstart_time 
        time_pass_str = time.strftime('%H:%M:%S', time.gmtime(time_pass));

        # if finished is True, prog might = False (if the iteration
        # didn't end on an update)
        if prog:
            self.last_prog = prog
        else:
            prog = (self.last_prog[1],self.last_prog[1])

        try:
            ratio = float(prog[0])/prog[1]
        except:
            print(f'\n MY ERROR prog: {prog}')
        est_time  = (1-ratio)*time_pass/ratio

        est_time_str = time.strftime('%H:%M:%S', time.gmtime(est_time));
        steps = list(prog)
        if finished:
            est_time_str = 'done!    '
            steps[0] = steps[1]
        msg = self.fmt_prog.format(prog=steps,ratio=(int(100*ratio)),
                                   time_pass=time_pass_str,
                                   est_time=est_time_str,
                                   Tfinal=self.Tfinal)
        if self.stdscr:
            self.stdscr.addstr(self.row,0, msg)
            self.stdscr.refresh()
            if finished:
                self.row_plus1()
        else:
            if finished:
                print(msg)
            else:
                print(msg,end='\r')
    def finish_Tprogress(self):
        pass

    def row_plus1(self):
        '''progress the column count. Loop if necessary'''
        self.row += 1
        if self.row > curses.LINES-3:
            self.stdscr.addstr(0,0,'... output looping from bottom of screen ...')
            self.stdscr.clrtoeol()
            self.row = 1

    #             print(string)
    def print_done(self):
        time_pass_str = time.strftime('%H:%M:%S', time.gmtime(time.time()-self.start_time));
        if self.stdscr:
            self.stdscr.clrtobot()
            self.stdscr.refresh()
            self.stdscr.addstr(self.row,0,f'Total run time: {time_pass_str} seconds.')
            self.row_plus1()
            self.stdscr.addstr(self.row, 0, '  Finished ising model program. Press any key to continue.  ',
                    curses.color_pair(1))
            self.stdscr.clrtoeol()
            # curses.flash()
            self.stdscr.getkey()
        else:
            print('\n---- Program Finished ----')
            print(f'Total run time: {time_pass_str} seconds.')

    def draw_lattice(self, lattice, step=0, nstep=0, phase='',T='', B=''):
        '''
        Print the lattice.
        If not in curses windows, use lattice.print_spins().
        If in curses window, print as:
       0         1         2         3         4         5         6         7
       01234567890123456789012345678901234567890123456789012345678901234567890123456789
        200x200 Step   12000/100000   [Anneal/Burning/Analysis]
        T 1.34    B 9.13   |E| 12.34   |M| -0.123  
                             1 2 3 4 5 6 7
        0 1 0 1 0 1 0  or  1()[]()[]()[]()
        1 1 0 0 0 0 0      2()()()()()()()
        1 0 1 0 0 1 1      3()[][][][]()()
        1 0 1 0 1 0 1      4()[]()[]()[]()
        0 0 0 1 0 1 0      5[][][]()[]()[]
        1 0 1 0 0 1 0      6()[]()[][]()[]
        0 1 0 1 0 0 1      7[]()[]()[][]()
        '''
        if not self.stdscr:
            print('Not using n-curses screen. Will not print 2D Map')
            return

        # check for screen large enough
        N = self.inp['N']
        # if curses.COLS < 2*(N+1):
        #     self.stdscr.addstr(1, 1, 'Screen not wide enough to show 2D-map',
        #             curses.color_pair(1))
        #     return

        self.stdscr.clear()
        self.stdscr.addstr(0,0,curses.COLS*' ',curses.color_pair(3))
        self.stdscr.addstr(1,0,curses.COLS*' ',curses.color_pair(3))
        self.stdscr.addstr(0,1,'%ix%i'%(N,N), curses.color_pair(3))
        if step and nstep:
            self.stdscr.addstr(0,9,'Step %i/%i'%(step,nstep), curses.color_pair(3))
        if phase:
            self.stdscr.addstr(0,31,phase, curses.color_pair(3))
        if T:
            self.stdscr.addstr(1,1,'T %.2f'%T, curses.color_pair(3))
        if B:
            self.stdscr.addstr(1,11,'B %.2f'%B, curses.color_pair(3))

        self.stdscr.addstr(1,20,'|E| %.2f'%lattice.get_E(),curses.color_pair(3))
        self.stdscr.addstr(1,32,'|M| %.2f'%lattice.get_M(),curses.color_pair(3))
        values=lattice.get_numpy_spin_matrix()
        self.stdscr.addstr(3,0,'%i %i'%(np.min(curses.LINES-2),N))

        # plot column indices
        row_labels = checkered_nums(self.stdscr)
        col_labels = checkered_nums(self.stdscr)
        flip = 0
        for row in range(np.min([curses.LINES-3,N])):
            flip += 1
            # if flip%2:
            #     curses.init_pair(2,curses.COLOR_RED, curses.COLOR_WHITE)
            # else:
            #     curses.init_pair(2,curses.COLOR_BLUE, curses.COLOR_BLACK)
            row_labels.print(row+3,0)
            for col in range(np.min([int(curses.COLS/2)-2, N])):
                if row == 1:
                    col_labels.print(2, col*2+2)
                if values[row,col]==1:
                    self.stdscr.addstr(row+3,col*2+2,'()',curses.color_pair(i_UP))
                else:
                    self.stdscr.addstr(row+3,col*2+2,'[]',curses.color_pair(i_DN))

        self.stdscr.addstr(curses.LINES-1, 0, 'press a key to continue',
                curses.color_pair(1))
        # self.stdscr.addstr(5,0,'%s'%pairs)
        self.stdscr.refresh()
        self.stdscr.getkey()  

def plot_graphs(data): #T, n_sample, E_mean,E_std,M_mean,M_std): #plot graphs at end
    dat = np.array(data)
    print("0: " + str(dat[0,0]) + "\n")
    print("1: " + str(dat[0,1]) + "\n")
    print("2: " + str(dat[0,2]) + "\n")
    print("3: " + str(dat[0,3]) + "\n")
    print("4: " + str(dat[0,4]) + "\n")
    plt.figure(1)
    # plt.ylim(0,1)
    plt.errorbar(dat[:,0], dat[:,2], yerr=dat[:,3], fmt='o')
    plt.xlabel('Temperature')
    plt.ylabel('Average Site Energy')
    plt.figure(2)
    plt.errorbar(dat[:,0], np.absolute(dat[:,4]), yerr=dat[:,5], uplims=True, lolims=True,fmt='o')
    plt.xlabel('Temperature')
    plt.ylabel('Aveage Site Magnetization')
    plt.show()

def print_results(inp, EM_data, SC_data):
    data_filename, corr_filename = get_filenames(inp)
    with open(data_filename,'w') as f_out:
        writer = csv.writer(f_out, delimiter=',', lineterminator='\n')

        # first echo all the 'inp' dictionary, with a few 'important' settings first
        # followed by the remaining settings
        first_keys = 'N T_min T_max T_spacing r_flip EM_samples SC_samples'.split()
        first_vals = [inp[K] for K in first_keys]
        rest_keys = [K for K in inp.keys() if K not in first_keys]
        rest_keys.remove('print') # print is a method passed by convenience. Don't print it.
        rest_keys.sort()
        rest_vals = [inp[K] for K in rest_keys]
        writer.writerow(first_keys+rest_keys)
        writer.writerow(first_vals+rest_vals)

        # write an empty row
        writer.writerow([' '])
        writer.writerow('Temp n_samples E_mean E_std M_mean M_std'.split())
        for entry in EM_data:
            writer.writerow(entry)

    with open(corr_filename,'w') as f_out:
        writer = csv.writer(f_out, delimiter=',', lineterminator='\n')
        first_keys = 'N T_min T_max T_spacing r_flip EM_samples SC_samples'.split()
        first_vals = [inp[K] for K in first_keys]
        rest_keys = [K for K in inp.keys() if K not in first_keys]
        rest_keys.remove('print')
        rest_keys.sort()
        rest_vals = [inp[K] for K in rest_keys]
        writer.writerow(first_keys+rest_keys)
        writer.writerow(first_vals+rest_vals)

        writer.writerow([' '])
        writer.writerow(['Temp','n_samples']
                       +['mean_d%i'%i for i in range(1,len(SC_data[0][2])+1)]
                       +['std_d%i'%i  for i in range(1,len(SC_data[0][2])+1)]
        )
        for entry in SC_data:
            writer.writerow(np.hstack(entry))

