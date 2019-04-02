import multiprocessing as mp
from time import sleep
from mod_ising import *

def multi_lattice(inp, T):
    print(f'Starting Temp {T:.2}')
    time_start = time.time()
    rval = run_ising_lattice(inp,T,updates=False)
    time_pass_str = time.strftime('%H:%M:%S', time.gmtime(time.time()-time_start));
    print(f'Finished Temp {T:.2}. Total time: {time_pass_str}')
    return rval

def run_multi_core(inp):
    print("\n2D Ising Model Simulation; multi-core\n")
    T_array = make_T_array(inp)

    pool = mp.Pool(mp.cpu_count())
    jobs = [pool.apply_async(multi_lattice,(inp,T)) for T in T_array]
    results = [x.get() for x in jobs]
    EM_data = [x[0] for x in results]
    SC_data = [x[1] for x in results]
    EM_data.sort()
    SC_data.sort()
    print_results(inp, EM_data, SC_data)

def run_single_core(inp):
    # sequentially run through the desired temperatures and collect the output for each temperature
    EM_data = []
    SC_data = []
    inp['print'].start_singlecorr()
    for T in make_T_array(inp):
        EM_vals, SC_vals = run_ising_lattice(inp, T)
        EM_data.append( EM_vals )
        SC_data.append( SC_vals )

    print_results(inp, EM_data, SC_data)

    if inp['plots']:
        plot_graphs(EM_data)

class run_wrapper:
    '''Wrap execution either in, or out, of a python cursor terminal ("curses")'''
    def __init__(self,inp):
        self.inp = inp
    def run(self, stdscr=None):
        self.inp['print'] = printer(self.inp, stdscr)
        if self.inp['multiprocess']:
            run_multi_core(self.inp)
        if not self.inp['multiprocess']:
            run_single_core(self.inp)
        self.inp['print'].print_done()

if __name__ == "__main__":
    """Main program: run Ising Lattice here"""
    inp = set_input(argv, parse_cmd_line_args)

    # using ncurses with multiprocess will cause the program to fail
    if inp['multiprocess']:
        inp['curses'] = False

    wrap = run_wrapper(inp)
    
    if inp['curses']:
        #this is the call that provides stdscr to run_wrapper.run
        curses.wrapper(wrap.run) 
    else:
        wrap.run()
