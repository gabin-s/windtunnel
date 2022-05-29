from operator import delitem
import matplotlib.pyplot as plt
import numpy as np

def read_data(fname, nparams=1):
    """
    fname: str
        filename from which to load the data
    nparams: int
        number of parameters columns
    """
    d = np.loadtxt(fname, delimiter=' ', comments='#')

    params = d[:, :nparams].T
    raw    = d[:, nparams:]
    
    means = np.mean(raw, axis=1)
    stds   = np.std(raw, axis=1)

    return (*params, raw.shape[1], means, stds)


def exp_variate_grid(nthreads=4):
    """
    Compare the runtime when changing the number of columns vs. when changing
    the number of rows.

    Description: 
        We vary the number of columns (respectively rows) while keeping
        a constant number of rows (respectively columns) of 1000. 
        The runtime is plotted against the number of columns (respectively rows).

        The number of threads is constant.

    Datafiles:
        variate_grid_{nthreads}t_cols.results
        variate_grid_{nthreads}t_rows.results

    Results:
        Scales better when increasing `nrows` compared to `ncols`: 
            low overhead when adding more than 1000 rows
            vs. runtime linearly tied to the number of columns

    """
    plt.figure()
    for lbl in ['cols', 'rows']:
        n, nruns, t, stds = read_data(f'variate_grid_{nthreads}t_{lbl}.results')
    
        stderr = stds / np.sqrt(nruns)

        plt.errorbar(n, t, yerr=stderr, label=lbl, capsize=2)
        
    plt.legend(title='parameter')
    plt.xlabel('number of columns/rows')
    plt.ylabel(f'Average runtime over {nruns} runs ($s$)')
    plt.tight_layout()

def exp_num_threads(lbl="25Melem"):
    """
    Runtime as a function of the number of threads

    Description: 
        We vary the number of threads and measure the runtime for a fixed grid
        size, without any particles

    Datafiles:
        num_threads_{lbl}.results

    Results:
        The runtime quickly drops as we increase the number of threads up to 
        ~100 threads. After this point it starts increasing.

    """
    plt.figure()
    nthreads, nruns, t, stds = read_data(f'num_threads_{lbl}.results')
    
    stderr = stds / np.sqrt(nruns)

    print(f'exp_num_threads({lbl}): max(t)={max(t)}')

    plt.errorbar(nthreads, t, yerr=stderr, capsize=2, elinewidth=1, ecolor='red')
    plt.xlabel('number of threads')
    plt.ylabel(f'Average runtime over {nruns} runs ($s$)')
    plt.tight_layout()

def exp_num_threads_speedup(lbl="25Melem"):
    """
    Speedup as a function of the number of threads

    Description: 
        We vary the number of threads

    Datafiles:
        num_threads_{lbl}.results

    Results:
        The runtime quickly drops as we increase the number of threads up to 
        ~100 threads. After this point it starts increasing.

    """
    plt.figure()
    d = np.loadtxt(f'num_threads_{lbl}.results', delimiter=' ', comments='#')

    nthreads = d[:, 0]
    raw      = d[:, 1:]
    nruns    = raw.shape[1]
    
    speedup = np.mean(raw[0]) / raw

    stderr  = np.std(speedup, axis=1) / np.sqrt(nruns)
    speedup = np.mean(speedup, axis=1)

    print(f'exp_num_threads_speedup({lbl}): max(speedup)={max(speedup)}')

    plt.errorbar(nthreads, speedup, yerr=stderr, capsize=2, elinewidth=1, ecolor='red')
    plt.xlabel('number of threads')
    plt.ylabel('speedup (relative to 1 thread)')
    plt.tight_layout()


def exp_num_threads_speedup_multi(variate="cols", param_vals=[1000, 5000, 10000]):
    plt.figure()

    ncols = nrows = 1000

    for c in param_vals:
        if variate == 'cols': ncols = c
        else:                 nrows = c

        d = np.loadtxt(f'num_threads_{nrows}x{ncols}_20runs.results', delimiter=' ', comments='#')

        nthreads = d[:, 0]
        raw      = d[:, 1:]
        nruns    = raw.shape[1]
        
        speedup = np.mean(raw[0]) / raw
            
        stderr  = np.std(speedup, axis=1) / np.sqrt(nruns)
        speedup = np.mean(speedup, axis=1)

        lbl = f'${nrows} \\times {ncols}$'

        plt.errorbar(nthreads, speedup, yerr=stderr, 
            capsize=2, elinewidth=1, label=lbl)
    
    plt.legend(title="Grid size ($r\\times c$)")
    plt.xlabel('number of threads')
    plt.ylabel('speedup (relative to 1 thread)')
    plt.tight_layout()



# -- Runtime=f(nrows, ncols) --
exp_variate_grid(nthreads=4)
plt.savefig('figures/variate_grid_4t.svg')

exp_variate_grid(nthreads=96)
plt.savefig('figures/variate_grid_96t.svg')

# -- Runtime=f(n_threads) --
exp_num_threads('25Melem_5runs')
plt.savefig('figures/num_threads_25Melem_5runs.svg')

exp_num_threads('25Melem_10runs')
plt.savefig('figures/num_threads_25Melem_10runs.svg')

exp_num_threads('25Melem_20runs')
plt.savefig('figures/num_threads_25Melem_20runs.svg')

exp_num_threads('1000x1000_20runs')
plt.savefig('figures/num_threads_1Melem_20runs.svg')

# -- Speedup=f(n_threads) --
exp_num_threads_speedup('25Melem_10runs')
plt.savefig('figures/num_threads_speedup_25Melem_10runs.svg')

exp_num_threads_speedup('25Melem_5runs')
plt.savefig('figures/num_threads_speedup_25Melem_5runs.svg')

exp_num_threads_speedup('25Melem_20runs')
plt.savefig('figures/num_threads_speedup_25Melem_20runs.svg')

exp_num_threads_speedup('1000x1000_20runs')
plt.savefig('figures/num_threads_speedup_1Melem_20runs.svg')

# -- Speedup, with multiple element numbers --
exp_num_threads_speedup_multi(variate='rows')
plt.savefig('figures/num_threads_speedup_multi-variate-rows_20runs.svg')

exp_num_threads_speedup_multi(variate='cols', param_vals=[1000, 5000, 10000])
plt.savefig('figures/num_threads_speedup_multi-variate-cols_20runs.svg')