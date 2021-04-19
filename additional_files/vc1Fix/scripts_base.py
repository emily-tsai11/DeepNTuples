import os

import itertools
from multiprocessing import Process, Pool, cpu_count

def do_for_all_parallel(func=print,
                        loops=[['el', 'mu'], ['2j1t', '3j1t', '3j2t']],
                        *args,
                        **kwargs):
    '''
    Same as do_for_all but all jobs are run in parallel.
    Jobs are niced to prevent possibly blocking the machines, if too many jobs are started at once.
    By default Pool should start as many jobs as there are CPUs on the machine.
    '''
    product = itertools.product(*loops)

    pool = Pool()
    for loop_args in product:
        arguments = list(loop_args) + list(args)
        pool.apply_async(func=func, args=arguments, kwds=kwargs)

    pool.close()
    pool.join()
    pool.terminate()

