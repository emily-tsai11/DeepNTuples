from pprint import pprint
#from btv_flightDistance_recalculation import C_recalculate_flightDistance, C_rename_flightDistance
import btv_vc1Fix

def recalc(name):
#    return C_recalculate_flightDistance(name)
    return btv_vc1Fix.rootf(name)

def rename(name):
#    return C_rename_flightDistance(name)
    return btv_vc1Fix.rootf_2(name)

import os
from multiprocessing import Pool
files = os.listdir('.')
# print(files)
root_files = []
for file_name in files:
    if 'ntuple' in file_name and file_name.endswith('root'):
        root_files.append(file_name)

from scripts_base import do_for_all_parallel
pprint(root_files)
do_for_all_parallel(recalc, [root_files])
do_for_all_parallel(rename, [root_files])

print('Done!')

