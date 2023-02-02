import os

def restart(outpath, restart_path, minima_path, is_master):
        
    if not os.path.exists(outpath):
        os.mkdir(outpath)
        is_output = False
    else:
        is_output = True
        
    if not os.path.exists(restart_path):
        os.mkdir(restart_path)
        is_restart = False
    else:
        is_restart = True

    if not os.path.exists(minima_path):
        os.mkdir(minima_path)

    if is_restart and is_output:
        is_restart = True
        is_restart = checkfiles(restart_path, is_master)
    else:
        is_restart = False
    
    return is_restart


def checkfiles(restart_path, is_master):
    is_files = True
    if is_master:
        is_database = os.path.exists(restart_path + "minima.pickle.shelve.dat") or os.path.exists(restart_path + "minima.pickle.shelve.dat.dat")
        if not is_database:
            is_files = False
    
    if not is_master:
        is_minimum = os.path.exists(restart_path + "poscur.extxyz")
        if not is_minimum:
            is_files = False

    is_params = os.path.exists(restart_path + "params.json")
    if not is_params:
        is_files = False

    return is_files

