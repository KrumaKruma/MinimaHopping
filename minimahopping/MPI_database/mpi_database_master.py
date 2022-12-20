from minimahopping.mh.database import Database
from mpi4py import MPI
import minimahopping.MPI_database.mpi_messages as message
import time
import numpy as np


def MPI_database_server_loop(energy_threshold, minima_threshold, output_n_lowest_minima, is_restart = False, outpath='./', minima_path= "lowest_minima/", write_graph_output = True, maxTimeHours = np.inf, totalWorkers=None):

    if totalWorkers is None:
        print("total number of workers must be given to mpi_database_serverloop. aborting...")
        quit()
    
    current_workers = totalWorkers

    maxTimeSeconds = maxTimeHours * 3600
    t_start = time.time()

    status = MPI.Status()
    with Database(energy_threshold, minima_threshold, output_n_lowest_minima, is_restart, outpath, minima_path, write_graph_output) as db:
        comm_world = MPI.COMM_WORLD
        process_time = 0
        wait_time = 0
        t1 = time.time()

        # this set contains the ranks of all clients that stopped working.
        # When the size of this set is totalWorkers the server can be stopped
        stoppedClients = set()
        while True:
            continueSimulation = time.time() - t_start < maxTimeSeconds

            print('server_efficiency', process_time / (process_time + wait_time+1e-5), file=open('efficiency.txt', mode='a'))
            t2 = time.time()
            process_time += t2 - t1
            t1 = time.time()
            message_tag, data = comm_world.recv(status=status)
            t2 = time.time()
            sender = status.Get_source()
            wait_time += t2 - t1

            t1 = time.time()

            if message_tag == message.addelement:
                n_visit, label = db.addElement(data)
                comm_world.send((n_visit, label), sender)
            elif message_tag == message.addElementandConnectGraph:
                n_visit, label, temp = db.addElementandConnectGraph(*data)
                comm_world.send((n_visit, label, continueSimulation), sender)
                if not continueSimulation:
                    if sender in stoppedClients:
                        print('Server thinks client stopped working, but client is still running.')
                        print("rank of bad client:", sender)
                        MPI.COMM_WORLD.Abort()
                    stoppedClients.add(sender)
                    current_workers = current_workers - 1
            elif message_tag == message.get_element_index:
                index = db.get_element_index(data)
                comm_world.send(index, sender)
            elif message_tag == message.get_element:
                requested_minimum = db.get_element(data)
                comm_world.send(requested_minimum, sender)
            elif message_tag == message.clientWorkDone:
                if not sender in stoppedClients: # client stopped running and was not stopped by server.
                    stoppedClients.add(sender)
                    current_workers = current_workers - 1
            else:
                print('tag not known, shutting down')
                quit()

            if current_workers <= 0: # Last process that is still alive. The simulation can be stopped.
                return