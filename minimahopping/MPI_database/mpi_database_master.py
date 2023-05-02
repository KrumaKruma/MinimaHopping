from minimahopping.mh.database import Database
import mpi4py
mpi4py.rc.recv_mprobe = False
from mpi4py import MPI
import minimahopping.MPI_database.mpi_messages as message
import time
import numpy as np
import logging


def MPI_database_server_loop(energy_threshold, minima_threshold, output_n_lowest_minima, is_restart = False, outpath='./', minima_path= "lowest_minima/", write_graph_output = True, maxTimeHours = np.inf, maxNumberOfMinima = 0):

    current_workers = 0
    clientHasLeft = False

    maxTimeSeconds = maxTimeHours * 3600
    t_start = time.time()

    status = MPI.Status()
    with Database(energy_threshold, minima_threshold, output_n_lowest_minima, is_restart, outpath, minima_path, write_graph_output, maxNumberOfMinima=maxNumberOfMinima) as db:
        comm_world = MPI.COMM_WORLD
        process_time = 0
        wait_time = 0
        t1 = time.time()

        # this set contains the ranks of all clients that stopped working.
        stoppedClients = set()
        while True:
            continueSimulation = time.time() - t_start < maxTimeSeconds and not clientHasLeft

            t2 = time.time()
            process_time += t2 - t1
            logging.info('Average server load: %f percent. Processing last structure took %f seconds.'%(100 * process_time / (process_time + wait_time+1e-6), t2 - t1))
            t1 = time.time()
            logging.debug("Listening for message from clients")
            message_tag, data = comm_world.recv(status=status)
            t2 = time.time()
            sender = status.Get_source()
            logging.info("Received message with tag: %s from sender %i"%(message_tag, sender))
            wait_time += t2 - t1

            t1 = time.time()

            if message_tag == message.addelement:
                n_visit, label, temp = db.addElement(data)
                comm_world.send((n_visit, label, continueSimulation), sender)
                if not continueSimulation:
                    logging.info("Sent shutdown message to client %i"%sender)
                    if sender in stoppedClients:
                        logging.error('Server thinks client stopped working, but client is still running.')
                        logging.error("rank of bad client: %i"%sender)
                        return
                    stoppedClients.add(sender)
                    current_workers = current_workers - 1
            elif message_tag == message.addElementandConnectGraph:
                n_visit, label, temp = db.addElementandConnectGraph(*data)
                comm_world.send((n_visit, label, continueSimulation), sender)
                if not continueSimulation:
                    logging.info("Sent shutdown message to client %i"%sender)
                    if sender in stoppedClients:
                        logging.error('Server thinks client stopped working, but client is still running.')
                        logging.error("rank of bad client: %i"%sender)
                        return
                    stoppedClients.add(sender)
                    current_workers = current_workers - 1
            elif message_tag == message.get_element_index:
                index = db.get_element_index(data)
                comm_world.send(index, sender)
            elif message_tag == message.get_element:
                requested_minimum = db.get_element(data)
                comm_world.send(requested_minimum, sender)
            elif message_tag == message.clientWorkDone:
                logging.info("Client %i has stopped working. Send a stop simulation signal to all clients on next request."%sender)
                clientHasLeft = True
                if not sender in stoppedClients: # client stopped running and was not stopped by server.
                    stoppedClients.add(sender)
                    current_workers = current_workers - 1
            elif message_tag == message.loginRequestFromClient:
                logging.info("Client %i has succesfully logged in"%sender)
                current_workers += 1
            else:
                logging.error('tag not known, shutting down')
                return

            if current_workers == 0: # Last process that is still alive. The simulation can be stopped.
                return