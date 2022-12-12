from minimahopping.mh.database import Database
from mpi4py import MPI
import minimahopping.MPI_database.mpi_messages as message


def MPI_database_server_loop(energy_threshold, minima_threshold, output_n_lowest_minima, is_restart = False, outpath='./', minima_path= "lowest_minima/", write_graph_output = True):
    status = MPI.Status()
    with Database(energy_threshold, minima_threshold, output_n_lowest_minima, is_restart, outpath, minima_path, write_graph_output) as db:
        comm_world = MPI.COMM_WORLD
        while True:
            message_tag, data = comm_world.recv(status=status)
            sender = status.Get_source()
            print('message received ', message_tag, sender)

            if message_tag == message.addelement:
                print('message', data, type(data))
                n_visit, label = db.addElement(data)
                comm_world.send((n_visit, label), sender)
            elif message_tag == message.addElementandConnectGraph:
                n_visit, label = db.addElementandConnectGraph(*data)
                comm_world.send((n_visit, label), sender)
            elif message_tag == message.get_element_index:
                index = db.get_element_index(data)
                comm_world.send(index, sender)
            elif message_tag == message.get_element:
                requested_minimum = db.get_element(data)
                comm_world.send(requested_minimum, sender)
            else:
                print('tag not known, shutting down')
                quit()
