import logging
import minimahopping.mh.parameters

logger = None

def setupLogger(logLevel: int, file: str = None):
    global logger
    logger = logging.getLogger('Minima hopping logger')
    logger.setLevel(logLevel)
    formatter = logging.Formatter('%(message)s')
    
    if file == None:
        # create console handler and set level to debug
        ch = logging.StreamHandler()
    else:
        ch = logging.FileHandler(file)
    ch.setLevel(logLevel)
    ch.setFormatter(formatter)
    logger.addHandler(ch)