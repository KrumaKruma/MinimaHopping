import logging
import minimahopping.mh.parameters

logger = logging.getLogger('Minima hopping logger')
logger.propagate = False  # Avoid duplicate logging from root
logger.setLevel(logging.DEBUG)  # Default level

def setupLogger(logLevel: int, file: str = None):
    global logger
    logger.setLevel(logLevel)
    formatter = logging.Formatter('%(message)s')

    if logger.hasHandlers():
        logger.handlers.clear()  # Avoid duplicate handlers
    if file == None:
        ch = logging.StreamHandler()
    else:
        ch = logging.FileHandler(file)
    ch.setLevel(logLevel)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

# add some tests
if __name__ == '__main__':
    logger.critical('Criticall message before setting up')
    logger.debug('debeug before setting up')
    setupLogger(logging.DEBUG)
    logger.debug('debug message to console')
    setupLogger(logging.INFO, 'toto')
    logger.info('write to file')