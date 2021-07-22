import logging

def make_logger(file, level = logging.INFO):
    fmt = logging.Formatter('[%(levelname)s] %(message)s')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(fmt)

    fmt = logging.Formatter('[%(asctime)s %(levelname)s] %(message)s', datefmt = '%Y-%m-%d %H:%M:%S')
    file_handler = logging.FileHandler(file)
    file_handler.setFormatter(fmt)
    
    logger = logging.getLogger('main')
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)
    logger.setLevel(level)
    return logger