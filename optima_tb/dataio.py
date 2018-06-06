import pickle
import dill
from gzip import GzipFile
from os import path
import logging
logger = logging.getLogger(__name__)

def saveObject(filename, obj, compresslevel=5,method='pickle'):
    ''' Save an object to file -- use compression 5, since more is much slower but not much smaller '''
    with GzipFile(filename, 'wb', compresslevel=compresslevel) as fileobj:
        if method == 'dill':
            fileobj.write(dill.dumps(obj, protocol=-1))
        else:
            try:
                fileobj.write(pickle.dumps(obj, protocol=-1))
            except Exception as e:
                logger.warn(str(e))
                fileobj.write(dill.dumps(obj, protocol=-1))
    logger.info('Object saved to "%s"' % filename)
    return path.abspath(filename)

def loadObject(filename):
    ''' Load a saved file '''
    # Handle loading of either filename or file object
    if isinstance(filename, basestring):
        argtype='filename'
    else:
        argtype = 'fileobj'

    kwargs = {'mode': 'rb', argtype: filename}
    with GzipFile(**kwargs) as fileobj:
        filestr = fileobj.read() # Convert it to a string
        try:
            obj = pickle.loads(filestr) # Actually load it
        except Exception as e:
            logger.warn(str(e))
            obj = dill.loads(filestr)

    logger.info('Object loaded from "%s"' % filename)
    return obj
