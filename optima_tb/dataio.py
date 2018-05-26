from optima_tb.utils import OptimaException

#try: import cPickle as pickle   # For Python 2 compatibility
#except: import pickle
import dill
from gzip import GzipFile
from os import path

def saveObject(filename, obj, compresslevel=5, verbose=True):
    ''' Save an object to file -- use compression 5, since more is much slower but not much smaller '''
    with GzipFile(filename, 'wb', compresslevel=compresslevel) as fileobj:
        fileobj.write(dill.dumps(obj, protocol=-1))
    if verbose: print('Object saved to "%s"' % filename)
    return path.abspath(filename)

def loadObject(filename, verbose=True):
    ''' Load a saved file '''
    # Handle loading of either filename or file object
    if isinstance(filename, basestring): argtype='filename'
    else: argtype = 'fileobj'
    kwargs = {'mode': 'rb', argtype: filename}
    with GzipFile(**kwargs) as fileobj:
        obj = loadPickle(fileobj)
    if verbose: print('Object loaded from "%s"' % filename)
    return obj

def loadPickle(fileobj, verbose=False):
    ''' Loads a pickled object -- need to define legacy classes here since they're needed for unpickling '''
    
    # Load the file string
    filestr = fileobj.read()
    obj = dill.loads(filestr)
    return obj
