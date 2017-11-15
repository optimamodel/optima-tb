#%% Imports
import logging, functools, traceback
from collections import OrderedDict
from numpy import array
from numbers import Number
import numpy as np
import multiprocessing as mp

logger = logging.getLogger(__name__)


#%% Timing functions

def tic():
    '''
    First in a pair of functions for calculating time difference, similar to MATLAB...
    t = tic()
    toc(t)
    '''
    from time import time
    return time()

def toc(start=0, label='', sigfigs=3):
    '''
    Second in a pair of functions for calculating time difference, similar to MATLAB...
    t = tic()
    toc(t)
    '''
    from time import time
    elapsed = time() - start
    if label=='': base = 'Elapsed time: '
    else: base = 'Elapsed time for %s: ' % label
    logger.debug(base + '%.*f s' % (sigfigs, elapsed))
    return elapsed



#%% A printing function for various levels of verbosity

def printv(string, thisverbose=1, verbose=2, newline=True, indent=False):
    '''
    Optionally print a message and automatically indent. The idea is that
    a global or shared "verbose" variable is defined, which is passed to
    subfunctions, determining how much detail to print out.

    The general idea is that verbose is an integer from 0-4 as follows:
        0 = no printout whatsoever
        1 = only essential warnings, e.g. suppressed exceptions
        2 = standard printout
        3 = extra debugging detail (e.g. printout on each iteration)
        4 = everything possible (e.g. printout on each timestep)
    
    Thus a very important statement might be...
        printv('WARNING, everything is wrong', 1, verbose)
    ...whereas a much less important message might be...
        printv('This is timestep %i' % i, 4, verbose)
    '''
    
    if thisverbose>4 or verbose>4: print('WARNING: Verbosity should be from 0-4. This message: %i. Current: %i.' % (thisverbose, verbose))
    if verbose>=thisverbose: # Only print if sufficiently verbose
        indents = '  '*thisverbose*bool(indent) # Create automatic indenting
        if newline: print(indents+str(string)) # Actually print
        else: print(indents+str(string)), # Actually print
    return None



#%% Functions for dealing with recursive structures

def flattenDict(input_dict, base_key, sub_keys = None, comp_list = None, key_list = None, limit = 100):
    '''
    A function for flattening out a recursive dictionary, with an optional list of sub-keys (ignored if non-existent).
    The flattened out structure is returned as comp_list. Values can be an object or a list of objects.
    All keys (including base_key) within the recursion are returned as key_list.
    
    Specifically, this function is intended for dictionaries of the form...
        input_dict[key1][sub_key[0]] = [a, key2, b]
        input_dict[key1][sub_key[1]] = [c, d]
        input_dict[key2][sub_key[0]] = e
        input_dict[key2][sub_key[1]] = [e, f, g]
    ...which, for this specific example, will output list...
        [a, e, e, f, g, h, b, c, d]
        
    There is a max-depth of limit for the recursion.
    '''
    
    if limit < 1:
        raise OptimaException('ERROR: A recursion limit has been reached when flattening a dictionary, stopping at key "%s".' % base_key)    
    
    if comp_list is None: comp_list = []
    if key_list is None: key_list = []
    key_list.append(base_key)

    if sub_keys is None: input_list = input_dict[base_key]
    else:
        input_list = []
        for sub_key in sub_keys:
            if sub_key in input_dict[base_key]:
                val = input_dict[base_key][sub_key]
                if isinstance(val, list):
                    input_list += val
                else:
                    input_list.append(val)      # Handle unlisted objects.
    
    for comp in input_list:
        if comp in input_dict.keys():
            flattenDict(input_dict = input_dict, base_key = comp, sub_keys = sub_keys, comp_list = comp_list, key_list = key_list, limit = limit - 1)
        else:
            comp_list.append(comp)
    return comp_list, key_list
    
#%% docstring methods
def createcollist(oldkeys, title, strlen = 30, ncol = 3):
    ''' Creates a string for a nice columnated list (e.g. to use in __repr__ method) '''
    from numpy import ceil
    nrow = int(ceil(float(len(oldkeys))/ncol))
    newkeys = []
    for x in xrange(nrow):
        newkeys += oldkeys[x::nrow]
    
    attstring = title + ':'
    c = 0    
    for x in newkeys:
        if c%ncol == 0: attstring += '\n  '
        if len(x) > strlen: x = x[:strlen-3] + '...'
        attstring += '%-*s  ' % (strlen,x)
        c += 1
    attstring += '\n'
    return attstring
    
def objectid(obj):
    ''' Return the object ID as per the default Python __repr__ method '''
    return '%s\n<%s.%s at %s>\n' % (obj.name, obj.__class__.__module__, obj.__class__.__name__, hex(id(obj)))


def objatt(obj, strlen = 30, ncol = 3):
    ''' Return a sorted string of object attributes for the Python __repr__ method '''
    oldkeys = sorted(obj.__dict__.keys())
    return createcollist(oldkeys, 'Attributes', strlen = 30, ncol = 3)


def objmeth(obj, strlen = 30, ncol = 3):
    ''' Return a sorted string of object methods for the Python __repr__ method '''
    oldkeys = sorted([method + '()' for method in dir(obj) if callable(getattr(obj, method)) and not method.startswith('__')])
    return createcollist(oldkeys, 'Methods', strlen = 30, ncol = 3)
        
def objrepr(obj, showid=True, showmeth=True, showatt=True):
    ''' Return useful printout for the Python __repr__ method '''
    divider = '====================================================================================\n'
    output = ''
    if showid:
        output += objectid(obj)
        output += divider
    if showmeth:
        output += objmeth(obj)
        output += divider
    if showatt:
        output += objatt(obj)
        output += divider
    return output
    
def defaultrepr(obj, maxlen=55):
    ''' Prints out the default representation of an object -- all attributes, plus methods and ID '''
    keys = sorted(obj.__dict__.keys())
    maxkeylen = max([len(key) for key in keys])
    if maxkeylen<maxlen: maxlen = maxlen - maxkeylen
    formatstr = '%'+ '%i'%maxkeylen + 's'
    output  = objrepr(obj, showatt=False)
    for key in keys:
        thisattr = str(getattr(obj, key))
        if len(thisattr)>maxlen: thisattr = thisattr[:maxlen] + ' [...]'
        output += formatstr%key + ': ' + thisattr + '\n'
    output += '====================================================================================\n'

    return output
    
#%% Homebrew odict class

class odict(OrderedDict):
    '''
    An ordered dictionary, like the OrderedDict class, but supporting list methods like integer referencing, slicing, and appending.
    Version: 2016sep14 (cliffk)
    '''

    def __slicekey(self, key, slice_end):
        shift = int(slice_end=='stop')
        if isinstance(key, Number): return key
        elif type(key) is str: return self.index(key)+shift # +1 since otherwise confusing with names (CK)
        elif key is None: return (len(self) if shift else 0)
        else: raise Exception('To use a slice, %s must be either int or str (%s)' % (slice_end, key))


    def __is_odict_iterable(self, v):
        return type(v)==list or type(v)==type(array([]))


    def __getitem__(self, key):
        ''' Allows getitem to support strings, integers, slices, lists, or arrays '''
        if isinstance(key, (str,tuple)):
            try:
                output = OrderedDict.__getitem__(self,key)
                return output
            except: # WARNING, should be KeyError, but this can't print newlines!!!
                if len(self.keys()): 
                    errormsg = 'odict key "%s" not found; available keys are:\n%s' % (str(key), 
                        '\n'.join([str(k) for k in self.keys()]))
                else: errormsg = 'Key "%s" not found since odict is empty'% key
                raise Exception(errormsg)
        elif isinstance(key, Number): # Convert automatically from float...dangerous?
            thiskey = self.keys()[int(key)]
            return OrderedDict.__getitem__(self,thiskey)
        elif type(key)==slice: # Handle a slice -- complicated
            try:
                startind = self.__slicekey(key.start, 'start')
                stopind = self.__slicekey(key.stop, 'stop')
                if stopind<startind:
                    print('Stop index must be >= start index (start=%i, stop=%i)' % (startind, stopind))
                    raise Exception
                slicevals = [self.__getitem__(i) for i in range(startind,stopind)]
                try: return array(slicevals) # Try to convert to an array
                except: return slicevals
            except:
                print('Invalid odict slice... returning empty list...')
                return []
        elif self.__is_odict_iterable(key): # Iterate over items
            listvals = [self.__getitem__(item) for item in key]
            try: return array(listvals)
            except: return listvals
        else: # Handle everything else
            return OrderedDict.__getitem__(self,key)
        
        
    def __setitem__(self, key, value):
        ''' Allows setitem to support strings, integers, slices, lists, or arrays '''
        if isinstance(key, (str,tuple)):
            OrderedDict.__setitem__(self, key, value)
        elif isinstance(key, Number): # Convert automatically from float...dangerous?
            thiskey = self.keys()[int(key)]
            OrderedDict.__setitem__(self, thiskey, value)
        elif type(key)==slice:
            startind = self.__slicekey(key.start, 'start')
            stopind = self.__slicekey(key.stop, 'stop')
            if stopind<startind:
                errormsg = 'Stop index must be >= start index (start=%i, stop=%i)' % (startind, stopind)
                raise Exception(errormsg)
            slicerange = range(startind,stopind)
            enumerator = enumerate(slicerange)
            slicelen = len(slicerange)
            if hasattr(value, '__len__'):
                if len(value)==slicelen:
                    for valind,index in enumerator:
                        self.__setitem__(index, value[valind])  # e.g. odict[:] = arr[:]
                else:
                    errormsg = 'Slice "%s" and values "%s" have different lengths! (%i, %i)' % (slicerange, value, slicelen, len(value))
                    raise Exception(errormsg)
            else: 
                for valind,index in enumerator:
                    self.__setitem__(index, value) # e.g. odict[:] = 4
        elif self.__is_odict_iterable(key) and hasattr(value, '__len__'): # Iterate over items
            if len(key)==len(value):
                for valind,thiskey in enumerate(key): 
                    self.__setitem__(thiskey, value[valind])
            else:
                errormsg = 'Keys "%s" and values "%s" have different lengths! (%i, %i)' % (key, value, len(key), len(value))
                raise Exception(errormsg)
        else:
            OrderedDict.__setitem__(self, key, value)
        return None
    
    
    def pop(self, key, *args, **kwargs):
        ''' Allows pop to support strings, integers, slices, lists, or arrays '''
        if type(key)==str:
            return OrderedDict.pop(self, key, *args, **kwargs)
        elif isinstance(key, Number): # Convert automatically from float...dangerous?
            thiskey = self.keys()[int(key)]
            return OrderedDict.pop(self, thiskey, *args, **kwargs)
        elif type(key)==slice: # Handle a slice -- complicated
            try:
                startind = self.__slicekey(key.start, 'start')
                stopind = self.__slicekey(key.stop, 'stop')
                if stopind<startind:
                    print('Stop index must be >= start index (start=%i, stop=%i)' % (startind, stopind))
                    raise Exception
                slicevals = [self.pop(i, *args, **kwargs) for i in range(startind,stopind)] # WARNING, not tested
                try: return array(slicevals) # Try to convert to an array
                except: return slicevals
            except:
                print('Invalid odict slice... returning empty list...')
                return []
        elif self.__is_odict_iterable(key): # Iterate over items
            listvals = [self.pop(item, *args, **kwargs) for item in key]
            try: return array(listvals)
            except: return listvals
        else: # Handle string but also everything else
            try:
                return OrderedDict.pop(self, key, *args, **kwargs)
            except: # WARNING, should be KeyError, but this can't print newlines!!!
                if len(self.keys()): 
                    errormsg = 'odict key "%s" not found; available keys are:\n%s' % (str(key), 
                        '\n'.join([str(k) for k in self.keys()]))
                else: errormsg = 'Key "%s" not found since odict is empty'% key
                raise Exception(errormsg)
    
    """"
    def __repr__(self, maxlen=None, spaces=True, divider=True):
        ''' Print a meaningful representation of the odict '''
         # Maximum length of string to display
        toolong = ' [...]'
        divider = '#############################################################\n'
        if len(self.keys())==0: 
            output = 'odict()'
        else: 
            output = ''
            hasspaces = 0
            for i in range(len(self)):
                if divider and spaces and hasspaces: output += divider
                thiskey = str(self.keys()[i]) # Probably don't need to cast to str, but just to be sure
                thisval = str(self.values()[i])
                if not(spaces):                    thisval = thisval.replace('\n','\\n') # Replace line breaks with characters
                if maxlen and len(thisval)>maxlen: thisval = thisval[:maxlen-len(toolong)] + toolong # Trim long entries
                if thisval.find('\n'): hasspaces = True
                output += '#%i: "%s": %s\n' % (i, thiskey, thisval)
        return output
    """
    def __repr__(self, maxlen=None, spaces=True, divider=True,indent=2,indentLevel=0):
        return self.recursive_repr(indent, indentLevel)
        
    
    def recursive_repr(self,indent=2,indentLevel=0):
        offset = ' '*indent*indentLevel
        if len(self.keys())==0: 
            output = offset + 'empty <odict()>'
        else: 
            output = ''
            for k,v in self.iteritems():
                output += offset + str(k) 
                if isinstance(v,odict):
                    output += '\n' + v.__repr__(indent=indent,indentLevel=indentLevel+1)
                else:
                    output += ':\t' + str(v)
                output += '\n'
        return output
    
    
    def disp(self, maxlen=55, spaces=False, divider=False):
        ''' Print out flexible representation, short by default'''
        print(self.__repr__(maxlen=maxlen, spaces=spaces, divider=divider))
    
    
    def _repr_pretty_(self, p, cycle):
        ''' Stupid function to fix __repr__ because IPython is stupid '''
        print(self.__repr__())
    
    
    def index(self, item):
        ''' Return the index of a given key '''
        return self.keys().index(item)
    
    
    def valind(self, item):
        ''' Return the index of a given value '''
        return self.items().index(item)
    
    
    def append(self, item):
        ''' Support an append method, like a list '''
        keyname = str(len(self)) # Define the key just to be the current index
        self.__setitem__(keyname, item)
        return None
    
    
    def rename(self, oldkey, newkey):
        ''' Change a key name -- WARNING, very inefficient! '''
        nkeys = len(self)
        if isinstance(oldkey, Number): 
            index = oldkey
            keystr = self.keys()[index]
        elif type(oldkey) is str: 
            index = self.keys().index(oldkey)
            keystr = oldkey
        else: raise Exception('Key type not recognized: must be int or str')
        self.__setitem__(newkey, self.pop(keystr))
        if index<nkeys-1:
            for i in range(index+1, nkeys):
                key = self.keys()[index]
                value = self.pop(key)
                self.__setitem__(key, value)
        return None
    
    
    def sort(self, sortby=None):
        ''' Return a sorted copy of the odict. 
        Sorts by order of sortby, if provided, otherwise alphabetical'''
        if not sortby: allkeys = sorted(self.keys())
        else:
            if not isinstance(sortby, list): raise Exception('Please provide a list to determine the sort order.')
            if all(isinstance(x,basestring) for x in sortby): # Going to sort by keys
                if not set(sortby)==set(self.keys()): 
                    errormsg = 'List of keys to sort by must be the same as list of keys in odict.\n You provided the following list of keys to sort by:\n'
                    errormsg += '\n'.join(sortby)
                    errormsg += '\n List of keys in odict is:\n'
                    errormsg += '\n'.join(self.keys())
                    raise Exception(errormsg)
                else: allkeys = sortby
            elif all(isinstance(x,int) for x in sortby): # Going to sort by numbers
                if not set(sortby)==set(range(len(self))):
                    errormsg = 'List to sort by "%s" is not compatible with length of odict "%i"' % (sortby, len(self))
                    raise Exception(errormsg)
                else: allkeys = [y for (x,y) in sorted(zip(sortby,self.keys()))]
            else: raise Exception('Cannot figure out how to sort by "%s"' % sortby)
        out = odict()
        for key in allkeys: out[key] = self[key]
        return out



#%% An exception wrapper class

class OptimaException(Exception):
    ''' A wrapper class to allow for Optima-specific exceptions. Can be expanded. '''
    def __init(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


# Wrapper for parallel tasks
def runParallel(func, params, num_threads=2):
    if num_threads == None:
        pool = mp.Pool()
    else:
        pool = mp.Pool(num_threads)
    pool.map_async(func, params)
    pool.close()
    pool.join()

# Allows stacktrace in threads. HOW TO USE: Decorate any function you wish to call with runParallel (the function func
# references) with '@trace_exception'. Whenever an exception is thrown by the decorated function when executed parallel,
# a stacktrace is printed; the thread terminates but the execution of other threads is not affected.
def trace_exception(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            traceback.print_exc()
    return wrapper



# Convenience data structure to store a simulation interval from 'start' to 'stop'. 'step' is interpreted as time steps
# per year, e.g. step=365 represents a daily time step whereas step=12 represents a monthly time step
class SimInt:
    def __init__(self, start, stop, step=365):
        self.start = int(start)
        self.stop = int(stop)
        self.step = int(step)

    # returns the number of time steps between start and stop
    def length(self):
        return self.stop - self.start

    # returns the number of years contained (rounded up) by this data structure or from any given starting point in time
    def numCycles(self, start=None):
        # NOTE: typecasting for this matter is superfluous in python 3.x
        if start == None:
            return int(np.ceil(float(self.length()) / float(self.step)))
        else:
            return int(np.ceil((float(self.stop - start)) / float(self.step)))
