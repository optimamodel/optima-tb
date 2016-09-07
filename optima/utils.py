#%% Imports

from collections import OrderedDict
from numpy import array
from numbers import Number


#%% Homebrew odict class

class odict(OrderedDict):
    '''
    An ordered dictionary, like the OrderedDict class, but supporting list methods like integer referencing, slicing, and appending.
    Version: 2016feb09 (cliffk)
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
        if isinstance(key, Number): # Convert automatically from float...dangerous?
            return self.values()[int(key)]
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
        else: # Handle string but also everything else
            try:
                output = OrderedDict.__getitem__(self,key)
                return output
            except: # WARNING, should be KeyError, but this can't print newlines!!!
                if len(self.keys()): 
                    errormsg = 'odict key "%s" not found; available keys are:\n%s' % (str(key), 
                        '\n'.join([str(k) for k in self.keys()]))
                else: errormsg = 'Key "%s" not found since odict is empty'% key
                raise Exception(errormsg)
        
        
    def __setitem__(self, key, value):
        ''' Allows setitem to support strings, integers, slices, lists, or arrays '''
        if type(key)==str:
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
