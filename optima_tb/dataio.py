from __future__ import print_function

from optima_tb.utils import odict as optimaOdict
from optima_tb.utils import OptimaException

import ast
import numpy
import datetime
import types
import uuid
import zlib
import math

from dateutil import parser, tz
from twisted.python.reflect import qual, namedAny

import logging
logger = logging.getLogger(__name__)

"""
Methods for importing and exporting objects. 

Wrapper functions: sjjarvis, 30jan2017

"""

COMPRESSION_EXTENSION = {'zlib' : '.Z',
                         'gzip' : '.gzip',
                         None   : '',
                         ''     : ''}

def exportObj(obj_to_export,filename=None,format='json',compression='zlib'):
    """
    
    Params:
        obj_to_export  object to be exported
        format         object format to save to (default = json)
        filename       name of file to save to (not including any compression extension)
        compression    compression algorithm. Currently, only zlib is supported.
        
    
    """
    
    # Determine object representation
    if format=='json':
        contents = dumps(obj_to_export)
    else:
        raise OptimaException("Format for export unknown (%s)"%format)
    
    # Determine compression algorithm
    if compression == 'zlib':
        compressed = zlib.compress(contents, 6)
    elif compression is None or compression == '':
        logger.debug("No compression specified / listed")
    else:
        raise OptimaException("Compression for export unknown (%s)"%compression)
       
    # Set up filename 
    filename = "%s%s"%(filename,COMPRESSION_EXTENSION[compression])

    # Finally, save the formatted object (using compression) to the file
    try:
        f = open(filename, 'wb')
        f.write(compressed)
        f.close() 
        return filename
    except Exception, e:
        logger.error("Error when trying to save file: %s"%(filename))
        logger.error(str(e))
    

def importObj(filename,format='json',compression='zlib'):
    """
    
    Params:
        filename       string (including extensions)
        format         object format to save to (default = json)
        compression    compression algorithm. Currently, only zlib is supported.
        
    """
    # Read the file
    try:
        blob = open(filename, 'rb').read()
    except Exception, e:
        logger.error("Error when trying to load file: %s"%(filename))
        raise OptimaException(str(e))
    
    
    # Determine compression algorithm
    if compression == 'zlib':
        uncompressed = zlib.decompress(blob)
    elif compression is None or compression == '':
        logger.debug("No compression specified / listed")
        uncompressed = blob
    else:
        raise OptimaException("Compression for export unknown (%s)"%compression)
 
    
    # Determine object representation
    if format=='json':
        return loads(uncompressed)
    else:
        raise OptimaException("Format for export unknown (%s)"%format)
    


def dumps(obj):
    """
    Dumps data object to json format.
    
    Taken from code originally written for Optima HIV (_serialize.py)
    
    @author: ?
    
    """

    obj_registry = {}
    id_number = [0]
    saved_types = set()

    def default(r):
        
        if isinstance(r, optimaOdict):
            o = {
                "obj": "odict",
                "val": [(default(x), default(y)) for x, y in r.iteritems()]}

        elif isinstance(r, uuid.UUID):
            o = {"obj": "UUID", "val": str(r.hex)}

        elif isinstance(r, numpy.ndarray):
            o = {"obj": "numpy.ndarray", "val": [default(x) for x in r]}

        elif isinstance(r, numpy.bool_):
            o = bool(r)

        elif r == numpy.nan:
            o = {"obj": "numpy.NaN"}

        elif isinstance(r, datetime.datetime):
            r.replace(tzinfo=tz.tzlocal())
            o = {"obj": "datetime", "val": r.isoformat(" ")}

        elif isinstance(r, tuple):
            o = {"obj": "tuple", "val": [default(x) for x in r]}

        elif isinstance(r, (str, unicode, int, long, types.NoneType, bool)):
            o = r

        elif isinstance(r, float):
            if math.isnan(r):
                o = {"obj": "float", "val": "nan"}
            else:
                o = r

        elif isinstance(r, list):
            o = [default(x) for x in r]

        elif isinstance(r, dict):
            o = {}
            for x,y in r.items():
                o[default(x)] = default(y)

        else:
            if not r in obj_registry:

                my_id = id_number[0]
                id_number[0] += 1

                obj_registry[r] = [my_id, None]

                try:
                    results = default(r.__getstate__())
                except AttributeError:
                    results = default(r.__dict__)

                q = {
                    "obj": "obj",
                    "type": qual(r.__class__),
                    "val": results,
                    "id": my_id
                }
                obj_registry[r][1] = q
                saved_types.add(qual(r.__class__))

            o = {
                "obj": "ref",
                "ref": obj_registry[r][0]
            }

        return o

    schema = default(obj)
    new_registry = {x[0]:x[1] for x in obj_registry.values()}

    dumped = repr({"registry": new_registry, "schema": schema})

    return dumped



def loads(input):
    """
    Loads json back in, to return data object
    
    Taken from code originally written for Optima HIV (_serialize.py)
    
    @author: ?
    
    """
    a = ast.literal_eval(input)

    registry = a["registry"]
    schema = a["schema"]
    loaded = {}

    def decode(o):
        if isinstance(o, list):
            return [decode(x) for x in o]

        elif isinstance(o, dict):

            if "obj" in o:

                obj = o["obj"]
                val = o.get("val")

                if obj == "ref":

                    ref = o["ref"]
                    return decode(registry[ref])

                elif obj == "obj":

                    if o["id"] in loaded:
                        return loaded[o["id"]]

                    newClass = namedAny(o["type"])

                    new = object.__new__(newClass)

                    loaded[o["id"]] = new

                    state = {x:decode(y) for x,y in val.items()}
                    try:
                        new.__setstate__(state)
                    except AttributeError:
                        new.__dict__ = state

                    return new

                elif obj == "datetime":
                    o = parser.parse(val)
                    return o

                elif obj == "UUID":
                    o = uuid.UUID(val)
                    return o

                elif obj == "numpy.ndarray":
                    o = numpy.array([decode(x) for x in val])
                    return o

                elif obj == "odict":

                    o = optimaOdict()

                    for i in val:
                        key = decode(i[0])
                        val = decode(i[1])

                        o[key] = val

                elif obj == "tuple":
                    o = tuple([decode(x) for x in val])

                elif obj == "float":
                    o = float(val)

                else:
                    assert False, str(o)[0:1000]
            else:
                return {decode(x):decode(y) for x,y in o.iteritems()}

            return o

        return o

    done = decode(schema)
    return done




