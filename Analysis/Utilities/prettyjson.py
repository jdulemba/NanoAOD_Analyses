'''
Few json functions made prettier
Author: Mauro Verzetti, edited by Joseph Dulemba to work with python3
'''
import json
from pdb import set_trace

def dumps(obj, indent = 4, separators = (',', ': '), **kwargs):
    return json.dumps(obj, indent=indent, separators=separators, **kwargs)

def convert(input):
    '''converts json unicode strings into bytecode strings'''
    #set_trace()
    if isinstance(input, dict):
        return dict([(convert(key), convert(value)) for key, value in input.items()])
    elif isinstance(input, list):
        return [convert(element) for element in input]
    #elif isinstance(input, str):
    #    return input.decode()
    #elif isinstance(input, unicode):
    #    return input.encode('utf-8')
    else:
        return input  


def loads( input, object_hook=convert, **kwargs):
    return json.loads( input, object_hook=object_hook, **kwargs )
