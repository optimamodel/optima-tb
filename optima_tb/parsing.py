import ast
import numpy as np

supported_functions = {'exp': np.exp, 'floor': np.floor}

def parse_function(fcn_str):
    # Returns (fcn,dep_list)
    # Where dep_list corresponds to a list of keys for
    # the dict that needs to be passed to fcn()
    # supported_functions is a dict mapping ast names to functors imported in the namespace of this file
    assert '__' not in fcn_str, 'Cannot use double underscores in functions'
    fcn_str = fcn_str.replace(':','___')
    fcn_str = fcn_str.replace('^','**') # Adapter for programs
    fcn_ast = ast.parse(fcn_str, mode='eval')
    dep_list = []
    for node in ast.walk(fcn_ast):
        if isinstance(node, ast.Name) and node.id not in supported_functions:
            dep_list.append(node.id)
        elif isinstance(node,ast.Call) and hasattr(node,'func') and hasattr(node.func,'id'):
            assert node.func.id in supported_functions, "Only calls to supported functions are allowed"
    compiled_code = compile(fcn_ast,filename="<ast>", mode="eval")
    def fcn(**deps):
        return eval(compiled_code,deps,supported_functions)
    return (fcn,dep_list)