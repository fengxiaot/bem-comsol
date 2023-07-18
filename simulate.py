import mph
import os
import numpy as np
from shutil import rmtree
from classes import *
from itertools import product
import re

def remove_illegal_chars(filename):
    illegal_chars = r'[\\/:"*?<>|]'
    sanitized_filename = re.sub(illegal_chars, '', filename)
    return sanitized_filename

def simulate(mphpath, datapath, *parameters):
    '''
    Input:
    - mphpath: path of COMSOL mph file.
    - datapath: path where data is exported. Warning: ALL files in that path will be cleared.
    - parameters(optional): Param objects. If they are numbers or strings, the program will directly change it when do simulations.
    If they are lists or strings, the program will iterate based on it.
    '''

    client = mph.start()
    model = client.load(mphpath)
    paramdict = model.parameters()

    list_parameters = []
    for param in parameters:
        if isinstance(param.expr,(list,tuple,np.ndarray)):
            # If param.expr is an array, save it
            list_parameters.append(param)
        else: 
            # If param.expr is a number or string, directly set it
            model.parameter(param.name, str(param.expr) + param.unit)
    
    # Extract the expression lists from array-like parameters
    expr_list = []
    for list_parameter in list_parameters:
        expr_list.append(list_parameter.expr)
    
    # Set the parameters
    number_of_list_parameter = len(list_parameters)
    for parameter_set in product(*expr_list):
        # Set the list parameters through Cartesian product
        dirname = ''
        for i in range(number_of_list_parameter):
            if isinstance(parameter_set[i],(int,np.integer)):
                dirname = dirname + list_parameters[i].name + '=' + str(parameter_set[i]) + ','
                model.parameter(list_parameters[i].name, str(parameter_set[i]) + list_parameters[i].unit)
            elif isinstance(parameter_set[i],(float,np.floating)):
                dirname = dirname + list_parameters[i].name + '=' '{:.3f}'.format(parameter_set[i]) + ','
                model.parameter(list_parameters[i].name, str(parameter_set[i]) + list_parameters[i].unit)
            elif isinstance(parameter_set[i],(str,np.character)):
                dirname = dirname + list_parameters[i].name + '=' + parameter_set[i] + ','
                model.parameter(list_parameters[i].name, parameter_set[i] + list_parameters[i].unit)
            else:
                raise TypeError('No type matching.')
        
        # Create the directory based on parameter sets
        dirname = remove_illegal_chars(dirname.rstrip(','))
        dirpath = os.path.join(datapath,dirname)
        if os.path.exists(dirpath):
            rmtree(dirpath)
        os.makedirs(dirpath)

        # Set DC voltage sequentially and simulate
        for param in paramdict:
            if param.startswith('DC'):
                for otherparam in paramdict:
                    if otherparam.startswith('DC'):
                        model.parameter(otherparam,'0[V]')
                model.parameter(param,'1[V]')
                print(model.parameters())
                model.solve()
                for node in model/'exports':
                    node.property('filename',os.path.join(os.getcwd(), dirpath, param + '.csv'))
                    node.run()

    model.clear()
    model.reset()
    client.clear()