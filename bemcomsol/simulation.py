import mph
import os
import numpy as np
from shutil import rmtree
from itertools import product
import pandas as pd
import re

class Parameter:
    def __init__(self, name, expression, unit=''):
        self.name = name
        self.expr = expression
        self.unit = unit

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

def load(path, datalabel, coordinates=3):
    '''
        Input
        - *path*: Directory path where you want to import CSV file from.
        - *datalabel*: Data labels of CSV file. The elements of dataLabel should match the number of columns, and seperated by comma ','. It should begin with coordinate variables, and then physics quantity. Example: 'x,y,z,esbeV' or 'x,y,Ex,Ey'.
        - *coordinates*: Number of coordinate variables. For example, 'x,y,Ex,Ey' should be set as 2. Default value is 3.

        Output: A list of object
        - For each coordinate variable, it creates a numpy array recording the coordinate data.
        - For each physics quantity, it creates a dict named with 'label' + 's', whose keys are 'DCi'.

        Example:
        If the data outputted by COMSOL contains five columns, the first three columns represent spatial coordinates,
        and the last two columns respectively record the potential and electric field, after loading it with 
        `load(path, 'x,y,z,V,E', coordinates=3)` , we will get 3 numpy arrays `x` , `y` and `z` recording the coordinates,
        and 2 dicts, `Vs` and `Es` . We can use `Vs['DC3']` to access the electric potential when DC3 is set to 1[V] and others set to 0[V].
    '''
    cwd = os.path.join(os.getcwd(), path)
    files = os.listdir(cwd)
    DCnum = len(files)
    labels = datalabel.split(',')
    for i in range(coordinates,len(labels)):
        vars = labels[i] + 's'
        exec(f"{vars} = dict()")
    for file in files:
        with open(os.path.join(cwd,file), 'r') as COMSOL: 
            lines = COMSOL.readlines()
        newlines = [row for row in lines if row[0]!='%'] # Delete annotation
        with open('output.csv', 'w') as csvdata:
            csvdata.write(datalabel+'\n')
            csvdata.writelines(newlines)
        data = pd.read_csv('output.csv')
        DCi = file.strip('.csv')
        for i in range(coordinates,len(labels)):
            var = labels[i]
            vars = var + 's'
            exec(f"{vars}['{DCi}'] = data['{var}'].values")
    for i in range(coordinates):
        coordinate = labels[i]
        exec(f"{coordinate} = data['{coordinate}'].values")
    datalist = []
    for i in range(len(labels)):
        if i < coordinates:
            vars = labels[i]
        else:
            vars = labels[i] + 's'
        exec(f'datalist.append({vars})')
    os.remove('output.csv')
    return datalist