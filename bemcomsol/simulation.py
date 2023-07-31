import mph
import os
import numpy as np
from shutil import rmtree
from itertools import product
import pandas as pd
import re

class Parameter:
    '''
    COMSOL Parameter with attributes `name` `expr` and `unit`.

    `unit` attribute default set to be none.
    '''
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
    Use COMSOL for simulation and export the results to the specified folder.

    Parameters
    ----------
    mphpath : str
        The path to the COMSOL mph file.
    datapath : str
        The folder for exporting data files. Warning: Overwrite mode enabled.
    *parameters : Parameter objects, optional
        Other parameters that need to be manipulated (other than DC voltages).
        For each Parameter object, 
        - if its `expr` attribute is a number or string,
            the program will directly change it when performing simulations.
        - if its `expr` attribute is a list of numbers, list of strings or a numpy array,
            the program will iterate based on it.
        - if there is more than one Parameter object with `expr` attributes as a list or ndarray,
            the program will calculate their Cartesian product and iterate and simulate them sequentially.
    
    Returns
    -------
    CSV data files. Named after the electrodes, representing the data when
    this electrode is set to 1[V] while the other electrodes are set to 0[V].

    Examples
    --------
    The following code shows how to simulate a blade trap with different lengths (L) and blade widths (d).

    >>> from bemcomsol.simulation import simulate
    >>> simulate('model/blade.mph','data',Parameter('L',[100,200,400],'[um]'),Parameter('d',['0.8*L','0.9*L']))
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
        Read the data files obtained from the COMSOL simulation and convert them into numpy arrays.

        Parameters
        ----------
        path : str
            Directory path from which the program import CSV file.
        datalabel : str
            Column labels of the data, seperate by comma `','`.
            It should begin with coordinate variables, and then physics quantities.
        coordinates : int
            Number of coordinate variables. Default value is 3.

        Returns
        -------
        For each coordinate variable, it creates a numpy array to record the coordinate values,
        and the array is named as specified in `datalabel`.

        For each physical quantity, it creates a dictionary named as the specified name in `datalabel`
        followed by a suffix 's'. The keys of the dictionary are electrode names like `'DCi'`.
        The corresponding values are the values of the physical quantity when this electrode
        is set to 1[V] while the other electrodes are set to 0[V].

        Examples
        --------
        Assuming that after the simulation, we obtained two CSV files in the "data" folder,
        namely `DC1.csv` and `DC2.csv`. These files record the potential and the y-direction
        electric field when two DC electrodes are set to 1[V], with x=0, y coordinates at -1, 0, and 1.

        >>> # DC1.csv                 # DC2.csv
        ... +----+----+----+----+    +----+----+----+----+
        ... |x   |y   |V   |Ey  |    |x   |y   |V   |Ey  |
        ... +====+====+====+====+    +====+====+====+====+
        ... |0   |-1  |0   |-1  |    |0   |-1  |1   |1   |
        ... +----+----+----+----+    +----+----+----+----+
        ... |0   |-1  |1   |0   |    |0   |-1  |0   |0   |
        ... +----+----+----+----+    +----+----+----+----+
        ... |0   |-1  |0   |1   |    |0   |-1  |0   |0   |
        ... +----+----+----+----+    +----+----+----+----+

        We could use the following code to load them

        >>> from bemcomsol.simulation import load
        >>> [x,y,Vs,Eys] = load('data','x,y,V,Ey',coordinates=2)

        NOTE: The `datalabel` does not necessarily have to be the same as the variable names used in COMSOL;
        you are free to specify it as you wish.

        To access the data, use

        >>> print(x)
        >>> # x = [0,0,0]
        >>> print(Vs['DC2'])
        >>> # Eys['DC2'] = [1,0,0]
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
        newlines = [row for row in lines if row[0]!='%'] # Delete annotation generate by COMSOL
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