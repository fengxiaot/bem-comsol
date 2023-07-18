import pandas as pd
import os

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
