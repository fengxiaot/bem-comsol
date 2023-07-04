import pandas as pd
import os

def load(dataLabel,param=3):
    '''
        Input
        - dataLabel: Data labels of CSV file. The elements of dataLabel should match the number of columns, and seperated by comma ','. It should begin with coordinate variables, and then physics quantity. Example: 'x,y,z,esbeV' or 'x,y,Ex,Ey'.
        - param: Number of coordinate variables. For example, 'x,y,Ex,Ey' should be set as 2. Default value is 3.
        Output: A list of object
        - For each coordinate variable, it creates a numpy array recording the data points
        - For each physics quantity, it creates a dict named with 'label' + 's', whose keys are 'DCi'.
    '''
    cwd = os.path.normpath(os.getcwd() +'/data/')
    files = os.listdir(cwd)
    DCnum = len(files)
    labels = dataLabel.split(',')
    for i in range(param,len(labels)):
        vars = labels[i] + 's'
        exec(f"{vars} = dict()")
    for file in files:
        with open(cwd + '/' + file, 'r') as COMSOL: 
            lines = COMSOL.readlines()
        newlines = [row for row in lines if row[0]!='%'] # Delete annotation
        with open('output.csv', 'w') as csvdata:
            csvdata.write(dataLabel+'\n')
            csvdata.writelines(newlines)
        data = pd.read_csv('output.csv')
        DCi = file.strip('.csv')
        for i in range(param,len(labels)):
            var = labels[i]
            vars = var + 's'
            exec(f"{vars}['{DCi}'] = data['{var}'].values")
    for i in range(param):
        coordinate = labels[i]
        exec(f"{coordinate} = data['{coordinate}'].values")
    datalist = []
    for i in range(len(labels)):
        if i < param:
            vars = labels[i]
        else:
            vars = labels[i] + 's'
        exec(f'datalist.append({vars})')
    return datalist
