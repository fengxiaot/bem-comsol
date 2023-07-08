import mph
import os
from globalvar import *
from shutil import rmtree

client = mph.start()
model = client.load(mph_file_path) # The path of mph file
if os.path.exists('data'):
    rmtree('data')
os.makedirs('data')
paramdict = model.parameters()

for param in paramdict:
    if param.startswith('DC'):
        for otherparam in paramdict:
            if otherparam.startswith('DC'):
                model.parameter(otherparam,'0[V]')
        model.parameter(param,'1[V]')
        print(model.parameters())
        model.solve()
        for node in model/'exports':
            node.property('filename',os.path.normpath(os.getcwd() +'/data/' + param + '.csv'))
            node.run()

model.clear()
model.reset()
client.clear()