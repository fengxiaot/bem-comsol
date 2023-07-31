from setuptools import setup

setup(
    name='bemcomsol',
    version='0.1.1',
    description='A python package that can simulate electrostatics for ion traps using boundary element method provided by COMSOL Multiphysics.',
    author='Xiaotian Feng',
    author_email='fxt2003@outlook.com',
    url='https://github.com/fengxiaot/bem-comsol',
    packages=['bemcomsol'],
    install_requires=[
        'numpy',
        'scipy',
        'mph',
        'pandas'
    ],
)
