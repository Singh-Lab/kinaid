from setuptools import setup, find_packages

setup(
    name='kinaid',
    version='1.0',
    find_packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'scikit-learn',
        'pandas',
        'tqdm',
        'tqdm-notebook',
        'plotly',
        'statsmodels',
        'dash_cytoscape',
        'mpire'
    ],
    author='Javed M. Aman',
)