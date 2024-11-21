from setuptools import setup, find_packages

setup(
    name='kinaid',
    version='0.1',
    packages=['kinaid'],
    package_dir={'kinaid': 'kinaid'},
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
    entry_points={
        'console_scripts': [
            'kinaid-install = kinaid:run_default_configuration',
            'kinaid-add = kinaid:run_add_organism',
        ]
    }
)