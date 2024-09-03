'''
Kinaid Dashboard

  A dashboard for analyzing and visualizing data from the phosphoproteomics experiments.
  Author: Javed M. Aman
  Date: 09/01/2024
  Version: 1.0
  License: MIT

'''

__author__ = 'Javed M. Aman'
__version__ = '1.0'
__license__ = 'MIT'
__email__ = 'javeda@princeton.edu'

from dash import Dash
import dash_bootstrap_components as dbc

#setup stylesheets
external_stylesheets = [dbc.themes.LITERA]
app = Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)

if __name__ == '__main__':
  app.run(host='0.0.0.0', port='8050', debug=True)