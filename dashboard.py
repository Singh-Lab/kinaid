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

from dash import Dash, dcc, html
import dash_bootstrap_components as dbc
import uuid

#setup stylesheets
external_stylesheets = [dbc.themes.LITERA]
app = Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)

app.title = 'KINAID Dashboard'
logo_image = 'assets/logo3.png'

'''
Styling
'''

# the style arguments for the sidebar. We use position:fixed and a fixed width
SIDEBAR_STYLE = {
    'position': 'fixed',
    'top': 0,
    'left': 0,
    'bottom': 0,
    'width': '13rem',
    'padding': '1rem 1rem',
    'background-color': '#f8f9fa',
    'overflow': 'scroll'

}

# the styles for the main content position it to the right of the sidebar and
# add some padding.
CONTENT_STYLE = {
    "margin-top": "2rem",
    "margin-left": "20%",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
    "display": "inline-block",
    "width" : "60%"
}

'''
Sidebar
'''

sidebar = html.Div(
    [
        html.Img(src=logo_image, alt='logo', className='logo', style={'height':'50px', 'width':'100%'}),
        html.Hr(),  
        dbc.Nav(
            [
                dbc.Label('Home', color='black', style={'font-weight':'bold'}),
                dbc.NavLink('About', href='#about', active='exact', style={'padding-top':'5px'}, external_link=True, id='about-link'),
                dbc.NavLink('Settings', href='#settings', active='exact', style={'padding-top':'5px', 'padding-bottom':'15px'}, external_link=True, id='settings-link'),
                dbc.Label('Figures', color='black', style={'font-weight':'bold'}),
                dbc.NavLink('Matches Dataframe', href='#results-placeholder', active='exact', style={'padding-top':'5px'}, external_link=True, id='matches-link'),
                dbc.NavLink('Match Counts', href='#barplot-item', active='exact', style={'padding-top':'5px'}, external_link=True, id='counts-link'),
                dbc.NavLink('Peptide Volcano', href='#peptide-scatter-item', active='exact', style={'padding-top':'5px'}, external_link=True, id='peptide-scatter-link'),
                dbc.NavLink('Match Heatmap', href='#heatmap-item', active='exact', style={'padding-top':'5px'}, external_link=True, id='heatmap-link'),
                dbc.Label('Kinase Enrichment', color='black'),
                dbc.NavLink('Activty z-score', href='#zscore-item', active='exact', style={'padding-top':'5px'}, external_link=True, id='z-score-link'),
                dbc.NavLink('Activity p-value vs log2FC', href='#kinase_scatter-item', active='exact', style={'padding-top':'5px'}, external_link=True, id='kinase-scatter-link'),
                dbc.Label('Network', color='black'),
                dbc.NavLink('Kinase Hub', href='#hub-item', active='exact', style={'padding-top':'5px', 'padding-bottom':'15px'}, external_link=True, id='kinase-hub-link'),
                dbc.NavLink('Full Kinase Network', href='#network-item', active='exact', style={'padding-top':'5px', 'padding-bottom':'15px'}, external_link=True, id='kinase-network-link'),
            ],
            vertical=True,
            pills=True
        )     
    ],

    style=SIDEBAR_STYLE,
)

def serve_layout() :
  '''
  Build Layout
  '''
  session_id = str(uuid.uuid4())
  return html.Div([
      html.P(id='placeholder', style={'display':'none'}),
      dcc.Store(data=session_id, id='session-id'),
      dcc.Store(data=False, id='session-started'),
      #dcc.Store(data=None, id='upload-excel-store'),
      #dcc.Store(data=None, id='upload-data-store'),
      sidebar,
#      content
  ])
app.layout = serve_layout


if __name__ == '__main__':
  app.run(host='0.0.0.0', port='8060', debug=True)