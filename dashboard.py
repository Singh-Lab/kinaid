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



from dash import Dash, dcc, html, dash_table, Input, Output, State, no_update, ALL, MATCH, callback_context
import dash_bootstrap_components as dbc
import uuid
import pandas as pd
from kinaid import ortholog
from cachetools import TTLCache
from threading import Lock
import base64
import io
import random
import os
from kinaid.matching import PWM_Matrices,Scoring,PeptideBackground
from kinaid.session import Session
from icecream import ic
from typing import List
import dash_cytoscape as cyto
import zipfile
import argparse



#setup stylesheets
external_stylesheets = [dbc.themes.LITERA]
app = Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)

logo_image = 'assets/logo3.png'

example_filename = 'Leutert_KC_example.csv'
example_df = pd.read_csv(example_filename)

organism_scientific_names = {'human' : 'h. sapiens',
                    'mouse' : 'm. musculus',
                    'fly' : 'd. melanogaster',
                    'yeast' : 's. cerevisiae',
                    'worm' : 'c. elegans',
                    'zebrafish' : 'd. rerio'}

columns = ['peptide', 'log2fc', 'id', 'site', 'dependent']
placeholder_dict = {k:'<'+k.capitalize()+' Column>' for k in columns}
placeholder_dict['sheet'] = '<Sheet Name>'

orthologs_dir = 'orthologs'
data_dir = './data'
orthologs_suffix = '_orthologs_final.tsv'

default_organism = 'fly'
default_ambiguous = True

human_kinase_file = os.path.join(data_dir,'human_kinases_final.tsv')

ortholog_manager = ortholog.OrthologManager(orthologs_dir, orthologs_suffix, organisms=organism_scientific_names.keys(), human_kinase_file=human_kinase_file)
default_organism_orthologs = ortholog_manager.get_orthologs(default_organism)
default_gene_id = 'GeneID'
default_kinase_symbols = list(default_organism_orthologs.get_all_kinase_symbols_for_gene_id(default_gene_id, default_ambiguous))

default_network_match_threshold = 99
default_match_threshold = 90


johnson_ST_matrices_file = os.path.join(data_dir,'ST-Kinases.xlsx')
johnson_Y_matrices_file = os.path.join(data_dir,'Y-Kinases.xlsx')
densitometry_file = os.path.join(data_dir,'ST-Kinases_densitometry.xlsx')

ST_matrices = PWM_Matrices(johnson_ST_matrices_file)
ST_matrices_wfav = PWM_Matrices(johnson_ST_matrices_file)

ST_matrices_wfav.add_densitometry(densitometry_file)

Y_matrices = PWM_Matrices(johnson_Y_matrices_file)

_scoring_ = {'ST': Scoring(ST_matrices), 'Y' :Scoring(Y_matrices)}
_scoring_wfav = {'ST': Scoring(ST_matrices_wfav), 'Y' :Scoring(Y_matrices)}


st_background_file = os.path.join(data_dir, 'johnson_ochoa_background.tsv')
st_background_wfav_file = os.path.join(data_dir, 'johnson_ochoa_background_wfav.tsv')
y_background_file = os.path.join(data_dir, 'johnson_tyrosine_background.tsv')

_background_ = {'ST' : PeptideBackground(st_background_file), 'Y' : PeptideBackground(y_background_file)}
_background_wfav = {'ST' : PeptideBackground(st_background_wfav_file), 'Y' : PeptideBackground(y_background_file)}

cyto.load_extra_layouts()

'''
Helper Functions
'''

def determine_id_type_using_decision_tree(id, organism) :
    #determine id type using decision tree

    if(organism == 'human') :
        if(id.startswith('ENSP') or id.startswith('ENST') or id.startswith('ENSG')) :
            return 'Ensembl'
        elif(id[0].isdigit()) :
            return 'GeneID'
        else :
            return 'UniProtKB'
    elif(organism == 'mouse') :
        if(id.startswith('ENSMUSP') or id.startswith('ENSMUST') or id.startswith('ENSMUSG')) :
            return 'Ensembl'
        elif(not id[0].isdigit()) :
            return 'UniProtKB'
        elif(len(id) < 7 and id[0].isdigit()) :
            return 'GeneID'
        else :
            return 'MGI'
    elif(organism == 'fly') :
        if(id.startswith('FBgn')) :
            return 'FlyBase'
        elif(id[0].isdigit()) :
            return 'GeneID'
        else :
            return 'UniProtKB'
    elif(organism == 'yeast') :
        if(id.startswith('S00') or id.startswith('SGD')) :
            return 'SGD'
        elif(id[0].isdigit()) :
            return 'GeneID'
        else :
            return 'UniProtKB'
    elif(organism == 'worm') :
        if(id.startswith('WBGene')) :
            return 'WormBase'
        elif(id[0].isdigit()) :
            return 'GeneID'
        else :
            return 'UniProtKB'
    elif(organism == 'zebrafish') :
        if(id.startswith('ZDB-GENE')) :
            return 'ZFIN'
        elif(id[0].isdigit()) :
            return 'GeneID'
        else :
            return 'UniProtKB'
        
    return None

def guess_id_type(ids, organism, number_to_sample=10) :
    num_ids = len(ids)
    ids_to_sample = min(num_ids, number_to_sample)

    #sample ids_to_sample ids from ids
    sampled_ids = random.sample(ids, ids_to_sample)

    #determine id type for each sampled id and store counts in a dictionary
    id_type_counts = {}
    for id in sampled_ids :
        id_type = determine_id_type_using_decision_tree(id, organism)
        if id_type not in id_type_counts :
            id_type_counts[id_type] = 0
        id_type_counts[id_type] += 1
    
    #sort id_type_counts by value
    id_type_counts = {k: v for k, v in sorted(id_type_counts.items(), key=lambda item: item[1], reverse=True)}

    #return the id type with the highest count
    return list(id_type_counts.keys())[0]

def fill_common_names(value_dict, column_names) :
  '''
  match column names to common names
  probably can do it faster with loops or something
  '''
  value_dict_copy = value_dict.copy()
  column_names_copy = column_names.copy()
  #peptide
  if 'motif' in column_names_copy :
    value_dict_copy['peptide'] = 'motif'
    column_names_copy.remove('motif')
  elif 'Motif' in column_names_copy :
    value_dict_copy['peptide'] = 'Motif'
    column_names_copy.remove('Motif')
  elif 'Peptide' in column_names_copy :
    value_dict_copy['peptide'] = 'Peptide'
    column_names_copy.remove('Peptide')
  elif 'peptide' in column_names_copy :
    value_dict_copy['peptide'] = 'peptide'
    column_names_copy.remove('peptide')
    
  #log2fc
  if 'log2fc' in column_names_copy :
    value_dict_copy['log2fc'] = 'log2fc'
    column_names_copy.remove('log2fc')
  elif 'log2FC' in column_names_copy :
    value_dict_copy['log2fc'] = 'log2FC'
    column_names_copy.remove('log2FC')
  elif 'Log2FC' in column_names_copy :
    value_dict_copy['log2fc'] = 'Log2FC'
    column_names_copy.remove('Log2FC')
  elif 'Log2fc' in column_names_copy :
    value_dict_copy['log2fc'] = 'Log2fc'
    column_names_copy.remove('Log2fc')
  #id
  if 'id' in column_names_copy :
    value_dict_copy['id'] = 'id'
    column_names_copy.remove('id')
  elif 'ID' in column_names_copy :
    value_dict_copy['id'] = 'ID'
    column_names_copy.remove('ID')
  elif 'Id' in column_names_copy :
    value_dict_copy['id'] = 'Id'
    column_names_copy.remove('Id')
  elif 'Uniprot' in column_names_copy :
    value_dict_copy['id'] = 'Uniprot'
    column_names_copy.remove('Uniprot')
  elif 'uniprot' in column_names_copy :
    value_dict_copy['id'] = 'uniprot'
    column_names_copy.remove('uniprot')

  #site
  if 'site' in column_names_copy :
    value_dict_copy['site'] = 'site'
    column_names_copy.remove('site')
  elif 'Site' in column_names_copy :
    value_dict_copy['site'] = 'Site'
    column_names_copy.remove('Site')
  elif 'position' in column_names_copy :
    value_dict_copy['site'] = 'position'
    column_names_copy.remove('position')
  elif 'Position' in column_names_copy :
    value_dict_copy['site'] = 'Position'
    column_names_copy.remove('Position')
  elif 'Site Position' in column_names_copy :
    value_dict_copy['site'] = 'Site Position'
    column_names_copy.remove('Site Position')
  elif 'site position' in column_names_copy :
    value_dict_copy['site'] = 'site position'
    column_names_copy.remove('site position')

  #pvalue
  if 'pvalue' in column_names_copy :
    value_dict_copy['dependent'] = 'pvalue'
    column_names_copy.remove('pvalue')
  elif 'p-value' in column_names_copy :
    value_dict_copy['dependent'] = 'p-value'
    column_names_copy.remove('p-value')
  elif 'P-value' in column_names_copy :
    value_dict_copy['dependent'] = 'P-value'
    column_names_copy.remove('P-value')
  elif 'Pvalue' in column_names_copy :
    value_dict_copy['dependent'] = 'Pvalue'
    column_names_copy.remove('Pvalue')

  #get all None value items from value_dict_copy
  none_values = [k for k,v in value_dict_copy.items() if v is None]


  for k in none_values :
    #if there are no more column names left, break
    if(len(column_names_copy) == 0) :
      break
    #place columns in order of appearance
    value_dict_copy[k] = column_names_copy[0]
    column_names_copy.remove(column_names_copy[0])
  
    
  return value_dict_copy

def process_df_columns(df, organism) :
  '''
  Figure out column options from dataframe
  '''
  column_names = list(df.columns)
  column_options = [{'label': col, 'value': col} for col in column_names]

  all_column_options = [[{'label': placeholder_dict[c], 'value': None}] + column_options for c in columns]
  disable_dropdown_dict = {k:False for k in columns}

  value_dict = {k:None for k in columns}
  value_dict = fill_common_names(value_dict, column_names)

  id_type_value = None
  if(value_dict['id'] is not None) :
    #make sure id column is string
    df[value_dict['id']] = df[value_dict['id']].astype(str)
    ids = df[value_dict['id']].tolist()
    id_type_value = guess_id_type(ids, organism)
  out = (
    all_column_options,
    [value_dict[c] for c in columns],
    [disable_dropdown_dict[c] for c in columns],
    id_type_value
  )
  return out


def write_archive(bytes_io, df_list, fileNameList):
  with zipfile.ZipFile(bytes_io, mode='w') as zf:
      for df,fn in zip(df_list,fileNameList):
          df_bytes=df.to_csv(index=False).encode('utf-8')
          zf.writestr(fn ,df_bytes)
  
'''
Create Cache
'''
timeout = 60 * 60 # 60 minutes
max_sessions = 100
cache = TTLCache(max_sessions, ttl=timeout)

cache_lock = Lock()

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
About Section
'''
about = html.Div([
        html.H2('KINAID', style={'text-align':'center'}),
        html.Hr(),
        html.P('This is a web application for the analysis of phosphoproteomics data. It uses the matrices from Johnson et al. (2023) and Yaron-Barir et al. (2024) to perform the matching between kinsases and substrates.', style={'text-align':'left'}),
        html.P('Refreshing or closing this page will clear and reset the data, generating a new session. While page is loaded, data will remain on the server at most 60 minutes since last activity', style={'font-weight':'bold'}),
        html.P('Please reference:', style={'text-align':'center', 'font-style':'oblique'}),
        html.P('Aman J, Zhu A, Wühr M, Shvartsman S, Singh M. (2024)', style={'text-align':'center', 'font-weight':'bold', 'font-style':'oblique'}),
        html.A(href='https://github.com/Singh-Lab/kinaid', children=[
          html.P('GitHub', style={'margin-left':'47%', 'font-weight':'bold', 'clickable':'true'})
        ]),
    ],
    id='about'
)

'''
Download example
'''
example = html.Div([
        html.P('The example is from Leutert et al. 2023:the osmotic stress condition in KC and expected to perturb the HOG pathway'),
        dbc.Accordion([
            dbc.AccordionItem(
                html.Div(
                    children = [
                      dash_table.DataTable(
                          data=example_df.to_dict('records'),
                          columns=[{'name': i, 'id': i} for i in example_df.columns],
                          id='example-table',
                          export_format='csv',
                          
                          style_table={
                              'overflowX': 'scroll',
                              'maxHeight': '500px',
                              },
                          style_cell={
                              'minWidth': '100px', 'width': '100px', 'maxWidth': '100px',
                              'overflow': 'hidden',
                              'textOverflow': 'ellipsis',
                              'textAlign': 'center'
                          },
                          style_header={
                              'backgroundColor': 'rgb(230, 230, 230)',
                              'fontWeight': 'bold',
                              'textAlign': 'center'
                          },
                          style_data_conditional=[{
                              'if': {'row_index': 'odd'},
                              'backgroundColor': 'rgb(248, 248, 248)'
                          }],
                      )
                    ],
                ),
                title='Example Data: Leutert et al. 2023',
            )
        ],
        start_collapsed=True, 
      ),
    ])

'''
Settings Section
'''

settings = html.Div(children=[
        html.H2('Settings', style={'text-align':'center'}),
        html.Hr(),
        dbc.Label('Organism (select first)', style={'margin-left': '1rem', 'margin-top': '2rem'}),
        html.Div([
            dbc.RadioItems(
                id='organism-radioitems',

                options=[{'label': o, 'value': s} for s,o in organism_scientific_names.items()],
                value=None,
                labelStyle={'display': 'inline-block'},
            style={

                'margin-left': '2rem',
            },
                inline=True
            ),
        ]),
        
        html.Hr(),
        html.P('Upload your phosphoproteomics data here (.csv, .tsv, .xlsx)', style={'text-align':'center'}),
                #choose organism to for id mapping
        dbc.Row(
            [
            dbc.Col(
                dbc.Input(
                    id='upload-path',
                    placeholder='{.csv|.tsv|.xlsx}',
                    type='text',
                    value='',
                    disabled=True
                ),
                style={'width':'100%', 'height':'3rem', 'padding-left':'10rem'},
            ),
            dbc.Col(
                  dcc.Upload(
                    dbc.Button('Upload File', color='primary', className='mr-1',
                                style={'width':'8rem',
                                       'font-weight':'bold',
                                       'height' : '2.8rem',
                                       'font-size':'14px',
                                       'vertical-align':'top'},
                                disabled=True,
                                id='upload-button'
                              ),
                                       
                    id='upload-data',
                    multiple=False,
                    disabled=True
                  ),
                  width={'size': 2, 'order': '2'}
            ),
            dbc.Col(
              dbc.Button('Example', color='info', className='mr-1',
                          style={'width':'8rem',
                                'font-weight':'bold',
                                'height' : '2.8rem',
                                'font-size':'14px',
                                'vertical-align':'top'},
                          id='example-button'
                        ),
                  width={'size': 2, 'order': '3'},
                  
            ),
            ],
            justify='center',
            className='g-0',
        ),
        dbc.Label('Sheet Name (Excel only)', style={'margin-left': '1rem', 'margin-top': '2rem'}),
        dbc.Select(
            id='sheet-dropdown',
            placeholder='Sheet',
            style={
                'width': '20rem',
                'margin-left': '2rem',
            },
            disabled=True,
            value=None
        ),
        dbc.Label('Column with peptides [Required]', style={'margin-left': '1rem', 'margin-top': '2rem'}),
        dbc.Select(
            id={'type': 'column-dropdown', 'index': 'peptide'},
            placeholder=placeholder_dict['peptide'],
            style={
                'width': '20rem',
                'margin-left': '2rem',
            },
            disabled=True,
        ),     
        dbc.Label('Column with log2FC [Required for enrichment analyses]', style={'margin-left': '1rem', 'margin-top': '2rem'}),
        dbc.Select(
            id={'type': 'column-dropdown', 'index': 'log2fc'},
            placeholder=placeholder_dict['log2fc'],
            style={
                'width': '20rem',
                'margin-left': '2rem',
            },
        disabled=True,
        ),
        dbc.Row([
            dbc.Col([
              dbc.Label('Column with ids [Required for heatmap and network]', style={'margin-left': '1rem', 'margin-top': '2rem'}),
              dbc.Select(
                  id={'type': 'column-dropdown', 'index': 'id'},
                  placeholder=placeholder_dict['id'],
                  style={
                      'width': '20rem',
                      'margin-left': '2rem',
                  },
                  disabled=True,
              ),
            ]),
            dbc.Col([
                dbc.Label('Type of ID', style={'margin-top': '2rem'}),
                dbc.Select(
                    id='id-type-dropdown',
                    placeholder='UniProtKB',
                    style={
                        'width': '20rem',
                        'margin-left': '1rem',
                    },
                ),
            ]),
        ]),
        dbc.Label('Column with site positions [Optional for naming peptides]', style={'margin-left': '1rem', 'margin-top': '2rem'}),
        dbc.Select(
            id={'type': 'column-dropdown', 'index': 'site'},
            placeholder=placeholder_dict['site'],
            style={
                'width': '20rem',
                'margin-left': '2rem',
            },
            disabled=True,
        ),

        dbc.Label('Column with mass-spec p-value [Optional for volcano]', style={'margin-left': '1rem', 'margin-top': '2rem'}),
        dbc.Select(
            id={'type': 'column-dropdown', 'index': 'dependent'},
            placeholder=placeholder_dict['dependent'],
            style={
                'width': '20rem',
                'margin-left': '2rem',
            },
        disabled=True,
        ),
        dbc.Accordion([
            dbc.AccordionItem(
                html.Div(
                    children = [
                        dbc.Row([
                            dbc.Col([
                                #toggle checkmark to use favorability or not
                                dbc.Label('Use S/T favorability?'),
                                html.Div([
                                    dbc.Checklist(
                                        id='favorability-check',
                                        options=[
                                            {'label': 'Use favorability (Recommended)', 'value': 'favorability'},
                                        ],
                                        style={
                                            'margin-left': '10px'
                                        },
                                    value=['favorability']
                                    )
                                ]),
                                dbc.Label('Allow ambiguous kinase mappings'),
                                html.Div([
                                    dbc.Checklist(
                                        id='ambiguous-check',
                                        options=[
                                            {'label': 'Yes (Recommended)', 'value': 'ambiguous'},
                                        ],
                                        style={
                                            'margin-left': '10px'
                                        },
                                        value=['ambiguous'] if default_ambiguous else [],
                                    )
                                ]),
                                #get threshold for ds
                                dbc.Label(f'Match threshold (default {default_match_threshold})'),

                                html.Div(
                                    [dbc.Input(
                                        id='threshold-input',
                                        type='number',
                                        placeholder=str(default_match_threshold),
                                        value=default_match_threshold,
                                        style={
                                            'margin-left': '10px',
                                            'width' : '100px',
                                            'margin-bottom': '10px'
                                        }
                                    ),]
                                ),
                                dbc.Label('Subset of Kinases'),
                                dcc.Upload(
                                    dbc.Button('Upload list of kinases', color='primary', className='mr-1',
                                                style={'width':'8rem',
                                                    'font-weight':'bold',
                                                    'height' : '3.5rem',
                                                    'font-size':'14px',
                                                    'vertical-align':'middle',
                                                    'margin-top': '10px',}),
                                    id='upload-kinases',
                                    multiple=False
                                )

                            ]),
                            dbc.Col([
                                dbc.Label('Kinases to match'),
                                html.Div([
                                    dcc.Dropdown(id='kinase-selection',
                                        options=[{'label': k, 'value': k} for k in default_kinase_symbols],
                                        multi=True,
                                        clearable=True,
                                        value=default_kinase_symbols,
                                        style={
                                            'maxHeight': '500px',
                                            'overflowY': 'scroll',
                                        },
                                        disabled=True
                                    )
                                                    
                                ])
                            ])
                        ])
                    ]
                ),
                title='Advanced Options',
            ),
        ],
        style={
            'margin': '10px',
            'width' : '65%'
        },
        start_collapsed=True,
        ),
        dbc.Button([dbc.Spinner(children= html.Div(id = 'loading'), size='sm'), 'Submit'],
                    size='lg', id='submit-button', color='primary', className='mr-1',
                    style={'width':'8rem',
                            'font-weight':'bold',
                            'height' : '2.8rem',
                            'font-size':'20px',
                            'margin-left':'47%',}),
    ],
    id='settings',
    style={'border': '1px solid #d3d3d3', 'border-radius': '5px', 'padding': '10px'}
)

'''
Time Out
'''
timeout_dialog = html.Div([
  dcc.ConfirmDialog(
    id='timeout-popup',
    message='Your session has timed out. Please refresh the page to start a new session.',
    displayed=False
    )
])

'''
Result Table
'''
result_table =  dash_table.DataTable(
    id='output-table',
    export_format='csv',
    columns=[{'name': i, 'id': i} for i in ['Peptide', 'log2FC', 'Kinase(s)']],
    style_table={
        'overflowX': 'scroll',
        'maxHeight': '500px',
        'overflowY': 'scroll'
        },
    style_cell={
        'minWidth': '100px', 'width': '100px', 'maxWidth': '100px',
        'overflow': 'hidden',
        'textOverflow': 'ellipsis',
        'textAlign': 'center'
    },
    style_header={
        'backgroundColor': 'rgb(230, 230, 230)',
        'fontWeight': 'bold',
        'textAlign': 'center'
    },
    style_data_conditional=[{
        'if': {'row_index': 'odd'},
        'backgroundColor': 'rgb(248, 248, 248)'
    }],
)

'''
Download Section
'''
download_section = dbc.Row([
          dbc.Col([
            dbc.Label('Download .csv files:', style={'font-weight':'bold', 'margin-left': '50%'}),
            dbc.RadioItems(
              id='download-radioitems',
              options=[
                {'label': 'Download percentiles', 'value': 'percentiles'},
                {'label': 'Download network', 'value': 'network'},
                {'label': 'Download kinase network', 'value': 'kinase_network'},
                {'label': 'Download kinase statistics', 'value': 'kinase_stats'},
              ],
              value='scores',
              labelStyle={'display': 'block'},
              style={
                'margin-left': '50%'
              },
            )
          ]),
          dbc.Col([
              dbc.Button('Download', size='lg', id='download-button', color='primary', className='mr-1',
                         style={'width':'8rem',
                            'font-weight':'bold',
                            'height' : '2.8rem',
                            'font-size':'20px',
                            'margin-left':'47%',}),
              dcc.Download(id='download-data')
          ]),
        ])

downloads = html.Div([
    html.H2('Download', style={'text-align':'center', 'margin-top':'10rem'}),
    html.Hr(),
    html.Div(
      children = [
        html.H4('(submit data to enable)', style={'text-align':'center'}),
      ],
      id='download-placeholder',
    )
])

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

figures = html.Div([
    html.H2('Figures', style={'text-align':'center', 'margin-top':'10rem'}),
    html.Hr(),
    html.Div(
      children = [
        html.H4('(submit data to see results)', style={'text-align':'center'}),
      ],
      id='results-placeholder',
    ),
    dbc.Accordion([
      dbc.AccordionItem(
        dbc.Spinner(
          children = [
            dcc.Graph(
              id = {'type': 'figure', 'index': 'barplot'},
              ),
            ],
        ),
        title='Match Counts',  
        id='barplot-item'
      )
      ],
      id={'type': 'figure-accordion', 'index': 'barplot'},
      start_collapsed=True
    ),
    dbc.Accordion([
      dbc.AccordionItem(   
        children = [
          dbc.Spinner(
            children = [
              dcc.Graph(
                id = {'type': 'figure', 'index': 'peptide_scatter'}
              ),
            ],
          ),
          html.Div(
            [
              html.H4('Options', style={'text-align':'center'}),
              dbc.Label('Kinase(s) to highlight'),
              html.Div([
                dcc.Dropdown(
                  id='scatter-kinase-selection',
                  options=[{'label': k, 'value': k} for k in default_kinase_symbols],
                  multi=True,
                  clearable=True,
                  value=[],
                  style={
                    'maxHeight': '3rem',
                    #'overflowY': 'scroll',
                  }
                )                   
              ])
            ],
            style={'border': '1px solid #d3d3d3', 'border-radius': '5px', 'padding': '30px', 'background-color': '#f4f4f4'}
          ),
        ],

        title='Peptide Volcano',
        id='peptide-scatter-item'      
      ),
    ],
      id={'type': 'figure-accordion', 'index': 'peptide_scatter'},
      start_collapsed=True,
    ),
    dbc.Accordion(
        [
          dbc.AccordionItem(
            dbc.Spinner(
              children = [
                dcc.Graph(
                  id = {'type': 'figure', 'index': 'heatmap'}
                ),
              ],
            ),
            title='Match Heatmap',
            style={
              'overflowY': 'scroll',
              'overflowX': 'scroll'
            },
            id='heatmap-item'
          )
        ],
        id={'type': 'figure-accordion', 'index': 'heatmap'},
        start_collapsed=True
    ),
])


enrichment_analysis = html.Div([
  dbc.Accordion([
    dbc.AccordionItem(   
      children = [
        html.Div(
          [
            html.H4('Options', style={'text-align':'center'}),
            dbc.Label('FDR threshold', style={'margin-left': '9rem'}),
            dbc.Input(
              id='fdr-input',
              type='number',
              value=0.05,
              min=0,
              max=1,
              step=0.01,
              style={
                'margin-left': '10rem',
                'width' : '100px'
              }
            ),
          ],
          style={'border': '1px solid #d3d3d3', 'border-radius': '5px', 'padding': '10px', 'background-color': '#f4f4f4'}
        ),
        dbc.Spinner(
          children = [
            dcc.Graph(
              id = {'type': 'figure', 'index': 'zscore'}
            ),
          ],
        ),
      ],
      title='Kinase Activity z-score',
      id='zscore-item'   
    ),
  ],
  id={'type': 'figure-accordion', 'index': 'zscore'},
  start_collapsed=True,
  ),
  dbc.Accordion([
    dbc.AccordionItem([
      'Scatter plot of log2FC vs -log10(adj p-value)',
      dbc.Spinner(
        children = [
          dcc.Graph(
            id = {'type': 'figure', 'index': 'kinase_scatter'}
          )
        ],
      )
      ],
      title='Kinase Activity p-value vs log2FC',
      id='kinase_scatter-item'
    )
    ],
    id={'type': 'figure-accordion', 'index': 'kinase_scatter'},
    start_collapsed=True
  ),


  
])
network_analysis = html.Div([
  #html.H2('Data Analysis', style={'text-align':'center', 'margin-top':'10rem'}),
  #html.Hr(),
  dbc.Accordion([ 
    dbc.AccordionItem(
      children = [
        html.Div(
        [
          html.H4('Options', style={'text-align':'center'}),
          dbc.Row([
            dbc.Col([
              dbc.Label('Kinase(s) to generate hub'),
              html.Div([
                dcc.Dropdown(
                  id='hub-kinase-selection',
                  options=[{'label': k, 'value': k} for k in default_kinase_symbols],
                  multi=True,
                  clearable=True,
                  value=[],
                  style={
                    'maxHeight': '3rem',
                    #'overflowY': 'scroll',
                  }
                )                   
              ])
            ]),
            dbc.Col([
              dbc.Label('match threshold', style={'margin-left': '9rem'}),
              dbc.Input(
                id='hub-match-threshold',
                type='number',
                value=default_match_threshold,
                min=0,
                max=100,
                step=0.5,
                style={
                  'margin-left': '10rem',
                  'width' : '100px',
                  'height' : '3rem'
                },
                debounce=True
              ),
            ]),
            dbc.Col([
              dbc.Label('Show only kinases?'),
              dbc.Checklist(
                id='hub-only-kinases',
                options=[
                  {'label': 'Yes', 'value': 'only_kinases'},
                ],
                style={
                  'margin-left': '10px'
                },
                value=['only_kinases']
              )
            ])
          ]),

          # dbc.Row([
          #   dbc.Col([
          #     #reload button centered
          #     dbc.Button('Reload', size='lg', id='hub-reload-button', color='primary', className='mr-1',
          #                 style={'width':'8rem',
          #                       'font-weight':'bold',
          #                       'height' : '2.8rem',
          #                       'font-size':'20px',
          #                       'margin-left':'45%',
          #                       'margin-top': '10px'}),
          #   ])
          # ])
          
        ],
        style={'border': '1px solid #d3d3d3', 'border-radius': '5px', 'padding': '30px', 'background-color': '#f4f4f4'}
        ),
        dbc.Spinner(
          children = [
            dcc.Graph()
          ],
          id='hub-spinner'
        ),
      ],
      title='Kinase Hub',
      style={
        'overflowY': 'scroll',
        'overflowX': 'scroll'            
      },
      id='hub-item',
    ),
  ],
  #id={'type': 'figure-accordion', 'index': 'hub'},
  id={'type': 'cyto-accordion', 'index': 'hub'},
  start_collapsed=True  
  ),
  dbc.Accordion([
    dbc.AccordionItem(children = [
      html.Div(
        [
          html.H4('Options', style={'text-align':'center'}),
          dbc.Label('match threshold', style={'margin-left': '9rem'}),
          dbc.Input(
            id='network-match-threshold',
            type='number',
            value=default_network_match_threshold,
            min=0,
            max=100,
            step=0.5,
            style={
              'margin-left': '10rem',
              'width' : '100px'
            },
            debounce=True
          ),
        ],
        style={'border': '1px solid #d3d3d3', 'border-radius': '5px', 'padding': '10px', 'background-color': '#f4f4f4'}
      ),
      'Network plot of kinases and their targets',
      dbc.Spinner(
        id = 'network-spinner',
        children = [dcc.Graph()],
      )
    ],
    title='Full Kinase Network',
    id = 'network-item'
    )
  ],
  id={'type': 'cyto-accordion', 'index': 'network'},
  start_collapsed=True
  ),
])


tooltips = html.Div([
  dbc.Tooltip(
    'A barplot of the number of matches for each kinase',
    target={'type': 'figure-accordion', 'index': 'barplot'},
    placement='top',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'A volcano plot of the mass-spec p-value vs the log2FC in phosphorylation',
    target={'type': 'figure-accordion', 'index': 'peptide_scatter'},
    placement='top',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'A heatmap of the matches between kinases and peptides',
    target={'type': 'figure-accordion', 'index': 'heatmap'},
    placement='top',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'A bar plot of the z-scores of the log2FCs for each kinase',
    target={'type': 'figure-accordion', 'index': 'zscore'},
    placement='top',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'A scatter plot of the -log adjusted p-value (of z-scores) vs mean log2FC of the matches for each kinase',
    target={'type': 'figure-accordion', 'index': 'kinase_scatter'},
    placement='top',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'Generate a network piecewise around selected kinases',
    target={'type': 'cyto-accordion', 'index': 'hub'},
    placement='top',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'A network of all kinases and their targets',
    target={'type': 'cyto-accordion', 'index': 'network'},
    placement='top',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'The name of the sheet in the Excel file that contains the phosphoprotemics data',
    target='sheet-dropdown',
    placement='right',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'The name of the column in the data set that contains the -5/+4 10-mer of the phosphorylation site',
    target={'type': 'column-dropdown', 'index': 'peptide'},
    placement='right',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'The name of the column that contains the change in phosphroylation of the site between conditions (usually log base 2)',
    target={'type': 'column-dropdown', 'index': 'log2fc'},
    placement='right',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'The name of the column that contains the identifier of the protein of the site',
    target={'type': 'column-dropdown', 'index': 'id'},
    placement='right',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'The name of the column that contains the position of the phosphorylation site',
    target={'type': 'column-dropdown', 'index': 'site'},
    placement='right',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'The name of the column that contains the mass-spec p-value of the phosphorylation site (or any other dependent variable)',
    target={'type': 'column-dropdown', 'index': 'dependent'},
    placement='right',
    style={'font-size':'1rem'}
  ),
  dbc.Tooltip(
    'The organism specific identifier of the proteins in the dataset (e.g. UniProtKB)',
    target='id-type-dropdown',
    placement='left',
    style={'font-size':'1rem'}
  ),

])

'''
All Content
'''
content= html.Div(
    [
       about,
       example,
       settings,
       timeout_dialog,
       downloads,
       figures,
       enrichment_analysis,
       network_analysis,
       tooltips
    ],
    id='main-content',
    style=CONTENT_STYLE
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
      content
  ])

'''
Callbacks
'''

@app.callback([Output('kinase-selection', 'value', allow_duplicate=True),
                Output('kinase-selection', 'options', allow_duplicate=True),
                Output('id-type-dropdown', 'options'),
                #Output('id-type-dropdown', 'value', allow_duplicate=True),
                Output('upload-button', 'disabled'),
                Output('upload-data', 'disabled'),
                Output('kinase-selection', 'disabled', allow_duplicate=True)],
              Input('organism-radioitems', 'value'),
              State('ambiguous-check', 'value'),
              prevent_initial_call=True

)
def click_organism(organism : str, ambiguous : List[str]):
  '''
  Callback for clicking on an organism
    
  Update the kinase selection dropdown and id types dropdown
  '''
  print(callback_context.triggered_prop_ids)
  available_id_types = list(ortholog_manager.get_orthologs(organism).get_available_id_types())
  symbol_names = list(ortholog_manager.get_orthologs(organism).get_all_kinase_symbols_for_gene_id(available_id_types[0], ambiguous[0] == 'ambiguous' if len(ambiguous) > 0 else False))
  symbol_options = [{'label': k, 'value': k} for k in symbol_names]
  return (symbol_names, symbol_options, available_id_types, False, False, False)


@app.callback([
                Output('kinase-selection', 'value', allow_duplicate=True),
                Output('kinase-selection', 'options', allow_duplicate=True),
              ],
                Input('upload-kinases', 'contents'),
              [
                State('upload-kinases', 'filename'),
                State('organism-radioitems', 'value'),
                State('ambiguous-check', 'value'),
                State('id-type-dropdown', 'value'),
              ],
              prevent_initial_call=True
)
def upload_kinases(contents : str, filename : str, organism : str, ambiguous : List[str], id_type : str) :
  '''
  Callback for uploading kinases
  '''
  _, content_string = contents.split(',')
  decoded = base64.b64decode(content_string)

  print('loaded file: %s' % filename)

  current_symbol_names_list = list(ortholog_manager.get_orthologs(organism).get_all_kinase_symbols_for_gene_id(id_type,  ambiguous[0] == 'ambiguous' if len(ambiguous) > 0 else False))

  current_symbol_names = set(current_symbol_names_list)
  #expect txt file

  if not filename.endswith('.txt') :
    print('File must be .txt')
    return no_update
  new_kinase_names = [l.strip() for l in decoded.decode('utf-8').split('\n') if l.strip() != '']

  new_kinase_names = set(new_kinase_names)

  #get kinases that in kinase_names but not in current_organism_kinases
  kinase_names_diff = new_kinase_names.difference(current_symbol_names)

  print('Missing kinases: %s' % ', '.join(kinase_names_diff))

  if (len(new_kinase_names) == len(kinase_names_diff)) :
    print('No new kinases so reverting to original list')
    return no_update
  new_kinase_names = new_kinase_names.intersection(current_symbol_names)
  kinase_names = list(new_kinase_names)

  kinase_names.sort()

  kinase_options = [{'label': k, 'value': k} for k in kinase_names]

  return kinase_names, kinase_options

@app.callback(
   [Output('sheet-dropdown', 'disabled'),
    Output('sheet-dropdown', 'options'),
    Output('sheet-dropdown', 'value'),
    Output({'type': 'column-dropdown', 'index': ALL}, 'options', allow_duplicate=True),
    Output({'type': 'column-dropdown', 'index': ALL}, 'value', allow_duplicate=True), 
    Output({'type': 'column-dropdown', 'index': ALL}, 'disabled', allow_duplicate=True),
    Output('id-type-dropdown', 'value', allow_duplicate=True),
    Output('upload-path', 'value')],
   Input('upload-data', 'contents'),
   [State('upload-data', 'filename'),
    State('organism-radioitems', 'value')],
    prevent_initial_call=True,
)
def upload_button(contents, filename, organism) :
  '''
  Callback for uploading data
  '''
  if contents is None :
    return no_update
  if not any(filename.endswith(ext) for ext in ['.csv', '.tsv', '.xlsx']) :
    return no_update
  _, content_string = contents.split(',')
  decoded = base64.b64decode(content_string)

  df = None
  sheet_disabled = True
  sheet_options = []

  if filename.endswith('.csv') :
    df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
    sheet_value = None
  elif filename.endswith('.tsv') :
    df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t')
    sheet_value = None
  elif filename.endswith('.xlsx') :
    excel_file = pd.ExcelFile(io.BytesIO(decoded))
    sheet_names = excel_file.sheet_names
    sheet_options = [{'label': sheet_name, 'value': sheet_name} for sheet_name in sheet_names]
    sheet_value = sheet_names[0]
    df = pd.read_excel(excel_file, sheet_name=sheet_value)
    sheet_disabled = False

  columns_out = process_df_columns(df, organism)
  out = (
    sheet_disabled,
    sheet_options,
  sheet_value) + columns_out

  out += tuple([filename])
    
  return out

@app.callback([
                Output({'type': 'column-dropdown', 'index': ALL}, 'options', allow_duplicate=True),
                Output({'type': 'column-dropdown', 'index': ALL}, 'value', allow_duplicate=True),
                Output({'type': 'column-dropdown', 'index': ALL}, 'disabled', allow_duplicate=True),
                Output('id-type-dropdown', 'value', allow_duplicate=True)],
              Input('sheet-dropdown', 'value'),
              [
                State('upload-data', 'contents'),
                State('upload-data', 'filename'),
                State('organism-radioitems', 'value')
              ],
              prevent_initial_call=True)
def change_sheet(sheet_name, excel_data, filename, organism) :
  '''
  Callback for changing sheet
  '''
  if not filename.endswith('.xlsx') :
    return no_update
  if excel_data is None :
    return no_update
  _, content_string = excel_data.split(',')
  decoded = base64.b64decode(content_string)
  excel_file = pd.ExcelFile(io.BytesIO(decoded))
  df = pd.read_excel(excel_file, sheet_name=sheet_name)
  columns_out = process_df_columns(df, organism)

  return columns_out


@app.callback([
                Output('loading', 'children'),
                Output('results-placeholder', 'children', allow_duplicate=True),
                Output('scatter-kinase-selection', 'options', allow_duplicate=True),
                Output('hub-kinase-selection', 'options', allow_duplicate=True),
                Output('hub-match-threshold', 'value', allow_duplicate=True),
                Output('download-placeholder', 'children'),
                Output('session-started', 'data')
              ],
                Input('submit-button', 'n_clicks'),
              [
                State('session-id', 'data'),
                State('organism-radioitems', 'value'),
                State('upload-data', 'contents'),
                State('upload-data', 'filename'),
                State('sheet-dropdown', 'value'),
                State('kinase-selection', 'value'),
                State({'type': 'column-dropdown', 'index': ALL}, 'id'),
                State({'type': 'column-dropdown', 'index': ALL}, 'value'),
                State('id-type-dropdown', 'value'),
                State('threshold-input', 'value'),
                State('ambiguous-check', 'value'),
                State('favorability-check', 'value'),
                State('sheet-dropdown', 'disabled'),
              ],
              prevent_initial_call=True
)
def press_submit(
    n_clicks : int,
    session : str,
    organism : str,
    filedata : str,
    filename : str,
    sheet_name : str,
    kinase_names : list,
    column_ids : list,
    column_values : list,
    id_type : str,
    threshold : float,
    ambiguous : List[str],
    favorability : str,
    sheet_disabled : bool
    ) :
    '''
    Callback for pressing submit
    '''
    if n_clicks is None or n_clicks == 0 :
        print('How did you get here? No clicks!')
        return no_update

    if cache.currsize > max_sessions :
        print('Users maxed out')
        return no_update
    
    if filedata is None :
        return no_update
  
    _, content_string = filedata.split(',')
    decoded = base64.b64decode(content_string)
    
    #determine if excel or not
    #if sheet disabled and file ends with .xlsx, then it is an excel file
    if not sheet_disabled and filename.endswith('.xlsx') :
        excel_file = pd.ExcelFile(io.BytesIO(decoded))
        df = pd.read_excel(excel_file, sheet_name=sheet_name)
    else :
        if filename.endswith('.csv') :
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
        elif filename.endswith('.tsv') :
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')), sep='\t')
    
    #column_id to column index dict
    column_id_to_index = {c['index']:i for i,c in enumerate(column_ids)}
    
    #get column names
    column_names_dict = {c:column_values[i] for c,i in column_id_to_index.items()}

    if (len(favorability) > 0) and (favorability[0] == 'favorability') :
      print('Using favorability')
      kinaid_session = Session(session, organism, df, column_names_dict, _scoring_wfav, _background_wfav, ortholog_manager, selected_symbols=kinase_names, id_type=id_type, ambiguous=(ambiguous[0] == 'ambiguous' if len(ambiguous) > 0 else False), debug=True)
    else :
      print('Not using favorability')
      kinaid_session = Session(session, organism, df, column_names_dict, _scoring_, _background_, ortholog_manager, selected_symbols=kinase_names, id_type=id_type, ambiguous=(ambiguous[0] == 'ambiguous' if len(ambiguous) > 0 else False), debug=True)
    
    cache_lock.acquire()
    cache[session] = kinaid_session
    cache_lock.release()
    
    table_df = kinaid_session.get_kinase_matches_df()
    result_table.data = table_df.to_dict('records')
    result_table.columns = [{'name': i, 'id': i} for i in table_df.columns]
    result_table.style_data_conditional = [{
      'if': {'row_index': 'odd'},
      'backgroundColor': 'rgb(248, 248, 248)'
    }]
    
    return ('Finished', result_table, kinase_names, kinase_names, threshold, download_section, True)

@app.callback(   
  [Output('organism-radioitems', 'value'),
   Output('upload-data', 'contents'),
   Output('upload-data', 'filename'),
   Output('upload-button', 'disabled',allow_duplicate=True)],
   Input('example-button', 'n_clicks'),
   prevent_initial_call=True
)
def click_example(n_clicks) :
  '''
  Click example button
  '''
  ic('click example')
  if n_clicks is None or n_clicks == 0 :
    return no_update
  
  df = pd.read_csv(example_filename)
  csv_string = df.to_csv(index=False)

  #make it look like an uploaded file with utf-8 encoding and base64 encoding
  csv_string = 'data:text/csv;charset=utf-8;base64,' + base64.b64encode(csv_string.encode()).decode()
  return 'yeast', csv_string, example_filename, False


@app.callback([
                Output({'type': 'figure', 'index': MATCH}, 'figure', allow_duplicate=True),
                Output({'type': 'figure', 'index': MATCH}, 'style', allow_duplicate=True),
              ],
                Input({'type': 'figure-accordion', 'index': MATCH}, 'active_item'),
              [
                State({'type': 'figure-accordion', 'index': MATCH}, 'id'),
                State('session-id', 'data'),
                State('session-started', 'data')
              ],
              prevent_initial_call=True
                      
)
def update_figure(active, id, session, started) :
  '''
  Update figure
  '''
  if active is None :
    return no_update
  
  if not started :
    return no_update
  
  #get session from cache
  if session not in cache :
    print('Session not in cache')
    return no_update, no_update
    
  cache_lock.acquire()
  kinaid_session = cache[session]
  cache[session] = kinaid_session
  cache_lock.release()
  
  plot_type = id['index']
  
  if active == 'item-0':
    print('updating: %s ...' % plot_type)
    fig = kinaid_session.get_figure_by_name(plot_type)
    
    style = kinaid_session.get_figure_style_by_name(plot_type)
    return fig, style

  return no_update

#update network graph
@app.callback([
                Output('network-spinner', 'children'),
                Output('timeout-popup', 'displayed', allow_duplicate=True),
              ],
              [
                Input({'type': 'cyto-accordion', 'index': 'network'}, 'active_item'),
                Input('network-match-threshold', 'value')
              ],
              [
                State('session-id', 'data'),
                State('session-started', 'data')
              ],
              prevent_initial_call=True
)
def update_network(active, threshold, session, started) :
  '''
  Draw network figure
  '''
  if active is None :
    return no_update
  #get session from cache
    
  if not started :
    return no_update
  
  if session not in cache :
    print('Session not in cache')
    return no_update, True
    
  cache_lock.acquire()
  kinaid_session = cache[session]
  cache[session] = kinaid_session
  cache_lock.release()

  if active is None :
    return no_update

  
  fig = kinaid_session.get_full_kinase_network_fig(threshold)
  
  return fig, False

@app.callback([
                Output('hub-spinner', 'children'),
                Output('timeout-popup', 'displayed', allow_duplicate=True),
              ],
              [
                Input({'type': 'cyto-accordion', 'index': 'hub'}, 'active_item'),
                Input('hub-kinase-selection', 'value'),
                Input('hub-match-threshold', 'value'),
                Input('hub-only-kinases', 'value')
              ],
              [
                State('session-id', 'data'),
                State('session-started', 'data')
              ],
              prevent_initial_call=True

)
def update_hub_figure(active, kinase_names, threshold, kinase_only, session, started) :
  '''
  Update hub figure
  '''
  if active is None :
    return no_update
  
  if not started :
    return no_update

  if session not in cache :
    print('Session not in cache')
    return no_update, True
    
  cache_lock.acquire()
  kinaid_session = cache[session]
  cache[session] = kinaid_session
  cache_lock.release()

  if len(kinase_only) > 0 :
    kinase_only = (kinase_only[0] == 'only_kinases')
  else:
    kinase_only = False
  
  fig = kinaid_session.get_kinase_hub_fig(kinase_names, threshold, kinase_only)

  return fig, False


@app.callback([
                Output('download-data', 'data'),
                Output('timeout-popup', 'displayed', allow_duplicate=True),
              ],
              Input('download-button', 'n_clicks'),
              [
                State('session-id', 'data'),
                State('session-started', 'data'),
                State('download-radioitems', 'value'),
                State('upload-data', 'filename')
              ],
              prevent_initial_call=True,
)
def download_csv(n_clicks, session, started, dataframe_name, filename) :
  if n_clicks is None or n_clicks == 0 :
    return no_update
  
  if not started :
    return no_update
  
  if session not in cache :
    print('Session not in cache')
    return no_update, True
    
  cache_lock.acquire()
  kinaid_session = cache[session]
  cache[session] = kinaid_session
  cache_lock.release()


  if dataframe_name == 'percentiles' :
    dfs_dict = kinaid_session.get_percentiles_dfs()
    df_list = list(dfs_dict.values())
    fileNameList = [f'{kt}_percentiles.csv'for kt in dfs_dict.keys()]
    out_file = filename.split('.')[0] + '_percentiles.zip'
    return dcc.send_bytes(lambda x: write_archive(x, df_list, fileNameList), out_file), False
  elif dataframe_name == 'network' :
    df = kinaid_session.get_network_df()
    out_file = filename.split('.')[0] + '_network.csv'
    return dcc.send_data_frame(df.to_csv, out_file), False
  elif dataframe_name == 'kinase_network' :
    df = kinaid_session.get_kinase_only_network_df()
    out_file = filename.split('.')[0] + '_kinase_network.csv'
    return dcc.send_data_frame(df.to_csv, out_file), False
  elif dataframe_name == 'kinase_stats' :
    df = kinaid_session.get_stat_df()
    out_file = filename.split('.')[0] + '_kinase_stats.csv'
    return dcc.send_data_frame(df.to_csv, out_file), False

  return no_update, False


@app.callback(
    [Output('kinase-selection', 'value', allow_duplicate=True),
     Output('kinase-selection', 'options', allow_duplicate=True),
     Output('kinase-selection', 'disabled', allow_duplicate=True)],
    Input('ambiguous-check', 'value'),
    [State('organism-radioitems', 'value'),
     State('id-type-dropdown', 'value')],
    prevent_initial_call=True
)
def update_kinase_selection(ambiguous : List[str], organism : str, id_type : str) :
  '''
  Update kinase selection
  '''
  if organism is None:
    return no_update
  symbol_names = list(ortholog_manager.get_orthologs(organism).get_all_kinase_symbols_for_gene_id(id_type, ambiguous=(ambiguous[0] == 'ambiguous' if len(ambiguous) > 0 else False)))
  symbol_options = [{'label': k, 'value': k} for k in symbol_names]
  return symbol_names, symbol_options, False


@app.callback([
                Output({'type': 'figure', 'index': 'peptide_scatter'}, 'figure', allow_duplicate=True),
                Output('timeout-popup', 'displayed', allow_duplicate=True),
              ],
                Input('scatter-kinase-selection', 'value'),
              [
                State('session-id', 'data'),
                State('session-started', 'data')
              ],
              prevent_initial_call=True
)
def update_peptide_scatter_figure(kinase_names, session, started) :
  '''
  Update peptide scatter figure
  '''
  
  print('Updating peptide scatter figure ---kinase')

  if not started :
    return no_update
  
  if session not in cache :
    print('Session not in cache')
    return no_update, True
    
  cache_lock.acquire()
  kinaid_session = cache[session]
  cache[session] = kinaid_session
  cache_lock.release()
  
  fig = kinaid_session.get_peptide_scatter_fig(kinase_names)
  
  return fig, False


@app.callback([
                Output({'type': 'figure', 'index' : 'zscore'}, 'figure', allow_duplicate=True),
                Output('timeout-popup', 'displayed', allow_duplicate=True)
              ],
                Input('fdr-input', 'value'),
              [
                State('session-id', 'data'),
                State('session-started', 'data')
              ],
              prevent_initial_call=True
)
def update_zscore_figure(fdr, session, started) :
    '''
    Update zscore figure
    '''
    if not started :
        return no_update, True
    
    if session not in cache :
        print('Session not in cache')
        return no_update
    
    cache_lock.acquire()
    kinaid_session = cache[session]
    cache[session] = kinaid_session
    cache_lock.release()
    
    fig = kinaid_session.get_zscore_fig(fdr)
    
    return fig, False
  
  
if __name__ == '__main__':
  argparser = argparse.ArgumentParser()
  argparser.add_argument('--port', type=int, default=8050)
  argparser.add_argument('--host', type=str, default='0.0.0.0')
  argparser.add_argument('--no_debug', action='store_false')
  args = argparser.parse_args()
  
  app.title = 'KINAID Dashboard'
  app.layout = serve_layout

  app.run(host=args.host, port=args.port, debug=args.no_debug)