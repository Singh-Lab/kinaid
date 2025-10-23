from typing import Final, Set, Dict
import pandas as pd
from .matching import MatchWithMapping, Scoring, PeptideBackground
from .ortholog import OrthologManager
import time
from time import perf_counter
import plotly.express as px
from typing import List
from scipy.cluster.hierarchy import linkage, dendrogram
import plotly.graph_objects as go
import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import fdrcorrection
import math
import dash_cytoscape as cyto
from plotly.subplots import make_subplots



class Session :
    ID_COLUMN: Final[str] = '__ENTRY_ID'
    CLEAN_PEPTIDE_COLUMN: Final[str] = '__CLEAN_PEPTIDE'
    cyto.load_extra_layouts()

    

    @staticmethod
    def handle_special_id_column(id_column : pd.Series, id_type : str, out_id_column : str = '__MOD_ID') -> pd.Series :
        id_column_new = id_column.copy()
        if id_type == 'SGD' :
            id_column_new = id_column.apply(lambda x: x.split(':')[1] if ':' in x else x)
            id_column_new = 'SGD:' + id_column_new
        return id_column_new
        
        
    @staticmethod
    def build_id_column(df : pd.DataFrame, column_names : Dict[str, str], id_column = 'id') -> pd.Series :
        ids_column = None
        if column_names[id_column] is None and column_names['site'] is None :
            ids_column = df.index
        elif column_names[id_column] is not None and column_names['site'] is None :
            ids_column = df[column_names[id_column]].astype(str) + '_' + df.index.astype(str)
        elif column_names[id_column] is not None and column_names['site'] is not None :
            ids_column = df[column_names[id_column]].astype(str) + '_' + df[column_names['site']].astype(str)
        return ids_column
    
    
    def __init__(self,
                 session_id : str,
                 organism : str,
                 df : pd.DataFrame,
                 column_names : Dict[str, str],
                 scoring : Dict[str , Scoring],
                 background : Dict[str , PeptideBackground],
                 ortholog_manager : OrthologManager,
                 selected_symbols : Set[str] = set(),
                 match_threshold : float = 90.0,
                 id_type : str = 'GeneID',
                 ambiguous : bool = True,
                 debug = False) :
        
        self._session_id = session_id
        self._organism = organism
        self._column_names = column_names
        self._network_df = None
                
        try :
            sequence_format = df.apply(lambda x: Scoring.get_sequence_format(x[column_names['peptide']], x.name), axis=1)
        except ValueError as e :
            raise ValueError(f'Error in peptide sequence format: {e}')
        
        if len(sequence_format.unique()) != 1 :
            raise ValueError(f'Peptide sequence format must be consistent: {sequence_format.unique()}')
            
        self._mode = sequence_format.unique()[0]

        phospho_column_types = df[self._column_names['peptide']].apply(lambda x: Scoring.get_phosphorylation_site_type(x, self._mode)) 
        
        self._orthologs = {}
        matching = {}
        relevant_columns = [c for c in column_names.values() if c is not None]
        self._dfs = {}
        organism_orthologs = ortholog_manager.get_orthologs(organism=organism)
        self._num_peptides = {}
        self._all_selected_symbols = selected_symbols.copy()
        self._supported_kinase_types = list()
        sorted(self._all_selected_symbols)
        self._selected_symbols = {}

        supported_symbols = set()
        for phospho_type in sorted(phospho_column_types.unique()) :
            if organism_orthologs.is_supported_kinase_type(phospho_type) and phospho_type in scoring and phospho_type in background :
                self._supported_kinase_types.append(phospho_type)
                self._orthologs[phospho_type] = organism_orthologs.get_orthologs(gene_id_type=id_type, kinase_type=phospho_type, ambiguous=ambiguous)
                matching[phospho_type] = MatchWithMapping(scoring=scoring[phospho_type], background=background[phospho_type], mapping=self._orthologs[phospho_type], selected_symbols=self._all_selected_symbols)
                self._selected_symbols[phospho_type] = matching[phospho_type].get_selected_symbols()
                supported_symbols |= self._selected_symbols[phospho_type]
                self._dfs[phospho_type] = df[phospho_column_types == phospho_type][relevant_columns].copy().reset_index(drop=True)
                self._dfs[phospho_type][Session.CLEAN_PEPTIDE_COLUMN] = self._dfs[phospho_type][self._column_names['peptide']].apply(lambda x: scoring[phospho_type].clean_sequence(x, self._mode))
                self._dfs[phospho_type][column_names['id']] = Session.handle_special_id_column(self._dfs[phospho_type][column_names['id']], id_type)
                self._dfs[phospho_type][Session.ID_COLUMN] = Session.build_id_column(self._dfs[phospho_type], column_names, 'id')
                duplicated_ids = self._dfs[phospho_type][Session.ID_COLUMN].duplicated()
                if duplicated_ids.any() :
                    print(f'Warning: Duplicated IDs in {phospho_type} data: {self._dfs[phospho_type][Session.ID_COLUMN][duplicated_ids]}')
                    #drop duplicates
                    self._dfs[phospho_type] = self._dfs[phospho_type][~duplicated_ids]
                self._dfs[phospho_type].set_index(Session.ID_COLUMN, inplace=True)
                self._num_peptides[phospho_type] = len(self._dfs[phospho_type])
        
        self._supported_kinase_types = sorted(self._supported_kinase_types)
        self._supported_kinase_types = set(self._supported_kinase_types)
        self._all_selected_symbols = supported_symbols
        
        if (debug) :
            start_time = perf_counter()
        
        self._percentiles = {}
        self._kinase_matches = {}
        self._peptide_matches = {}
        
        for phospho_type in self._supported_kinase_types :
            self._percentiles[phospho_type] = matching[phospho_type].get_percentiles_for_selected_kinases(self._dfs[phospho_type][Session.CLEAN_PEPTIDE_COLUMN], self._mode)
            self._kinase_matches[phospho_type] = MatchWithMapping.get_kinase_matches_for_peptides(self._num_peptides[phospho_type], self._percentiles[phospho_type], match_threshold)
            self._peptide_matches[phospho_type] = dict(sorted(MatchWithMapping.get_peptide_matches_for_kinases(self._percentiles[phospho_type], match_threshold).items()))
        
        if (debug) :
            print(f'Elapsed time for percentiles and matches: {perf_counter() - start_time:.2f} seconds')
        
        self._last_accessed = time.time()
    
    def get_percentiles_df(self, kinase_type : str) -> pd.DataFrame :
        if kinase_type in self._supported_kinase_types :
            df = pd.DataFrame(self._percentiles[kinase_type], index = self._dfs[kinase_type].index)
            return df
        else :
            raise ValueError(f'Invalid kinase type: {kinase_type}')
        
    def get_percentiles_dfs(self) -> Dict[str, pd.DataFrame] :
        return {kt : self.get_percentiles_df(kt) for kt in self._supported_kinase_types}
    
    def get_kinase_matches_df(self) -> pd.DataFrame :
        matches = {}
        for kinase_type in self._supported_kinase_types :
            matches[kinase_type] = [','.join(self._kinase_matches[kinase_type][i]) for i in range(len(self._kinase_matches[kinase_type]))]
            matches[kinase_type] = zip(self._dfs[kinase_type].index, self._dfs[kinase_type][Session.CLEAN_PEPTIDE_COLUMN], matches[kinase_type])
        
        all_matches = []
        for kinase_type in self._supported_kinase_types :
            all_matches += list(matches[kinase_type])
        
        return pd.DataFrame(all_matches, columns=[Session.ID_COLUMN, 'peptide', 'kinase_matches'])
    
    def get_counts_barplot_fig(self) :
        counts = [(kt,k,self._orthologs[kt].handle_multispecificity_long(k),len(pl)) for kt in self._supported_kinase_types for k,pl in self._peptide_matches[kt].items()]
        df = pd.DataFrame(counts, columns=['kinase_type', 'kinase', 'kinase_short', 'count'])
        df = df.sort_values(by=['kinase_type', 'count'], ascending=[False, True])
        #color_map = {'ST' : 'lightblue', 'Y' : 'lightcoral'}
        fig = px.bar(df, x='count', y='kinase_short', color='kinase_type', hover_data={'kinase_type' : False,
                                                                                       'kinase' : True,
                                                                                       'kinase_short' : False,
                                                                                       'count': True})
        fig.update_layout(xaxis_title_text='# peptides', yaxis_title_text='kinase')
        fig.update_layout(xaxis_tickangle=-45)
        fig.update_layout(showlegend=True)
        fig.update_yaxes(showgrid=True)
        
        return fig

    @staticmethod
    def __wrap__text(kinase_matches : Set[str], width : int = 20, sep=',') -> str :
        output_string = ''
        line_length = 0
        kinase_matches_exploded = [j for k in kinase_matches for j in k.split(sep)]
        for k in kinase_matches_exploded :
            if line_length + len(k) > width :
                output_string += '<br>'
                line_length = 0
            output_string += k + ','
            line_length += len(k) + 1
        return output_string[:-1]
    
    def get_peptide_scatter_fig(self, selected_kinases : Set[str] = set()) -> go.Figure :
        if self._column_names['dependent'] is None :
            raise ValueError('Dependent column must be specified')
        if self._column_names['log2fc'] is None :
            raise ValueError('log2fc column must be specified')
        
        peptide_scatter_dfs = []
        for kinase_type in self._supported_kinase_types :
            peptide_scatter_df = self._dfs[kinase_type][[self._column_names['dependent'], self._column_names['log2fc']]].copy()
            peptide_scatter_df['id'] = peptide_scatter_df.index
            peptide_scatter_df['kinase_matches'] = self._kinase_matches[kinase_type]
            peptide_scatter_dfs.append(peptide_scatter_df)
        
        peptide_scatter_df = pd.concat(peptide_scatter_dfs)
        peptide_scatter_df['kinase'] = peptide_scatter_df['kinase_matches'].apply(lambda x : Session.__wrap__text(x))
        peptide_scatter_df['match'] = peptide_scatter_df['kinase_matches'].apply(lambda x : any(k in selected_kinases for k in x))
        
        selected_kinases = set(sorted(list(selected_kinases)))
        match_string = Session.__wrap__text(selected_kinases)
        peptide_scatter_df['matched'] = peptide_scatter_df['match'].apply(lambda x : match_string if x else 'unmatched')
        peptide_scatter_df['matched'] = peptide_scatter_df['matched'].astype('category')

        fig = px.scatter(peptide_scatter_df, x=self._column_names['log2fc'], y=self._column_names['dependent'], color='matched', hover_data={'matched': False,
                                                                                                                                             self._column_names['log2fc']: True,
                                                                                                                                             self._column_names['dependent']: True,
                                                                                                                                             'id':True,
                                                                                                                                             'kinase':True})
        fig.update_layout(xaxis_title_text='log2fc', yaxis_title_text=self._column_names['dependent'])
        
        return fig
    
    def get_heatmap_fig(self, kinase_type : str, coloraxis : str ='coloraxis') -> go.Figure :
        if kinase_type not in self._supported_kinase_types :
            raise ValueError(f'Invalid kinase type: {kinase_type}')
        percentiles = self._percentiles[kinase_type]
        ids = self._dfs[kinase_type].index
        
        heatmap_df = pd.DataFrame.from_dict(percentiles, orient='index')
        heatmap_df.columns = ids
        
        Z = linkage(heatmap_df, method='ward', metric='euclidean')
        row_order = dendrogram(Z, no_plot=True)['leaves']
        
        heatmap_df = heatmap_df.iloc[row_order]
        
        heatmap_t_df = heatmap_df.transpose()
        Z = linkage(heatmap_t_df, method='ward', metric='euclidean')
        col_order = dendrogram(Z, no_plot=True)['leaves']
        
        heatmap_df = heatmap_df.iloc[:, col_order]

        hover_data = [[f"kin: {kin}<br>pep: {pep}<br>per: {heatmap_df.loc[kin, pep]:.2f}" for pep in heatmap_df.columns] for kin in heatmap_df.index]
        hover_df = pd.DataFrame(hover_data)
        hover_df.columns = heatmap_df.columns
        
        heatmap_df.index = [self._orthologs[kinase_type].handle_multispecificity_long(k) for k in heatmap_df.index]
        hover_df.index = heatmap_df.index
        
        fig = go.Figure(data=go.Heatmap(
            z=heatmap_df.values,
            x=heatmap_df.columns,
            y=heatmap_df.index,
            coloraxis=coloraxis,
            colorscale='Viridis',
            hoverinfo='text',
            text=hover_df.values)
        )
        
        fig.update_xaxes(showticklabels=False)

        #fig.update_traces(colorbar_orientation='h')
        #fig.update_layout(coloraxis = {'colorscale':'viridis', colorbar=dict(title='percentile')},)
        fig.update_layout(coloraxis = dict(colorscale='viridis',
                                           colorbar=dict(title='percentile', orientation='h'),
                                           cmin=0,
                                           cmax=100))
        
        return fig

    def get_all_heatmap_figs(self) -> go.Figure :
        kinase_types = sorted(self._supported_kinase_types)
        figure_heights = [len(self._selected_symbols[kt]) for kt in kinase_types]
        total_height = sum(figure_heights)
        figure_heights = [h / total_height for h in figure_heights]
                
        fig = make_subplots(rows=len(kinase_types), cols=1, shared_xaxes=False, vertical_spacing=0.02, row_heights=figure_heights)
        for row,kt in enumerate(kinase_types) :
            hfig = self.get_heatmap_fig(kt)
            trace = hfig.data[0]
            fig.add_trace(trace, row=row+1, col=1)
            fig.update_xaxes(showticklabels=False, row=row+1, col=1)
        
        fig.update_layout(coloraxis = dict(colorscale='viridis',
                                           colorbar=dict(title='percentile',
                                                         orientation='h'),
                                           cmin=0,
                                           cmax=100),
                          margin=dict(t=100, b=20, l=20, r=20))
        #fig.update_traces(colorbar_orientation='h')

        return fig
    
    def get_stat_df(self, combine_populations : bool = True) -> pd.DataFrame :
        log2fcs = {kt:self._dfs[kt][self._column_names['log2fc']] for kt in self._supported_kinase_types}
        log2fcs = {kt:l[l.notnull()] for kt,l in log2fcs.items()}
        means = {kt:np.mean(l) for kt,l in log2fcs.items()}
        stds = {kt:np.std(l) for kt,l in log2fcs.items()}
        total_n = {kt:len(l) for kt,l in log2fcs.items()}
        
        log2fc_dict = {kt : {k: self._dfs[kt].iloc[self._peptide_matches[kt][k]][self._column_names['log2fc']] for k in self._peptide_matches[kt].keys()} for kt in self._supported_kinase_types}
                
        log2fc_tuples = {kt : [(self._orthologs[kt].handle_multispecificity_long(k, use_short=False), self._orthologs[kt].handle_multispecificity_long(k), len(l), np.mean(l), np.std(l)) for k,l in ld.items() if len(l) > 0] for kt,ld in log2fc_dict.items()}

        log2fc_dfs = {kt : pd.DataFrame(log2fc_tuples[kt], columns=['kinase', 'short', 'n', 'mean', 'std']) for kt in self._supported_kinase_types}
        
        if combine_populations :
            combined_total_n = sum(total_n.values())
            mean = sum([m * total_n[kt] for kt,m in means.items()]) / combined_total_n
            std = np.sqrt(sum([total_n[kt] * (stds[kt] ** 2 + (means[kt] - mean) ** 2) for kt in self._supported_kinase_types]) / combined_total_n)
            means = {kt:mean for kt in self._supported_kinase_types}
            stds = {kt:std for kt in self._supported_kinase_types}
            total_n = {kt:combined_total_n for kt in self._supported_kinase_types}
        
        for kt in self._supported_kinase_types :
            log2fc_dfs[kt]['zscore'] = (log2fc_dfs[kt]['mean'] - means[kt]) / (stds[kt] / np.sqrt(log2fc_dfs[kt]['n']))
            log2fc_dfs[kt]['p'] = log2fc_dfs[kt]['zscore'].apply(lambda x: 1 - norm.cdf(x) if x > 0 else norm.cdf(x))
            log2fc_dfs[kt]['p_adj'] = fdrcorrection(log2fc_dfs[kt]['p'])[1]
        
        log2fc_df = pd.concat([df for df in log2fc_dfs.values()], ignore_index=True)

        if combine_populations :
            log2fc_df['p_adj'] = fdrcorrection(log2fc_df['p'])[1]
        
        return log2fc_df
        
    def get_zscore_fig(self, fdr_threshold : float = 0.05, combine_populations : bool = True) -> go.Figure :
        zscore_df = self.get_stat_df(combine_populations).copy()
        zscore_df['zscore_sig'] = zscore_df['p_adj'].apply(lambda x: x <= fdr_threshold)
        
        zscore_df.sort_values(by='zscore', ascending=False, inplace=True)
        fig = px.bar(zscore_df, x='zscore', y='short', hover_data={'zscore_sig': False, 'short' : False, 'mean': True, 'zscore' : True, 'p_adj' : True, 'kinase' : True}, color='zscore_sig')
        
        fig.update_layout(yaxis_categoryorder='total ascending')
        fig.update_layout(legend_title_text='FDR <= %0.2f' % fdr_threshold)
        fig.update_yaxes(showgrid=True)
        
        return fig
    
    def get_kinase_scatter_fig(self, combine_populations : bool = True) :
        kinase_scatter_df = self.get_stat_df(combine_populations).copy()
        kinase_scatter_df.sort_values(by='mean', ascending=False, inplace=True)
        smallest_non_zero = kinase_scatter_df[kinase_scatter_df['p_adj'] > 0]['p_adj'].min()
        neglog_smallest_non_zero = -math.log10(smallest_non_zero)
        kinase_scatter_df['-log10(p_adj)'] = kinase_scatter_df['p_adj'].apply(lambda x: -math.log10(x) if x > 0 else neglog_smallest_non_zero)
        
        fig = px.scatter(kinase_scatter_df, x='mean', y='-log10(p_adj)', hover_data={'kinase' : True, 'mean' : True,  'p_adj' : True})
        fig.update_layout(yaxis_title_text='-log10(adj p-value)')
        fig.update_layout(xaxis_title_text='mean log2FC')
        fig.update_layout(showlegend=False)
        fig.update_yaxes(showgrid=True)
        
        return fig
    def convert_network_to_cytoscape(self,
                                     network_df : pd.DataFrame,
                                     node1_col : str = 'substrate_id',
                                     node2_col : str = 'kinase',
                                     kinase_to_kinase_col : str = 'kinase_to_kinase',
                                     highlighted_symbols : Set[str] = set()) -> cyto.Cytoscape :

        network_df = network_df.copy()
        kinase_nodes = set(network_df[node2_col]) | set(network_df[network_df[kinase_to_kinase_col]][node1_col])
        nodes = set(network_df[node1_col]) | kinase_nodes
        stylesheet = [
            {
                'selector': 'node',
                'style': {
                    'label': 'data(id)'
                }
            },
            {
                'selector' : '.kinase',
                'style' : {'background-color' : 'blue'}
            },
            {
                'selector' : '.highlighted',
                'style' : {'background-color' : 'orange'}
            },
            {
                'selector' : '.substrate',
                'style' : {'background-color' : 'grey'}
            },
            {
                'selector' : 'edge',
                'style' : {'target-arrow-shape' : 'triangle', 'target-arrow-color' : 'green', 'target-arrow-fill' : 'filled','curve-style': 'bezier'}
            }

        ]
        
        nodes = [{'data' : {'id' : n, 'label' : n}, 'classes' : 'highlighted' if n in highlighted_symbols else 'kinase' if n in kinase_nodes else 'substrate'} for n in nodes]
        edges = [{'data' : {'source' : s, 'target' : t}} for s,t in zip(network_df['kinase'], network_df['substrate_id'])]
        
        elements = nodes + edges
        cyto_plot = cyto.Cytoscape(
            id='kinase_network',
            layout={'name': 'klay', 'animate': True},
            style={'width': '100%', 'height': '1000px'},
            stylesheet=stylesheet,
            elements=elements
        )

        return cyto_plot
    
    def get_network_df(self) -> pd.DataFrame :
        if self._network_df is not None :
            return self._network_df.copy()
        
        network_dfs = []
        for kinase_type in self._supported_kinase_types :
            percentiles_df = pd.DataFrame.from_dict(self._percentiles[kinase_type], orient='columns')
            if len(percentiles_df) != len(self._dfs[kinase_type]) :
                raise ValueError('Percentiles and source dataframes must have the same length')
            percentiles_df['substrate_id'] = self._dfs[kinase_type][self._column_names['id']].reset_index(drop=True)
            network_df = pd.melt(percentiles_df, id_vars=['substrate_id'], value_vars=self._percentiles[kinase_type].keys(), var_name='kinase', value_name='percentile')
            network_dfs.append(network_df)
        network_df = pd.concat(network_dfs)
        
        network_df['kinase_to_kinase'] = False
        for kinase_type in self._supported_kinase_types :
            long_names = set(self._orthologs[kinase_type].get_long_names())
            long_name_from_id_dict = self._orthologs[kinase_type].get_long_name_from_id_dict()
            network_df['substrate_id'] = network_df['substrate_id'].apply(lambda x: long_name_from_id_dict.get(x, x))
            network_df['kinase_to_kinase'] = network_df.apply(lambda x: (x['substrate_id'] in long_names) or x['kinase_to_kinase'], axis=1)
            
        network_df.sort_values(by='percentile', ascending=False, inplace=True)
        network_df.drop_duplicates(subset=['substrate_id', 'kinase'], inplace=True)
        
        network_df['kinase'] = network_df['kinase'].map(lambda x: self._orthologs[kinase_type].get_longest_long_dict().get(x, x))
        network_df['substrate_id'] = network_df['substrate_id'].map(lambda x: self._orthologs[kinase_type].get_longest_long_dict().get(x, x))
        
        self._network_df = network_df.copy()
        
        return network_df
    
    def get_kinase_only_network_df(self) -> pd.DataFrame :
        network_df = self.get_network_df()
        network_df = network_df[network_df['kinase_to_kinase']]
        return network_df
    
    def get_kinase_hub_network_df(self, selected_symbols_longest : Set[str] = set()) -> pd.DataFrame :
        network_df = self.get_network_df()
        network_df = network_df[network_df['kinase'].isin(selected_symbols_longest) | network_df['substrate_id'].isin(selected_symbols_longest)]

        return network_df
    
    def get_kinase_hub_fig(self, selected_symbols : Set[str] = set(), threshold : float = 90.0, kinase_only : bool = True) -> cyto.Cytoscape :
        selected_symbols_longest = {self._orthologs[kt].get_longest_long_dict().get(k, k) for kt in self._supported_kinase_types for k in selected_symbols}
        hub_network_df = self.get_kinase_hub_network_df(selected_symbols_longest)
        hub_network_df = hub_network_df[hub_network_df['percentile'] >= threshold]
        if kinase_only :
            hub_network_df = hub_network_df[hub_network_df['kinase_to_kinase']]
        
        highlighted_symbols_short = {self._orthologs[kt].get_short(k) for kt in self._supported_kinase_types for k in selected_symbols}
        
        for kinase_type in self._supported_kinase_types :
            hub_network_df['kinase'] = hub_network_df['kinase'].map(lambda x: self._orthologs[kinase_type].get_short(x))
            hub_network_df['substrate_id'] = hub_network_df['substrate_id'].map(lambda x: self._orthologs[kinase_type].get_short(x))
        
        return self.convert_network_to_cytoscape(hub_network_df, highlighted_symbols=highlighted_symbols_short)
    
    def get_size_of_kinase_network(self, threshold : float = 99.0) -> int :
        full_network_df = self.get_kinase_only_network_df()
        full_network_df = full_network_df[full_network_df['percentile'] >= threshold]
        return len(full_network_df)
    
    def get_full_kinase_network_fig(self, threshold : float = 99.0) -> cyto.Cytoscape :
        full_network_df = self.get_kinase_only_network_df()
        full_network_df = full_network_df[full_network_df['percentile'] >= threshold]

        for kinase_type in self._supported_kinase_types :
            full_network_df['kinase'] = full_network_df['kinase'].map(lambda x: self._orthologs[kinase_type].get_short(x))
            full_network_df['substrate_id'] = full_network_df['substrate_id'].map(lambda x: self._orthologs[kinase_type].get_short(x))
        return self.convert_network_to_cytoscape(full_network_df)

    def get_figure_by_name(self, name : str):
        self.last_accessed = time.time()
        
        if name == 'barplot':
            return self.get_counts_barplot_fig()
        elif name == 'peptide_scatter' or name == 'peptide_volcano':
            return self.get_peptide_scatter_fig()
        elif name == 'heatmap':
            return self.get_all_heatmap_figs()
        elif name == 'zscore':
            return self.get_zscore_fig()
        elif name == 'kinase_scatter':
            return self.get_kinase_scatter_fig()
        elif name == 'network':
            return self.get_full_kinase_network_fig()
        elif name == 'hub':
            return self.get_kinase_hub_fig()
        else:
            return None
    
    def get_figure_style_by_name(self, name : str) :
        self.last_accessed = time.time()
        #get number of kinases from self.kinase_peptides
        num_kinases = len(self._all_selected_symbols)
        plot_height = 1.2 * num_kinases
        plot_style = {'height': '%drem' % plot_height}

        if name in ['zscore', 'barplot', 'heatmap'] :
            return plot_style
        else :
            return {'height': '50rem'}
    