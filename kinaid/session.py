from typing import Final, Set, Dict
import pandas as pd
from .matching import MatchWithMapping, Scoring
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




class Session :
    ID_COLUMN: Final[str] = '__ENTRY_ID'
    CLEAN_PEPTIDE_COLUMN: Final[str] = '__CLEAN_PEPTIDE'
    

    @staticmethod
    def handle_special_id_column(df : pd.DataFrame, in_id_column : str, id_type : str, out_id_column : str = '__MOD_ID') -> pd.Series :
        if in_id_column not in df.columns :
            return None

        if id_type == 'SGD' :
            df[out_id_column] = df[in_id_column].apply(lambda x: x.split(':')[1] if ':' in x else x)
            df[out_id_column] = 'SGD:' + df[out_id_column]
        else :
            out_id_column = in_id_column
        return out_id_column
        
        
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
                matching_st : MatchWithMapping,
                matching_y : MatchWithMapping,
                selected_kinases : Set[str] = None,
                match_threshold : float = 90.0,
                id_type : str = 'GeneID',
                debug : bool = False
                ) :
        self._session_id = session_id
        self._organism = organism
        self._column_names = column_names
        self._matching_st = matching_st
        self._matching_y = matching_y
        
        available_kinases = set(matching_st._mapping.keys()) | set(matching_y._mapping.keys())
        self._dual_specificity_kinases = set(matching_st._mapping.keys()) & set(matching_y._mapping.keys())
        
        shortened_kinase_names_st = {k : k if '+' not in k else f'({v}-like)' for k,v in matching_st._mapping.items()}
        shortened_kinase_names_y = {k : k if '+' not in k else f'({v}-like)' for k,v in matching_y._mapping.items()}
        
        shortened_kinase_names_st = {k if k not in self._dual_specificity_kinases else f'{k}(ST)':v if k not in self._dual_specificity_kinases else f'{v}(ST)' for k,v in shortened_kinase_names_st.items()}
        shortened_kinase_names_y = {k if k not in self._dual_specificity_kinases else f'{k}(Y)':v if k not in self._dual_specificity_kinases else f'{v}(Y)' for k,v in shortened_kinase_names_y.items()}
        
        self._shortened_kinase_names = shortened_kinase_names_st | shortened_kinase_names_y
        


        if selected_kinases is not None and not selected_kinases.issubset(available_kinases) :
            raise ValueError(f'Selected kinases must be subset of available kinases: {available_kinases}')
        
        if selected_kinases is not None :
            self._selected_kinases = selected_kinases
        else :
            self._selected_kinases = available_kinases
        
        if debug :
            print(f'selected kinases: {self._selected_kinases}')
        try :
            sequence_format = df[self._column_names['peptide']].apply(lambda x : Scoring.get_sequence_format(x))
        except ValueError as e :
            raise ValueError(f'Error in peptide sequence format: {e}')
        
        if len(sequence_format.unique()) != 1 :
            raise ValueError(f'Peptide sequence format must be consistent: {sequence_format.unique()}')
        
        self._mode = sequence_format.unique()[0]
        
        phospho_column_types = df[self._column_names['peptide']].apply(lambda x : Scoring.get_phosphorylation_site(x, self._mode))
        
        #make sure df[Session.PHOSPHO_TYPE_COLUMN] is uppercase
        phospho_column_types = phospho_column_types.str.upper()
        
        relevant_columns = [c for c in column_names.values() if c is not None]
        self._st_df = df[(phospho_column_types == 'S') | (phospho_column_types == 'T')][relevant_columns].copy().reset_index(drop=True)
        self._y_df = df[phospho_column_types == 'Y'][relevant_columns].copy().reset_index(drop=True)
        
        self._st_df[Session.CLEAN_PEPTIDE_COLUMN] = self._st_df[self._column_names['peptide']].apply(lambda x : matching_st._scoring.clean_sequence(x, self._mode))
        self._y_df[Session.CLEAN_PEPTIDE_COLUMN] = self._y_df[self._column_names['peptide']].apply(lambda x : matching_y._scoring.clean_sequence(x, self._mode))
        
        self._st_ids = Session.build_id_column(self._st_df, column_names, 'id')
        self._y_ids = Session.build_id_column(self._y_df, column_names, 'id')
        
        num_st_peptides = len(self._st_df)
        num_y_peptides = len(self._y_df)
                
        if (debug) :
            start_time = perf_counter()
        
        sorted(self._selected_kinases)

        self._st_percentiles = matching_st.get_percentiles_for_selected_kinases(self._st_df[Session.CLEAN_PEPTIDE_COLUMN], self._selected_kinases, self._mode)
        self._y_percentiles = matching_y.get_percentiles_for_selected_kinases(self._y_df[Session.CLEAN_PEPTIDE_COLUMN], self._selected_kinases, self._mode)
        
        if (debug) :
            print(f'Number of ST kinases tested : {len(self._st_percentiles)}')
            print(f'Number of Y kinases tested : {len(self._y_percentiles)}')
        
        self._st_kinase_matches = MatchWithMapping.get_kinase_matches_for_peptides(num_st_peptides, self._st_percentiles, match_threshold)
        self._y_kinase_matches = MatchWithMapping.get_kinase_matches_for_peptides(num_y_peptides, self._y_percentiles, match_threshold)
        
        self._st_peptide_matches = dict(sorted(MatchWithMapping.get_peptide_matches_for_kinases(self._st_percentiles, match_threshold).items()))
        self._y_peptide_matches = dict(sorted(MatchWithMapping.get_peptide_matches_for_kinases(self._y_percentiles, match_threshold).items()))
        
        if (debug) :
            print(f'Elapsed time for percentiles and matches: {perf_counter() - start_time:.2f} seconds')
            
        self._last_accessed = time.time()
    
    def get_percentiles_df(self, kinase_type : str) -> pd.DataFrame :
        if kinase_type == 'ST' :
            df = pd.DataFrame(self._st_percentiles, index = self._st_ids)
            return df
        elif kinase_type == 'Y' :
            df = pd.DataFrame(self._y_percentiles, index = self._y_ids)
            return df
        else :
            raise ValueError(f'Invalid kinase type: {kinase_type}')
    
    def get_kinase_matches_df(self) -> pd.DataFrame :
        st_matches = [','.join(self._st_kinase_matches[i]) for i in range(len(self._st_kinase_matches))]
        y_matches = [','.join(self._y_kinase_matches[i]) for i in range(len(self._y_kinase_matches))]
        
        st_matches = zip(self._st_ids, self._st_df[Session.CLEAN_PEPTIDE_COLUMN], st_matches)
        y_matches = zip(self._y_ids, self._y_df[Session.CLEAN_PEPTIDE_COLUMN], y_matches)
        
        all_matches = list(st_matches) + list(y_matches)
        
        return pd.DataFrame(all_matches, columns=[Session.ID_COLUMN, 'peptide', 'kinase_matches'])
    
    def get_barplot_fig(self) :
        st_counts = [(k,len(pl)) for k,pl in self._st_peptide_matches.items()]
        y_counts = [(k,len(pl)) for k,pl in self._y_peptide_matches.items()]
        
        st_counts = [('ST', k, f'({self._matching_st._mapping[k]}-like)', n) if '+' in k else ('ST', k, k, n) for k,n, in st_counts if n > 0]
        y_counts = [('Y', k, f'({self._matching_y._mapping[k]}-like)', n) if '+' in k else ('Y', k, k, n) for k,n, in y_counts if n > 0]
                
        
        all_counts = st_counts + y_counts
        all_counts = [(t,k,sk,n) if k not in self._dual_specificity_kinases else (t,f'{k}({t})',f'{sk}({t})',n) for t,k,sk,n in all_counts]
        
        df = pd.DataFrame(all_counts, columns=['kinase_type', 'kinase', 'kinase_short', 'count'])

        #sort by kinase type then count
        df = df.sort_values(by=['kinase_type', 'count'], ascending=[False, False])
        
        #map color ST -> light blue, Y -> light orange
        color_map = {'ST' : 'lightblue', 'Y' : 'lightcoral'}
        fig = px.bar(df, x='count', y='kinase_short', color='kinase_type', color_discrete_map=color_map, hover_data={'kinase_type' : False,
                                                                                                                     'kinase' : True,
                                                                                                                     'kinase_short' : False,
                                                                                                                     'count': True})
        fig.update_layout(xaxis_title_text='# peptides', yaxis_title_text='kinase')
        fig.update_layout(xaxis_tickangle=-45)
        fig.update_layout(showlegend=True)
        fig.update_yaxes(showgrid=True)
        
        return fig
    
    @staticmethod
    def kinase_match_string_wrapper(kinase_matches : List[str], width : int = 20) -> str :
        output_string = ''
        line_length = 0
        kinase_matches_exploded = [j for k in kinase_matches for j in k.split('+')]
        sorted(kinase_matches_exploded)
        for k in kinase_matches_exploded :
            if line_length + len(k) > width :
                output_string += '<br>'
                line_length = 0
            output_string += k + ','
            line_length += len(k) + 1
        return output_string[:-1]
    
    def get_peptide_scatter_fig(self, selected_kinases : List[str] = []) -> go.Figure :
        if self._column_names['dependent'] is None :
            raise ValueError('Dependent column must be specified')
        if self._column_names['log2fc'] is None :
            raise ValueError('log2fc column must be specified')
        
        peptide_scatter_st_df = self._st_df[[self._column_names['dependent'], self._column_names['log2fc']]].copy()
        peptide_scatter_st_df['id'] = self._st_ids
        peptide_scatter_st_df['kinase_matches'] = self._st_kinase_matches
        
        peptide_scatter_y_df = self._y_df[[self._column_names['dependent'], self._column_names['log2fc']]].copy()
        peptide_scatter_y_df['id'] = self._y_ids
        peptide_scatter_y_df['kinase_matches'] = self._y_kinase_matches
        
        
        peptide_scatter_df = pd.concat([peptide_scatter_st_df, peptide_scatter_y_df])        
        
        #set_union_lambda = lambda x: set.union(*x)
        #peptide_scatter_df = peptide_scatter_df.groupby('id').agg({self._column_names['dependent'] : 'first',
        #                                                                               self._column_names['log2fc'] : 'first',
        #                                                                               'kinase_matches' : set_union_lambda}).reset_index()
        
        peptide_scatter_df['match'] = peptide_scatter_df['kinase_matches'].apply(lambda x : any(k in selected_kinases for k in x))
        #peptide_scatter_df['kinase'] = peptide_scatter_df['kinase_matches'].apply(lambda x : ', '.join(x))
        peptide_scatter_df['kinase'] = peptide_scatter_df['kinase_matches'].apply(lambda x : Session.kinase_match_string_wrapper(x))

        match_string = Session.kinase_match_string_wrapper(selected_kinases)
        
        peptide_scatter_df['matched'] = peptide_scatter_df['match'].apply(lambda x : match_string if x else 'unmatched')
        peptide_scatter_df['matched'] = peptide_scatter_df['matched'].astype('category')
        
        fig = px.scatter(peptide_scatter_df, x=self._column_names['log2fc'], y=self._column_names['dependent'], color='matched', hover_data={'matched': False,
                                                                                                                                             self._column_names['log2fc']: True,
                                                                                                                                             self._column_names['dependent']: True,
                                                                                                                                             'id':True,
                                                                                                                                             'kinase':True})
        fig.update_layout(xaxis_title_text='log2fc', yaxis_title_text=self._column_names['dependent'])
        
        return fig
    
    def get_heatmap_fig(self, kinase_type : str) -> go.Figure :
        if kinase_type == 'ST' :
            peptide_matches = self._st_percentiles
            ids = self._st_ids
            mapping = self._matching_st._mapping
        elif kinase_type == 'Y' :
            peptide_matches = self._y_percentiles
            ids = self._y_ids
            mapping = self._matching_y._mapping
        else :
            raise ValueError(f'Invalid kinase type: {kinase_type}')
        
        heatmap_df = pd.DataFrame.from_dict(peptide_matches, orient='index')
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
        hover_df.index = heatmap_df.index

        heatmap_df.index = [f'({mapping[k]}-like)' if '-like)' in k else k for k in heatmap_df.index]
        hover_df.index = heatmap_df.index

        fig = go.Figure(data=go.Heatmap(
            z=heatmap_df.values,
            x=heatmap_df.columns,
            y=heatmap_df.index,
            colorscale='Viridis',
            hoverinfo='text',
            text=hover_df.values)
        )
        
        fig.update_xaxes(showticklabels=False)

        fig.update_traces(colorbar_orientation='h')
        
        return fig
    
    def get_stat_df(self) -> pd.DataFrame :
        population_log2fc = pd.concat([self._st_df[self._column_names['log2fc']],self._y_df[self._column_names['log2fc']]], ignore_index=True)
        population_log2fc = population_log2fc[population_log2fc.notnull()]
        population_log2fc_mean = np.mean(population_log2fc)
        population_log2fc_std = np.std(population_log2fc)
        
        log2fc_st_dict = {k : self._st_df.loc[self._st_peptide_matches[k], self._column_names['log2fc']] for k in self._st_peptide_matches.keys()}
        log2fc_y_dict = {k : self._y_df.loc[self._y_peptide_matches[k], self._column_names['log2fc']] for k in self._y_peptide_matches.keys()}
        
        log2fc_st_dict = {k if k not in self._dual_specificity_kinases else f'{k}(ST)':l  for k,l in log2fc_st_dict.items()}
        log2fc_y_dict = {k if k not in self._dual_specificity_kinases else f'{k}(Y)':l  for k,l in log2fc_y_dict.items()}
        
        log2fc_dict = log2fc_st_dict | log2fc_y_dict
        
        log2fc_tuples = [(k, len(l), np.mean(l), np.std(l)) for k,l in log2fc_dict.items() if len(l) > 0]
        log2fc_df = pd.DataFrame(log2fc_tuples, columns=['kinase', 'n', 'mean', 'std'])
        log2fc_df['zscore'] = (log2fc_df['mean'] - population_log2fc_mean) / (population_log2fc_std / np.sqrt(log2fc_df['n']))
        log2fc_df['p'] = log2fc_df['zscore'].apply(lambda x: 1 - norm.cdf(x) if x > 0 else norm.cdf(x))
        log2fc_df['p_adj'] = fdrcorrection(log2fc_df['p'])[1]
        
        return log2fc_df
    
    def get_zscore_fig(self, fdr_threshold : float = 0.05) -> go.Figure :
        zscore_df = self.get_stat_df().copy()
        zscore_df['zscore_sig'] = zscore_df['p_adj'].apply(lambda x: x <= fdr_threshold)
        
        zscore_df.sort_values(by='zscore', ascending=False, inplace=True)
        zscore_df['kinase_short'] = zscore_df['kinase'].map(self._shortened_kinase_names)
        fig = px.bar(zscore_df, x='zscore', y='kinase_short', hover_data={'zscore_sig': False, 'kinase_short' : False, 'zscore' : True, 'p_adj' : True, 'kinase' : True}, color='zscore_sig')
        
        fig.update_layout(yaxis_categoryorder='total ascending')
        fig.update_layout(legend_title_text='FDR <= %0.2f' % fdr_threshold)
        fig.update_yaxes(showgrid=True)
        
        return fig
    
    def get_kinase_scatter_fig(self) :
        kinase_scatter_df = self.get_stat_df().copy()
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





        
        
        
        
        
        


        
        