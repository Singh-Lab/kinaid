from typing import Final, Set, Dict
import pandas as pd
from .matching import MatchWithMapping, Scoring
import time
from time import perf_counter
import plotly.express as px



class Session :
    ID_COLUMN: Final[str] = '__ENTRY_ID'
    CLEAN_PEPTIDE_COLUMN: Final[str] = '__CLEAN_PEPTIDE'
    
    @staticmethod
    def build_id_column(df : pd.DataFrame, column_names : Dict[str, str]) -> pd.DataFrame :
        id_column = None
        if column_names['id'] is None and column_names['site'] is None :
            id_column = df.index
        elif column_names['id'] is not None and column_names['site'] is None :
            id_column = df[column_names['id']].astype(str) + '_' + df.index.astype(str)
        elif column_names['id'] is not None and column_names['site'] is not None :
            id_column = df[column_names['id']].astype(str) + '_' + df[column_names['site']].astype(str)
        return id_column
            
    def __init__(self,
                session_id : str,
                organism : str,
                df : pd.DataFrame,
                column_names : Dict[str, str],
                matching_st : MatchWithMapping,
                matching_y : MatchWithMapping,
                selected_kinases : Set[str] = None,
                match_threshold : float = 90.0,
                debug : bool = False
                ) :
        self._session_id = session_id
        self._organism = organism
        self._column_names = column_names
        self._matching_st = matching_st
        self._matching_y = matching_y
        
        available_kinases = set(matching_st._mapping.keys()) | set(matching_y._mapping.keys())
        self._dual_specificity_kinases = set(matching_st._mapping.keys()) & set(matching_y._mapping.keys())
        
        
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
        
        self._st_ids = Session.build_id_column(self._st_df, column_names)
        self._y_ids = Session.build_id_column(self._y_df, column_names)
        
        num_st_peptides = len(self._st_df)
        num_y_peptides = len(self._y_df)
        
        print(num_y_peptides)
        
        if (debug) :
            start_time = perf_counter()

        self._st_percentiles = matching_st.get_percentiles_for_selected_kinases(self._st_df[Session.CLEAN_PEPTIDE_COLUMN], self._selected_kinases, self._mode)
        self._y_percentiles = matching_y.get_percentiles_for_selected_kinases(self._y_df[Session.CLEAN_PEPTIDE_COLUMN], self._selected_kinases, self._mode)
        
        if (debug) :
            print(f'Number of ST kinases tested : {len(self._st_percentiles)}')
            print(f'Number of Y kinases tested : {len(self._y_percentiles)}')
        
        self._st_kinase_matches = sorted(MatchWithMapping.get_kinase_matches_for_peptides(num_st_peptides, self._st_percentiles, match_threshold))
        self._y_kinase_matches = sorted(MatchWithMapping.get_kinase_matches_for_peptides(num_y_peptides, self._y_percentiles, match_threshold))
        
        self._st_peptide_matches = MatchWithMapping.get_peptide_matches_for_kinases(self._st_percentiles, match_threshold)
        self._y_peptide_matches = MatchWithMapping.get_peptide_matches_for_kinases(self._y_percentiles, match_threshold)
        
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
        
        st_counts = [('ST',k,n) if k not in self._dual_specificity_kinases else ('ST',k+'(ST)',n) for k,n in st_counts if n > 0]
        y_counts = [('Y',k,n) if k not in self._dual_specificity_kinases else ('Y',k+'(Y)',n) for k,n in y_counts if n > 0]
        
        
        all_counts = st_counts + y_counts
        
        df = pd.DataFrame(all_counts, columns=['kinase_type', 'kinase', 'count'])
        print(df.head())

        #sort by kinase type then count
        df = df.sort_values(by=['kinase_type', 'count'], ascending=[False, False])
        
        #map color ST -> light blue, Y -> light orange
        color_map = {'ST' : 'lightblue', 'Y' : 'lightcoral'}
        fig = px.bar(df, x='count', y='kinase', color='kinase_type', color_discrete_map=color_map, hover_data=['count'])
        fig.update_layout(xaxis_title_text='# peptides', yaxis_title_text='kinase')
        fig.update_layout(xaxis_tickangle=-45)
        fig.update_layout(showlegend=True)
        fig.update_yaxes(showgrid=True)
        
        return fig
        
        
        
        