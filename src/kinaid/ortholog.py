import os
from typing import List, Dict, Set
import pandas as pd

class OrthologsWithGeneType :

    
    def __init__(self, df : pd.DataFrame,
                 organism : str,
                 gene_id_type: str,
                 kinase_type: str = 'ST',
                 ambiguous: bool = False,
                 multispecificity_symbols : set = {},
                 longest_long : Dict[str, str] = {},
                 symbol_column: str = 'symbol',
                 kinase_name_column: str = 'kinase_name',
                 gene_id_type_column: str = 'gene_id_type',
                 gene_id_column: str = 'gene_id',
                 kinase_type_column : str = 'kinase_type',
                 short_column : str = 'short',
                 long_column : str = 'long',
                 ambiguous_column : str = 'ambiguous') :
        
        self._organism = organism
        self._kinase_type = kinase_type
        self._gene_id_type = gene_id_type
        self._ambiguous = ambiguous
        self._multispecificity = multispecificity_symbols
        df_copy = df.copy()
        df_copy = df_copy[df_copy[gene_id_type_column] == gene_id_type]
        df_copy = df_copy[df_copy[kinase_type_column] == kinase_type]
        if not ambiguous :
            df_copy = df_copy[~df_copy[ambiguous_column]]
        short_to_long_df = df_copy.groupby(short_column)[long_column].agg(set)
        valid_mapping = short_to_long_df.apply(lambda x: len(x) == 1)
        if not all(valid_mapping) :
            print(short_to_long_df[~valid_mapping])
            raise ValueError('Some short names map to multiple long names')
        
        self._long_to_short_dict = df_copy.set_index(long_column)[short_column].to_dict()
        
        self._covered_kinases = df_copy[kinase_name_column].unique()
        self._covered_symbols = df_copy[symbol_column].unique()
        
        self._long_to_kinase_name_dict = df_copy[[long_column, kinase_name_column]].groupby(long_column).agg('first')[kinase_name_column].to_dict()
        
        self._long_to_gene_id_dict = df_copy[[long_column, gene_id_column]].groupby(long_column).agg(list)[gene_id_column].to_dict()

        self._id_to_long_dict = df_copy.set_index(gene_id_column)[long_column].to_dict()
        
        self._longest_long_dict = longest_long
    
    
    def get_longest_long_dict(self) -> Dict[str, str] :
        return self._longest_long_dict.copy()
    
    def get_long_to_short_dict(self) -> Dict[str, str] :
        return self._long_to_short_dict.copy()
    
    def get_short(self, long_name : str) -> str :
        return self._long_to_short_dict.get(long_name, long_name)
    
    def handle_multispecificity_long(self, long_name : str, use_short = True) -> str :
        return f'({self._kinase_type}){self.get_short(long_name) if use_short else long_name}' if long_name in self._multispecificity else (self.get_short(long_name) if use_short else long_name)
    
    def get_long_names(self) -> List[str] :
        return list(self._long_to_short_dict.keys())
    
    def get_kinase_names(self) -> List[str] :
        return list(self._long_to_kinase_name_dict.values())
    
    def get_long_name_from_id_dict(self) -> Dict[str, str] :
        return self._id_to_long_dict.copy()
    
    def get_long_name_from_id(self, gene_id : str, ignore_warning : bool = False, debug : bool = False) -> str :
        if(gene_id in self._id_to_long_dict) :
            return self._id_to_long_dict[gene_id]
        else :
            if not ignore_warning :
                print(f'Gene ID {gene_id} not found in the orthologs')
                return None
            else :
                return gene_id
    
    def get(self, symbol : str) -> str :
        return self._long_to_kinase_name_dict[symbol]
    
    def print_stats(self) :
        print('\n')
        print(f'Organism: {self._organism}')
        print(f'Kinase type: {self._kinase_type}')
        print(f'Gene ID type: {self._gene_id_type}')
        print(f'Ambiguous: {self._ambiguous}')
        print(f'Number of covered kinases: {len(self._covered_kinases)}')
        print(f'Number of covered symbols: {len(self._covered_symbols)}')
        
class OrganismOrthologs:
    def __init__(self, organism : str, ortholog_file : str, debug = False) :
        self.organism = organism
        self.ortholog_file = ortholog_file
        self.ortholog_df = pd.read_csv(ortholog_file, sep = '\t')
        
        self._gene_id_types = self.ortholog_df['gene_id_type'].unique()
        kinase_types = self.ortholog_df['kinase_type'].unique()
        
        self._kinase_types = set(kinase_types)
        
        multispecificity_dict = {g:self.ortholog_df[self.ortholog_df['gene_id_type'] == g][['kinase_type', 'long']].groupby('long').agg(set)['kinase_type'].to_dict() for g in self._gene_id_types}
        self._multispecificity_longs = {g:{k for k, v in d.items() if len(v) > 1} for g,d in multispecificity_dict.items()}
        
        symbols_to_long = {g:self.ortholog_df[self.ortholog_df['gene_id_type'] == g][['kinase_type', 'symbol', 'long']].groupby('symbol').agg(set)['long'].to_dict() for g in self._gene_id_types}
        multispecificity_symbols = {g:{s:sorted(v, key=len, reverse=True) for s, v in d.items() if len(v) > 1} for g,d in symbols_to_long.items()}
        
        longest_long = {g:{s:v[0] for s,v in d.items()} for g,d in multispecificity_symbols.items()}
        self.orthologs = {}
        for gene_id_type in self._gene_id_types :
            for kinase_type in kinase_types :
                for ambiguous in [True, False] :
                    ortholog = OrthologsWithGeneType(self.ortholog_df,
                                                     organism,
                                                     gene_id_type,
                                                     kinase_type,
                                                     ambiguous,
                                                     self._multispecificity_longs[gene_id_type],
                                                     longest_long[gene_id_type])
                    self.orthologs[(gene_id_type, kinase_type, ambiguous)] = ortholog
                    if debug :
                        ortholog.print_stats()
                    
    
    def get_orthologs(self, gene_id_type: str, kinase_type: str, ambiguous: bool) -> OrthologsWithGeneType :
        #print(gene_id_type, kinase_type, ambiguous)
        return self.orthologs[(gene_id_type, kinase_type, ambiguous)]
        
    def is_supported_kinase_type(self, kinase_type : str) -> bool :
        return kinase_type in self._kinase_types
    
    def get_available_id_types(self) -> List[str] :
        return self._gene_id_types
    
    def get_supported_kinase_types(self) -> Set[str] :
        return self._kinase_types
    
    def get_all_kinase_symbols_for_gene_id(self, gene_id_type : str, ambiguous : bool = True) -> List[str] :
        
        kinase_symbols_long = set([ln for kinase_type in self._kinase_types for ln in self.get_orthologs(gene_id_type, kinase_type, ambiguous).get_long_names()])
        return kinase_symbols_long
      
class OrthologManager:
        
    def __init__(self, orthology_dir : str,
                 suffix : str = '_orthologs_final.tsv',
                 organisms : List[str] = [],
                 human_kinase_file : str = None,
                 debug = False) :
        #get the list of files in the directory with the '_orthologs_final.tsv' extension
        self.ortholog_dir = orthology_dir
        if organisms is None or len(organisms) == 0 :
            ortholog_files = [f for f in os.listdir(orthology_dir) if f.endswith(suffix)]
            self.organism_list = [f.removesuffix(suffix) for f in ortholog_files]
        else :
            ortholog_files = [os.path.join(orthology_dir, organism + suffix) for organism in organisms if organism != 'human']
            if not all([os.path.exists(f) for f in ortholog_files]) :
                raise ValueError('Not all ortholog files exist')
            self.organism_list = list(organisms)
            
            
        self._organism_ortho_dict = {}
        for organism in self.organism_list:
            if debug: 
                print(organism)
            if organism != 'human' :
                self._organism_ortho_dict[organism] = OrganismOrthologs(organism, os.path.join(orthology_dir, organism + suffix), debug=debug)
        
        if human_kinase_file is not None :
            if debug:
                print('Human')
            human_orthologs = OrganismOrthologs('human', human_kinase_file)
            self.organism_list.append('human')
            self._organism_ortho_dict['human'] = human_orthologs
            
    
    def get_orthologs(self, organism : str, gene_id_type : str = None, kinase_type : str = None, ambiguous : bool = True) -> OrthologsWithGeneType :
        if gene_id_type is not None and kinase_type is not None :
            return self._organism_ortho_dict[organism].get_orthologs(gene_id_type, kinase_type, ambiguous)
        else :
            return self._organism_ortho_dict[organism]
if __name__ == '__main__' :
    ortholog_manager = OrthologManager('orthologs')