import os
from typing import List, Dict
import pandas as pd

class OrthologsWithGeneType :
    def __init__(self, df : pd.DataFrame,
                 organism : str,
                 gene_id_type: str,
                 kinase_type: str = 'ST',
                 ambiguous: bool = False,
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
        
        self._long_to_kinase_name_dict = df_copy[[long_column, kinase_name_column]].groupby(long_column).agg(list)[kinase_name_column].to_dict()
        
        self._long_to_gene_id_dict = df_copy[[long_column, gene_id_column]].groupby(long_column).agg(list)[gene_id_column].to_dict()
        
    def print_stats(self) :
        print('\n')
        print(f'Organism: {self._organism}')
        print(f'Kinase type: {self._kinase_type}')
        print(f'Gene ID type: {self._gene_id_type}')
        print(f'Ambiguous: {self._ambiguous}')
        print(f'Number of covered kinases: {len(self._covered_kinases)}')
        print(f'Number of covered symbols: {len(self._covered_symbols)}')
        
class SpeciesOrthologs:
    def __init__(self, species : str, ortholog_file : str) :
        self.species = species
        self.ortholog_file = ortholog_file
        self.ortholog_df = pd.read_csv(ortholog_file, sep = '\t')
        
        gene_id_types = self.ortholog_df['gene_id_type'].unique()
        kinase_types = self.ortholog_df['kinase_type'].unique()
        
        self.orthologs = {}
        for gene_id_type in gene_id_types :
            for kinase_type in kinase_types :
                for ambiguous in [True, False] :
                    ortholog = OrthologsWithGeneType(self.ortholog_df,
                                                     species,
                                                     gene_id_type,
                                                     kinase_type,
                                                     ambiguous)
                    self.orthologs[(gene_id_type, kinase_type, ambiguous)] = ortholog
                    ortholog.print_stats()
    
    def get_orthologs(self, gene_id_type: str, kinase_type: str, ambiguous: bool) -> OrthologsWithGeneType :
        return self.orthologs[(gene_id_type, kinase_type, ambiguous)]
        
class OrthologManager:
        
    def __init__(self, orthology_dir : str, suffix : str = '_orthologs_final.tsv', debug : bool = False) :
        #get the list of files in the directory with the '_orthologs_final.tsv' extension
        self.ortholog_dir = orthology_dir
        self.ortholog_files = [f for f in os.listdir(orthology_dir) if f.endswith(suffix)]
        self.organism_list = [f.split('_')[0] for f in self.ortholog_files]
        
        self._species_ortho_dict = {}
        for organism in self.organism_list:
            print(organism)
            self._species_ortho_dict[organism] = SpeciesOrthologs(organism, os.path.join(orthology_dir, organism + suffix))
    
    def get_orthologs(self, species : str, gene_id_type : str, kinase_type : str, ambiguous : bool ) -> SpeciesOrthologs :
        return self._species_ortho_dict[species].get_orthologs(gene_id_type, kinase_type, ambiguous)
    
if __name__ == '__main__' :
    ortholog_manager = OrthologManager('orthologs')