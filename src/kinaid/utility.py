import requests
import time
import pandas as pd
import os
import numpy as np
from mpire import WorkerPool
from .matching import PWM_Matrices, Scoring, PeptideBackground
import argparse
from Bio import Entrez, SeqIO
from typing import List, Tuple, Set
from urllib.error import HTTPError, URLError
from tqdm import tqdm

tqdm.pandas()


Entrez.email = "javeda@princeton.edu"  # Replace with your email


class Utility :
    __MAPPING_API = 'https://rest.uniprot.org/idmapping/run/'
    __map_confidence = {'high': 3, 'moderate': 2, 'low': 1}

    __METHOD_WEIGHTS = {'Compara' : 0.93,
                        'Domainoid' : 1.0,
                        'ECGA' : 1.5,
                        'eggNOG' : 0.9,
                        'HGNC' : 1.5,
                        'Hieranoid' : 1.0,
                        'Homologene' : 1.0,
                        'Inparanoid' : 1.05,
                        'Isobase' : 0.95,
                        'OMA' : 1.01,
                        'OrthoDB' : 1.01,
                        'OrthoFinder' : 1.0,
                        'OrthoInspector' : 1.0,
                        'orthoMCL' : 0.9,
                        'Panther' : 1.1,
                        'Phylome' : 0.91,
                        'RoundUp' : 1.03,
                        'SGD' : 1.5,
                        'SonicParanoid' : 1.0,
                        'TreeFam' : 0.96,
                        'ZFIN' : 1.5,
                        'NCBI' : 1.0,
                        'SwiftOrtho' : 1.0}

    @staticmethod
    def download_file(url : str, filename : str):
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36'
            }
        with requests.get(url, headers=headers, stream=True) as response:
            response.raise_for_status()
            with open(filename, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192): 
                    if chunk: 
                        file.write(chunk)    
    
    @staticmethod
    def make_map_ids_job(list_of_ids : list, from_db : str ='UniProtKB_AC-ID', to_db : str = 'GeneID') -> str:
        if from_db == to_db:
            return None
        if from_db == 'UniProtKB' :
            from_db2 = 'UniProtKB_AC-ID'
        else :
            from_db2 = from_db
        
        ids = ','.join(list_of_ids)
        url = Utility.__MAPPING_API
        params = {'from': from_db2, 'to': to_db, 'ids': ids}
        response = requests.post(url, data=params)

        if not response.ok:
            print('UniProt job server responded:', response.status_code)
            print(response.text)
            return None
        else :
            job_id = response.text.split(':')[1].replace('\"', '').replace('}','')
            return job_id
    
    @staticmethod
    def get_map_ids_results(job_id : str,
                            url : str = 'https://rest.uniprot.org/idmapping/stream/',
                            exponential_backoff : bool = True,
                            max_tries : int = 10) -> dict:
        
        url = url + job_id

        for i in range(max_tries):
            response = None
            try: 
                response = requests.get(url)
            except requests.ConnectTimeout:
                print('Connection timeout')
            except requests.ConnectionError:
                print('Connection error with UniProt, waiting for 60 seconds')
                time.sleep(60)
            if response is None or not response.ok:
                if response is not None:
                    print('UniProt map server responded:', response.status_code)
                    print(response.text)
                if exponential_backoff:
                    time.sleep(2**i)
                else:
                    time.sleep(1)
            else :
                break
            print(f'Retry {i+1}')
        if response == None or not response.ok:
            print('Failed to get mapping results from UniProt, check connection')
            return {}
        
        result =  response.json()
        id_mapping_df = pd.DataFrame.from_records(data=result['results'])
        if(id_mapping_df.empty):
            return {}
        id_mapping_df.columns = ['from', 'to']
        id_mapping_group_df = id_mapping_df.groupby('from')['to'].apply(list)
        return id_mapping_group_df.to_dict()

    @staticmethod
    def get_orthologs(gene : str,
                    in_taxon : int = 9606,
                    out_taxon : int = 7227,
                    url : str = 'https://www.flyrnai.org/tools/diopt/web/diopt_api/v9/get_orthologs_from_entrez',
                    endpoint : str =  'none',
                    best_score : bool = False,
                    best_score_rev : bool = False,
                    species_specific : bool =True,
                    confidence = 'moderate',
                    exponential_backoff : bool = True,
                    max_tries : int = 10,
                    debug = False) -> list :
        id = str(gene)
        url = f'{url}/{in_taxon}/{id}/{out_taxon}/{endpoint}'
        if(debug):
            print(url)
        for i in range(max_tries) :
            response = None
            try: 
                response = requests.post(url)
            except requests.ConnectTimeout:
                print('Connection timeout')
            except requests.ConnectionError:
                print('Connection error with DIOPT, waiting for 60 seconds')
                time.sleep(60)
            if response == None or not response.ok:
                if response != None:
                    print('DIOPT server responded:', response.status_code)
                    print(response.text)
                if exponential_backoff:
                    time.sleep(2**i)
                else:
                    time.sleep(1)
            else :
                break
        if response is None or not response.ok:
            print('Failed to get orthologs from DIOPT, check connection')
            return [(-1, None, None, None, None, None, None, None)]

        confidence_filter = Utility.__map_confidence[confidence]
                        
        result =  response.json()
        if id not in result['results']:
            return [(-1, None, None, None, None, None, None, None)]
        else :
            id_block = result['results'][id]
            matches = []
            for match in id_block.keys() :
                best_match = int(match)
                result_best_score = id_block[match]['best_score']
                result_best_score_rev = id_block[match]['best_score_rev']
                result_confidence = Utility.__map_confidence[id_block[match]['confidence']]
                result_symbol = id_block[match]['symbol']
                #result_symbol = result['search_details']['gene_details'][0]['symbol']
                #if species_specific:
                #        species_specific_geneid_type = result['search_details']['gene_details'][0]['species_specific_geneid_type']
                #        species_specific_geneid = result['search_details']['gene_details'][0]['species_specific_geneid']
                if best_score and result_best_score == 'No':
                    best_match = -1
                    #print('score')
                if best_score_rev and result_best_score_rev == 'No':
                    #print('rev')
                    best_match = -1
                match_confidence = result_confidence
                if match_confidence < confidence_filter:
                    #print(f'confidence: {match_confidence}')
                    best_match = -1
                if best_match > 0:
                    species_specific_geneid_type = None
                    species_specific_geneid = None
                    if species_specific:
                        species_specific_geneid_type = id_block[match]['species_specific_geneid_type']
                        species_specific_geneid = str(id_block[match]['species_specific_geneid'])
                    weighted_score = sum([Utility.__METHOD_WEIGHTS[m] for m in id_block[match]['methods']])
                    matches.append((best_match, species_specific_geneid_type, species_specific_geneid, match_confidence, weighted_score, result_best_score == 'Yes', result_best_score_rev == 'Yes', result_symbol))
            if len(matches) == 0:
                #print('no matches')
                matches.append((-1, None, None, None, None, None, None, None))
            return matches

    
    @staticmethod
    def get_ortholog_from_human(
        tax_id : int,
        human_entrez_id : int,
        kinase_type : str,
        best_score : bool = False,
        best_score_rev : bool = False,
        confidence : str = 'moderate',
        species_specific : bool = True) -> List[Tuple[str, int, int, str, str, int, float, bool, bool, str]] :
        __human_tax_id = '9606'
        orthologs = Utility.get_orthologs(human_entrez_id,
                                          in_taxon=__human_tax_id,
                                          out_taxon=tax_id,
                                          best_score=best_score,
                                          best_score_rev=best_score_rev,
                                          confidence=confidence,
                                          species_specific=species_specific)
        
        results = [(kinase_type, human_entrez_id) + o for o in orthologs]
        
        if len(results) == 0:
            return None
        else:
            return results

    @staticmethod
    def __job2__(shared_objects, entrez_id) :
        in_taxon = shared_objects['in_taxon']
        out_taxon = shared_objects['out_taxon']
        best_score = shared_objects['best_score']
        best_score_rev = shared_objects['best_score_rev']
        confidence = shared_objects['confidence']
        species_specific = shared_objects['species_specific']
        endpoint = shared_objects['endpoint']
        
        results = Utility.get_orthologs(entrez_id,
                                        in_taxon=in_taxon,
                                        out_taxon=out_taxon,
                                        best_score=best_score,
                                        best_score_rev=best_score_rev,
                                        confidence=confidence,
                                        species_specific=species_specific,
                                        endpoint=endpoint)
        
        results = [o for o in results if o[0] != -1]
        results = [(entrez_id,) + o for o in results]
        if len(results) == 0:
            return None
        else :
            return results
        
    @staticmethod
    def get_kinases_with_high_confidence_not_in_list(query_ids : Set[str], kinase_ids : Set[str], in_taxon, out_taxon, threads = 8) :
        shared_objects = {'in_taxon' : in_taxon,
                          'out_taxon' : out_taxon,
                          'best_score' : True,
                          'best_score_rev' : True,
                          'confidence' : 'high',
                          'species_specific' : False,
                          'endpoint' : 'best_match'}
        with WorkerPool(threads, shared_objects=shared_objects, start_method='spawn') as pool:
            results = pool.map_unordered(Utility.__job2__, query_ids, progress_bar=True)
        
        #remove None values from the list
        results = [x for x in results if x is not None]
        results = [item for sublist in results for item in sublist]
        
        #remove result that has best_match in the kinase_ids
        results = [r for r in results if r[1] not in kinase_ids]
        
        #get the list of entrez_ids from results
        results = set([r[0] for r in results])
        
        return results
        
    @staticmethod
    def __job__(shared_objects, entrez_id, kinase_type) :
        taxon_id = shared_objects['taxon_id']
        best_score = shared_objects['best_score']
        best_score_rev = shared_objects['best_score_rev']
        confidence = shared_objects['confidence']
        species_specific = shared_objects['species_specific']
        return Utility.get_ortholog_from_human(taxon_id,
                                               entrez_id,
                                               kinase_type,
                                               best_score,
                                               best_score_rev,
                                               confidence,
                                               species_specific)
        
        
    @staticmethod
    def build_ortholog_database_for_organism(
        human_kinase_entrez_ids_st : set,
        human_kinase_entrez_ids_y : set,
        organism_name : str,
        taxon_id : int,
        best_score : bool = False,
        best_score_rev : bool = False,
        confidence : str = 'moderate',
        species_specific_geneid : bool = True,
        threads = 2) -> pd.DataFrame :
    
        shared_objects = {'taxon_id' : taxon_id,
                          'best_score' : best_score,
                          'best_score_rev' : best_score_rev,
                          'confidence' : confidence,
                          'species_specific' : species_specific_geneid}
        entrez_ids_st = [(e,'ST') for e in human_kinase_entrez_ids_st]
        entrez_ids_y = [(e,'Y') for e in human_kinase_entrez_ids_y]
        entrez_ids = entrez_ids_st + entrez_ids_y
        
        with WorkerPool(threads, shared_objects=shared_objects, start_method='spawn') as pool:
            results = pool.map_unordered(Utility.__job__, entrez_ids, progress_bar=True)

        #remove None values from the list
        results = [x for x in results if x is not None]

        #flatten the list
        results = [item for sublist in results for item in sublist]
        
        print('Finished for species: ' + organism_name)
        
        result_df = pd.DataFrame(results, columns=['kinase_type', 'human_entrez_id', 'species_entrez_id', 'species_specific_geneid_type', 'species_specific_geneid', 'match_confidence', 'weighted_score', 'result_best_score', 'result_best_score_rev', 'symbol'])
        
        result_df['species_entrez_id'] = result_df['species_entrez_id'].astype(np.int64) 
        result_df['human_entrez_id'] = result_df['human_entrez_id'].astype(np.int64)


        #get list of species_entrez_ids with confidence < 3
        low_confidence_ids = set(result_df[result_df['match_confidence'] < 3]['species_entrez_id'].unique())
        high_confidence_ids = set(result_df[result_df['match_confidence'] == 3]['species_entrez_id'].unique())
        
        low_confidence_ids = low_confidence_ids - high_confidence_ids
        
        to_be_removed = Utility.get_kinases_with_high_confidence_not_in_list(low_confidence_ids, {k for _,k in entrez_ids}, taxon_id, 9606, threads)
        
        result_df['should_be_removed'] = result_df['species_entrez_id'].isin(to_be_removed)
        
        
        
        #special cases for yeast
        if organism_name == 'yeast' :
            result_df['species_specific_geneid'] = result_df['species_specific_geneid'].apply(lambda x: f'SGD:{x}')
                
        return result_df   

    @staticmethod
    def geneid_to_proteins(gene_id : str) -> bool:
        
        try:
            handle = Entrez.elink(dbfrom='gene', db='protein', id=gene_id)
            record = Entrez.read(handle)
            handle.close()
        except (HTTPError, URLError) as e:
            print(f"Error fetching links for GeneID {gene_id}: {e}")
            return None

        if not record or 'LinkSetDb' not in record[0] or not record[0]['LinkSetDb']:
            return None

        protein_ids = []
        for linkset in record[0]['LinkSetDb']:
            if linkset['LinkName'] == 'gene_protein':
                protein_ids = [link['Id'] for link in linkset['Link']]
                break

        if not protein_ids:
            return None
        
        return protein_ids
    
    @staticmethod
    def geneid_is_pseudogene(gene_id : str) -> bool:
        proteins = Utility.geneid_to_proteins(gene_id)
        if proteins is None:
            return True
        return False
    
    @staticmethod
    def remove_pseudogenes(df: pd.DataFrame, debug=False) :
        df = df.copy()
        df['is_pseudogene'] = df['species_entrez_id'].progress_apply(Utility.geneid_is_pseudogene)
        df = df[~df['is_pseudogene']]
        return df
    
    @staticmethod
    def refactor_ortholog_file(ortholog_file : str,
                            entrez_to_human_kinase_dict : dict,
                            output_file : str
                        ) :
        df = pd.read_csv(ortholog_file, sep='\t')
        df = df[~df['should_be_removed']]
        df.drop(columns=['should_be_removed'], inplace=True)
        
        df['match_confidence'] = df.apply(lambda x: 2.5 if ((not x['result_best_score']) and x['result_best_score_rev']) else x['match_confidence'], axis=1)
        df['match_confidence'] = df.apply(lambda x: 2.25 if (x['result_best_score'] and (not x['result_best_score_rev'])) else x['match_confidence'], axis=1)
        
        unique_kinase_types = df['kinase_type'].unique()
        df['human_entrez_id'] = df['human_entrez_id'].astype(np.int64)
        df['species_entrez_id'] = df['species_entrez_id'].astype(np.int64)
        
        unique_kinase_types_human_dict = {}
        unique_kinase_types_species_dict = {}
        
        df['kinase_name'] = df['human_entrez_id'].map(entrez_to_human_kinase_dict)

        for unique_kinase_type in unique_kinase_types:
            unique_kinase_type_df = df[df['kinase_type'] == unique_kinase_type]
            unique_kinase_type_group_human_df = unique_kinase_type_df.groupby('kinase_name')['match_confidence'].max()
            unique_kinase_types_human_dict[unique_kinase_type] = unique_kinase_type_group_human_df.to_dict()
            unique_kinase_type_group_species_df = unique_kinase_type_df.groupby('species_entrez_id')['match_confidence'].max()
            unique_kinase_types_species_dict[unique_kinase_type] = unique_kinase_type_group_species_df.to_dict()

        df['max_confidence_human'] = df.apply(lambda row: unique_kinase_types_human_dict[row['kinase_type']][row['kinase_name']], axis=1)
        df['max_confidence_species'] = df.apply(lambda row: unique_kinase_types_species_dict[row['kinase_type']][row['species_entrez_id']], axis=1)

        df = df[df['match_confidence'] == df['max_confidence_human']]
        df = df[df['match_confidence'] == df['max_confidence_species']]
        
        df['gene_id_type'] = 'GeneID'
        df['gene_id'] = df['species_entrez_id']
        
        df.sort_values(by=['kinase_type', 'symbol','weighted_score'],
                       ascending = [True,True,False], inplace=True)
        
        #create a column called is_max_weighted_score that is True if the weighted_score is the max for the kinase_type and symbol
        df['is_max_weighted_score'] = df.groupby(['kinase_type', 'symbol'])['weighted_score'].transform('max') == df['weighted_score']
        
        df_group_by_max_weighted_score = df[df['is_max_weighted_score']].groupby(['kinase_type', 'symbol'])['weighted_score'].agg(list).reset_index()
        df_group_by_max_weighted_score['ambiguous'] = df_group_by_max_weighted_score.apply(lambda x: len(x['weighted_score']) > 1, axis=1)
        
        print(f'{ortholog_file} : {len(df_group_by_max_weighted_score[df_group_by_max_weighted_score["ambiguous"]])}')
        
        df_group = df.groupby(['kinase_type', 'symbol'])['kinase_name'].agg(list).reset_index()
        df_group['ambiguous'] = df_group.apply(lambda x: len(x['kinase_name']) > 1, axis=1)
        
        #make a nested dictionary with kinase_type and symbol as keys and if ambiguous as value
        ambiguous_df = df_group[['kinase_type', 'symbol', 'ambiguous']].copy()
        #ambiguous_df = ambiguous_df[ambiguous_df['ambiguous']]
        ambiguous_dict = ambiguous_df.set_index(['kinase_type', 'symbol'])['ambiguous'].to_dict()

        df.drop_duplicates(subset=['kinase_type', 'symbol'], keep='first', inplace=True)
        
        df['ambiguous'] = df.apply(lambda x: ambiguous_dict[(x['kinase_type'], x['symbol'])], axis=1)
        
        df_entrez = df[['kinase_name', 'kinase_type', 'gene_id_type', 'gene_id', 'symbol', 'ambiguous']].copy()
        df_specific = df[['kinase_name', 'kinase_type', 'species_specific_geneid_type', 'species_specific_geneid', 'symbol', 'ambiguous']].copy()
        
        #remove rows with NaN values in gene_id or geneid_type
        df_specific = df_specific.dropna(subset=['species_specific_geneid', 'species_specific_geneid_type'])
        
        #remove rows with empty strings in gene_id or geneid_type
        df_specific = df_specific[(df['species_specific_geneid'] != '') & (df_specific['species_specific_geneid_type'] != '')]
        
        #remove rows with '-' in gene_id or geneid_type
        df_specific = df_specific[(df['species_specific_geneid'] != '-') & (df_specific['species_specific_geneid_type'] != '-')]  
        
        #try to species_specific_geneid to np.int64, okay if it fails
        try :
            df_specific['species_specific_geneid'] = df_specific['species_specific_geneid'].astype(np.int64)
        except ValueError:
            pass         
        
        
        #rename species_specific_geneid_type to geneid_type and species_specific_geneid to gene_id in df_specific
        df_specific.rename(columns={'species_specific_geneid_type': 'gene_id_type', 'species_specific_geneid': 'gene_id'}, inplace=True)
        
        df_out = pd.concat([df_entrez, df_specific])     
        
        df_out.to_csv(output_file, sep='\t', index=False)
    
    """        
    @staticmethod
    def refactor_ortholog_file(ortholog_file : str,
                            organism : str,
                            entrez_to_human_kinase_dict : dict,
                            output_file : str,
                            match_confidence : str  = 'moderate',
                            result_best_score : bool = False,
                            result_best_score_rev : bool = False) :
        df = pd.read_csv(ortholog_file, sep='\t')

        confidence = Utility.__map_confidence[match_confidence]


        df = df[~df['should_be_removed']]
        df.drop(columns=['should_be_removed'], inplace=True)
        
        df['match_confidence'] = df.apply(lambda x: 2.5 if ((not x['result_best_score']) and x['result_best_score_rev']) else x['match_confidence'], axis=1)
        df['match_confidence'] = df.apply(lambda x: 2.25 if (x['result_best_score'] and (not x['result_best_score_rev'])) else x['match_confidence'], axis=1)

        unique_kinase_types = df['kinase_type'].unique()
        df['human_entrez_id'] = df['human_entrez_id'].astype(np.int64)
        df['species_entrez_id'] = df['species_entrez_id'].astype(np.int64)
        
        unique_kinase_types_human_dict = {}
        unique_kinase_types_species_dict = {}

        df['kinase_name'] = df['human_entrez_id'].map(entrez_to_human_kinase_dict)

        for unique_kinase_type in unique_kinase_types:
            unique_kinase_type_df = df[df['kinase_type'] == unique_kinase_type]
            unique_kinase_type_group_human_df = unique_kinase_type_df.groupby('human_entrez_id')['match_confidence'].max()
            unique_kinase_types_human_dict[unique_kinase_type] = unique_kinase_type_group_human_df.to_dict()
            unique_kinase_type_group_species_df = unique_kinase_type_df.groupby('species_entrez_id')['match_confidence'].max()
            unique_kinase_types_species_dict[unique_kinase_type] = unique_kinase_type_group_species_df.to_dict()

        df['max_confidence_human'] = df.apply(lambda row: unique_kinase_types_human_dict[row['kinase_type']][row['human_entrez_id']], axis=1)
        df['max_confidence_species'] = df.apply(lambda row: unique_kinase_types_species_dict[row['kinase_type']][row['species_entrez_id']], axis=1)

        #remove rows with confidence less than the max confidence
        df = df[df['match_confidence'] == df['max_confidence_human']]
        df = df[df['match_confidence'] == df['max_confidence_species']]

        df['kinase_name'] = df['human_entrez_id'].map(entrez_to_human_kinase_dict)
        df['organism'] = organism
        df['gene_id_type'] = 'GeneID'
        df['gene_id'] = df['species_entrez_id']
        df = df[df['match_confidence'] >= confidence]
        
        
        #remove duplicates with the same kinase_type, symbol, and kinase_name
        df = df.drop_duplicates(subset=['kinase_type', 'symbol', 'kinase_name'], keep='first')
        
        if result_best_score:
            df = df[df['result_best_score']]
        if result_best_score_rev:
            df = df[df['result_best_score_rev']]
            
        #df.drop(columns=['species_specific_geneid', 'human_entrez_id', 'species_specific_geneid_type', 'species_entrez_id', 'match_confidence', 'result_best_score', 'result_best_score_rev'], inplace=True)

        df.sort_values(by=['kinase_type', 'organism', 'gene_id_type', 'kinase_name', 'match_confidence', 'result_best_score_rev', 'symbol'], inplace=True)
        df_group = df.groupby(['kinase_type', 'symbol']).agg(list).reset_index()
        df_group['ambiguous'] = df_group.apply(lambda x: len(x['kinase_name']) > 1, axis=1)
        
        df = df_group.explode(['human_entrez_id',
                               'species_entrez_id',
                               'species_specific_geneid_type',
                               'species_specific_geneid',
                               'match_confidence',
                               'result_best_score',
                               'result_best_score_rev',
                               'max_confidence_human',
                               'max_confidence_species',
                               'kinase_name',
                               'organism',
                               'gene_id_type',
                               'gene_id']).reset_index(drop=True)
        
                         
        df_entrez = df[['organism', 'kinase_name', 'kinase_type', 'gene_id_type', 'gene_id', 'symbol', 'ambiguous']].copy()
        df_specific = df[['organism', 'kinase_name', 'kinase_type', 'species_specific_geneid_type', 'species_specific_geneid', 'symbol', 'ambiguous']].copy()
        
        #remove rows with NaN values in gene_id or geneid_type
        df_specific = df_specific.dropna(subset=['species_specific_geneid', 'species_specific_geneid_type'])
        
        #remove rows with empty strings in gene_id or geneid_type
        df_specific = df_specific[(df['species_specific_geneid'] != '') & (df_specific['species_specific_geneid_type'] != '')]
        
        #remove rows with '-' in gene_id or geneid_type
        df_specific = df_specific[(df['species_specific_geneid'] != '-') & (df_specific['species_specific_geneid_type'] != '-')]  
        
        #try to species_specific_geneid to np.int64, okay if it fails
        try :
            df_specific['species_specific_geneid'] = df_specific['species_specific_geneid'].astype(np.int64)
        except ValueError:
            pass         
        
        
        #rename species_specific_geneid_type to geneid_type and species_specific_geneid to gene_id in df_specific
        df_specific.rename(columns={'species_specific_geneid_type': 'gene_id_type', 'species_specific_geneid': 'gene_id'}, inplace=True)
        
        df_out = pd.concat([df_entrez, df_specific])     
        
        df_out.to_csv(output_file, sep='\t', index=False)
    """

    @staticmethod
    def add_supported_id_by_geneid(df, geneid_type, original_geneid_type = 'GeneID', unmapped_kinases = None) :

        if geneid_type == original_geneid_type :
            return df_copy, set()
        

        df_copy = df[df['gene_id_type'] == original_geneid_type].copy()
        #convert gene_id column to str
        if original_geneid_type == 'GeneID':
            df_copy['gene_id'] = df_copy['gene_id'].astype(np.int64)
        df_copy['gene_id'] = df_copy['gene_id'].astype(str)

        kinase_group_df = df_copy[['kinase_name', 'gene_id']].groupby('kinase_name')['gene_id'].apply(set)
        kinase_to_gene_id_dict = kinase_group_df.to_dict()
        


        if unmapped_kinases is None:
            unmapped_kinases = kinase_to_gene_id_dict.keys()

        genes_to_map = {g for k in unmapped_kinases if k in kinase_to_gene_id_dict for g in kinase_to_gene_id_dict[k]}

        gene_ids_job = Utility.make_map_ids_job(genes_to_map, original_geneid_type, geneid_type)
        time.sleep(3)
        if gene_ids_job is None:
            return df_copy, unmapped_kinases
        gene_ids_dict = Utility.get_map_ids_results(gene_ids_job)

        unmapped_kinases = {k for k in unmapped_kinases if k in kinase_to_gene_id_dict and all(g not in gene_ids_dict or len(gene_ids_dict[g]) == 0 for g in kinase_to_gene_id_dict[k])}
        
        df_new = df_copy.copy()
        del df_copy
        df_new['gene_id'] = df_new['gene_id'].map(gene_ids_dict)
        df_new = df_new[~df_new['gene_id'].isnull()]
        df_new['gene_id_type'] = geneid_type

        df_new = df_new.explode('gene_id').reset_index(drop=True)

        return df_new, unmapped_kinases
    
    @staticmethod
    def get_human_kinase_ids(ST_matrix_to_uniprot, Y_matrix_to_uniprot, remove_tyr = True, ST_matrices_exist : Set[str] = None, Y_matrices_exist : Set[str] = None) :
        st_kinases_df = pd.read_excel(ST_matrix_to_uniprot, sheet_name='Table S1 Data')
        st_kinase_dict = dict(zip(st_kinases_df['Matrix_name'], st_kinases_df['Uniprot id']))
        
        if ST_matrices_exist is not None:
            does_not_exist = {m for m in st_kinase_dict if m not in ST_matrices_exist}
            if(len(does_not_exist) > 0):
                print(f'Warning: The following ST matrices do not exist: {",".join(does_not_exist)}')
            st_kinase_dict = {k: v for k, v in st_kinase_dict.items() if k in ST_matrices_exist}
            
        y_kinases_df = pd.read_excel(Y_matrix_to_uniprot, sheet_name='Table_S1_Data')
        
        #remove rows that have entry in the SUBTYPE column
        y_kinases_df = y_kinases_df[~y_kinases_df['SUBTYPE'].isnull()]
        y_kinase_dict = dict(zip(y_kinases_df['MATRIX_NAME'], y_kinases_df['UNIPROT_ID']))
        
        if remove_tyr:
            y_kinase_dict = {k.split('_TYR')[0]: u for k, u in y_kinase_dict.items()}
            
        if Y_matrices_exist is not None:
            does_not_exist = {m for m in y_kinase_dict if m not in Y_matrices_exist}
            if(len(does_not_exist) > 0):
                print(f'Warning: The following Y matrices do not exist: {",".join(does_not_exist)}')
            y_kinase_dict = {k: v for k, v in y_kinase_dict.items() if k in Y_matrices_exist}
        
        st_kinases_uniprot = set(st_kinase_dict.values())
        y_kinases_uniprot = set(y_kinase_dict.values())
        
        return st_kinases_uniprot, y_kinases_uniprot, st_kinase_dict, y_kinase_dict
    
    @staticmethod
    def get_human_kinase_geneIDs(kinases_uniprot : set, kinase_dict : dict) :
        human_uniprot_to_entrez_job = Utility.make_map_ids_job(kinases_uniprot, 'UniProtKB_AC-ID', 'GeneID')
        time.sleep(5)
        human_kinases_uniprot_to_entrez_dict = Utility.get_map_ids_results(human_uniprot_to_entrez_job, max_tries=10)
        
        #print the differences between the uniprot ids and the ids that were mapped
        unmapped_ids = kinases_uniprot - set(human_kinases_uniprot_to_entrez_dict.keys())
        print(f'The following kinases were not mapped to a GeneID: {",".join(unmapped_ids)}')
        
        uniprot_to_human_kinase_dict = {v: k for k, v in kinase_dict.items()}
        
        human_kinases_entrez_to_uniprot_dict = {int(e): uniprot for uniprot, entrez_id_list in human_kinases_uniprot_to_entrez_dict.items() for e in entrez_id_list}

        human_entrez_ids = set(human_kinases_entrez_to_uniprot_dict.keys())
        
        return human_kinases_uniprot_to_entrez_dict, human_kinases_entrez_to_uniprot_dict, uniprot_to_human_kinase_dict, human_entrez_ids, unmapped_ids
    
    @staticmethod
    def build_human_kinase_database(data_dir : str, output_file : str, ST_matrices : PWM_Matrices = None, Y_matrices : PWM_Matrices = None) :
        ST_matrix_to_uniprot_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM3_ESM.xlsx'
        ST_matrix_to_uniprot = os.path.join(data_dir,'ST-Kinases_to_Uniprot.xlsx')

        Y_matrix_to_uniprot_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM3_ESM.xlsx'

        ST_matrix_to_uniprot = os.path.join(data_dir,'ST-Kinases_to_Uniprot.xlsx')
        Y_matrix_to_uniprot = os.path.join(data_dir,'Y-Kinases_to_Uniprot.xlsx')

        if not os.path.exists(ST_matrix_to_uniprot):
            Utility.download_file(ST_matrix_to_uniprot_url, ST_matrix_to_uniprot)

        if not os.path.exists(Y_matrix_to_uniprot):
            Utility.download_file(Y_matrix_to_uniprot_url, Y_matrix_to_uniprot)
        
        if ST_matrices is not None:
            ST_matrices_exist = set(ST_matrices.get_kinase_names())
        else:
            ST_matrices_exist = None
        
        if Y_matrices is not None:
            Y_matrices_exist = set(Y_matrices.get_kinase_names())
        else:
            Y_matrices_exist = None
        st_kinases_uniprot, y_kinases_uniprot, st_kinase_dict, y_kinase_dict = Utility.get_human_kinase_ids(ST_matrix_to_uniprot, Y_matrix_to_uniprot, ST_matrices_exist=ST_matrices_exist, Y_matrices_exist=Y_matrices_exist)
        human_kinases_uniprot_to_entrez_st_dict,_,_,_,unmapped_st = Utility.get_human_kinase_geneIDs(st_kinases_uniprot, st_kinase_dict)        
        human_kinases_uniprot_to_entrez_y_dict,_,_,_,unmapped_y = Utility.get_human_kinase_geneIDs(y_kinases_uniprot, y_kinase_dict)

        if len(human_kinases_uniprot_to_entrez_st_dict) == 0 or len(human_kinases_uniprot_to_entrez_y_dict) == 0:
            raise ValueError('No human kinases could be mapped from UniProt to GeneID: perhaps UniProt is down or the IDs are incorrect')
        
        st_kinase_geneIDs = [('ST', 'human', k, 'GeneID', g, k) for k, u in st_kinase_dict.items() if u not in unmapped_st for g in human_kinases_uniprot_to_entrez_st_dict[u]]
        st_kinase_UniProt = [('ST', 'human', k, 'UniProtKB', u, k) for k, u in st_kinase_dict.items()]
        
        y_kinase_geneIDs = [('Y', 'human', k, 'GeneID', g, k) for k, u in y_kinase_dict.items() if u not in unmapped_y for g in human_kinases_uniprot_to_entrez_y_dict[u]]
        y_kinase_UniProt = [('Y', 'human', k, 'UniProtKB', u, k) for k, u in y_kinase_dict.items()]
        
        df_data = st_kinase_geneIDs + y_kinase_geneIDs + st_kinase_UniProt + y_kinase_UniProt
        
        df = pd.DataFrame(df_data, columns=['kinase_type', 'organism', 'kinase_name', 'gene_id_type', 'gene_id', 'symbol'])
        
        #sort df by kinase_type, gene_id_type, and kinase_name
        df = df.sort_values(by=['kinase_type', 'gene_id_type', 'kinase_name'])
        df['short'] = df['symbol']
        df['long'] = df['symbol']
        df['ambiguous'] = False
        
        df.to_csv(output_file, sep='\t', index=False)
    
    @staticmethod
    def load_human_kinases_database(human_kinases_database_file : str, kinase_type : str = 'ST') :
        df = pd.read_csv(human_kinases_database_file, sep='\t')
        df = df[df['kinase_type'] == kinase_type]
        
        df_uniprot = df[df['gene_id_type'] == 'UniProtKB'].copy()
        df_entrez = df[df['gene_id_type'] == 'GeneID'].copy()
        df_entrez['gene_id'] = df_entrez['gene_id'].astype(np.int64)
        kinase_to_uniprot_dict = dict(zip(df_uniprot['kinase_name'], df_uniprot['gene_id']))
        entrez_to_kinase_dict = dict(zip(df_entrez['gene_id'], df_entrez['kinase_name']))
        
        return kinase_to_uniprot_dict, entrez_to_kinase_dict

    @staticmethod 
    def build_final_orthologs_database(df_final : pd.DataFrame, output_file:str) :


        df_final = df_final.sort_values(by=['gene_id_type', 'kinase_type', 'kinase_name'])
        
        df_final['symbol'] = df_final['symbol'].astype(str)
        df_final['ambiguous'] = df_final['ambiguous'].astype(bool)
        df_final['ambiguous_kinase'] = df_final['ambiguous']
        df_final.drop(columns=['ambiguous'], inplace=True)
        df_group = df_final.groupby(['kinase_type', 'gene_id_type', 'kinase_name']).agg(list).reset_index()
        
        #df_group['all_not_ambiguous'] = df_group.apply(lambda x: (not any(x['ambiguous']) and len(x['ambiguous']) > 1), axis=1)
        #df_group['symbol'] = df_group.apply(lambda x: [s for s, a in zip(x['symbol'], x['ambiguous']) if not a] if ((len(x['symbol']) > 1) and (x['all_not_ambiguous'])) else x['symbol'], axis=1)
        #df_group['gene_id'] = df_group.apply(lambda x: [g for g, a in zip(x['gene_id'], x['ambiguous']) if not a] if ((len(x['gene_id']) > 1) and (x['all_not_ambiguous'])) else x['gene_id'], axis=1)
        #df_group['ambiguous'] = df_group.apply(lambda x: [a for a in x['ambiguous'] if not a] if ((len(x['ambiguous']) > 1) and (x['all_not_ambiguous'])) else x['ambiguous'], axis=1)
        
        df_group['symbol_set'] = df_group['symbol'].apply(set)
        df_group['ambiguous_symbol'] = df_group['symbol_set'].apply(lambda x: len(x) > 1)
        df_group['symbol_set'] = df_group['symbol_set'].apply(list)
        df_group['symbol_set'] = df_group['symbol_set'].apply(sorted)
        df_group['short'] = df_group.apply(lambda x: f'~{x["symbol_set"][0]}-like' if x['ambiguous_symbol'] else x['symbol_set'][0], axis=1)
        df_group['long'] = df_group['symbol_set'].apply(lambda x: f'{",".join(sorted(x))}')
        #df_group['ambiguous'] = df_group['ambiguous'].apply(lambda x: all(x))
        #df_group.drop(columns=['symbol_set', 'all_not_ambiguous'], inplace=True)
        df_group.drop(columns=['symbol_set'], inplace=True)
        df_final = df_group.explode(['gene_id', 'symbol', 'ambiguous_kinase']).reset_index(drop=True)
        df_final['ambiguous'] = df_final['ambiguous_symbol'] | df_final['ambiguous_kinase']
        
        #output_file = output_file.replace('.tsv', '_2.tsv')
        df_final.to_csv(output_file, sep='\t', index=False)

    @staticmethod        
    def DefaultConfiguration(threads : int = 8, data_dir = './data', orthologs_dir = './orthologs') :
        
        if not os.path.exists(data_dir):
            print('Creating data directory')
            os.makedirs(data_dir)

        johnson_ST_matrices_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM4_ESM.xlsx'
        johnson_ST_matrices_original_file = os.path.join(data_dir,'johnson_ST_matrices.xlsx')

        if not os.path.exists(johnson_ST_matrices_original_file):
            print('Downloading ST matrices')
            Utility.download_file(johnson_ST_matrices_url, johnson_ST_matrices_original_file)
            
        johnson_ST_matrices_file = os.path.join(data_dir,'ST-Kinases.xlsx')
        if not os.path.exists(johnson_ST_matrices_file):
            print('Rearranging ST matrices')
            PWM_Matrices.rearrange_matrices(johnson_ST_matrices_original_file, sheet_name = 'ser_thr_all_norm_scaled_matrice', output_file=johnson_ST_matrices_file)

        densitometry_file = os.path.join(data_dir,'ST-Kinases_densitometry.xlsx')
        if not os.path.exists(densitometry_file):
            print('Rearranging densitometry matrices')
            PWM_Matrices.rearrange_matrices(johnson_ST_matrices_original_file, sheet_name = 'ser_thr_all_raw_matrices', output_file=densitometry_file)

        johnson_Y_matrices_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM4_ESM.xlsx'
        johnson_Y_matrices_original_file = os.path.join(data_dir,'johnson_Y_matrices.xlsx')

        if not os.path.exists(johnson_Y_matrices_original_file):
            print('Downloading Y matrices')
            Utility.download_file(johnson_Y_matrices_url, johnson_Y_matrices_original_file)

        johnson_Y_matrices_file = os.path.join(data_dir,'Y-Kinases.xlsx')
        if not os.path.exists(johnson_Y_matrices_file):
            print('Rearranging Y matrices')
            PWM_Matrices.rearrange_matrices(johnson_Y_matrices_original_file, sheet_name = 'tyrosine_all_norm_scaled_matric', pos = ['-5', '-4', '-3', '-2', '-1', '1', '2', '3', '4', '5'], output_file = johnson_Y_matrices_file, remove_suffix = '_TYR')

        ST_matrix_to_uniprot_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM3_ESM.xlsx'
        ST_matrix_to_uniprot = os.path.join(data_dir,'ST-Kinases_to_Uniprot.xlsx')

        if not os.path.exists(ST_matrix_to_uniprot):
            Utility.download_file(ST_matrix_to_uniprot_url, ST_matrix_to_uniprot)
            
        Y_matrix_to_uniprot_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM3_ESM.xlsx'
        Y_matrix_to_uniprot = os.path.join(data_dir,'Y-Kinases_to_Uniprot.xlsx')

        if not os.path.exists(Y_matrix_to_uniprot):
            Utility.download_file(Y_matrix_to_uniprot_url, Y_matrix_to_uniprot)

        ochoa_background_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM5_ESM.xlsx'


        ochoa_background_original_file = os.path.join(data_dir,'ochoa_background.xlsx')


        if not os.path.exists(ochoa_background_original_file):
            print('Downloading Ochoa background')
            Utility.download_file(ochoa_background_url, ochoa_background_original_file)


        tyrosine_background_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM5_ESM.xlsx'

        tyrosine_background_original_file = os.path.join(data_dir,'tyrosine_background.xlsx')

        if not os.path.exists(tyrosine_background_original_file):
            print('Downloading tyrosine background')
            Utility.download_file(tyrosine_background_url, tyrosine_background_original_file)

        print('Loading ST matrices')
        ST_matrices = PWM_Matrices(johnson_ST_matrices_file)
        #ST_matrices.add_densitometry(densitometry_file)
        
        print('Loading ST matrices (w/ favorability)')
        ST_matrices_wfav = PWM_Matrices(johnson_ST_matrices_file)
        ST_matrices_wfav.add_densitometry(densitometry_file)

        print('Loading Y matrices (w/ non-canonical)')
        Y_matrices = PWM_Matrices(johnson_Y_matrices_file)

        print('Creating scoring objects')
        st_scoring = Scoring(ST_matrices)
        st_scoring_wfav = Scoring(ST_matrices_wfav)
        y_scoring = Scoring(Y_matrices)
        
        ochoa_background_wfav_file = os.path.join(data_dir, 'johnson_ochoa_background_wfav.tsv')

        if not os.path.exists(ochoa_background_wfav_file) :
            print('Building Ochoa background (w/ favorability)')
            PeptideBackground.background_factory(ochoa_background_original_file, 'Supplementary Table 3', st_scoring_wfav, ochoa_background_wfav_file, progress='terminal')
        
        ochoa_background_file = os.path.join(data_dir, 'johnson_ochoa_background.tsv')
        if not os.path.exists(ochoa_background_file) :
            print('Building Ochoa background')
            PeptideBackground.background_factory(ochoa_background_original_file, 'Supplementary Table 3', st_scoring, ochoa_background_file, progress='terminal')

        tyrosine_background_original_file = os.path.join(data_dir,'tyrosine_background.xlsx')
        tyrosine_background_file = os.path.join(data_dir, 'johnson_tyrosine_background.tsv')

        if not os.path.exists(tyrosine_background_file) :
            #doesn't matter if canonical or not because only using peptide sequences
            print('Building tyrosine background')
            PeptideBackground.background_factory(tyrosine_background_original_file, 'Annotation - Canonical only', y_scoring, tyrosine_background_file, progress='terminal')


        human_kinases_database_file = os.path.join(data_dir, 'human_kinases_final.tsv')
        if not os.path.exists(human_kinases_database_file):
            print('Creating id mapping of human kinases')
            Utility.build_human_kinase_database(data_dir, human_kinases_database_file, ST_matrices, Y_matrices)
        
        
        kinase_to_uniprot_st_dict, entrez_to_kinase_st_dict  = Utility.load_human_kinases_database(human_kinases_database_file, 'ST')
        kinase_to_uniprot_y_dict, entrez_to_kinase_y_dict  = Utility.load_human_kinases_database(human_kinases_database_file, 'Y')
        
        st_kinases_uniprot = set(kinase_to_uniprot_st_dict.values())
        y_kinases_uniprot = set(kinase_to_uniprot_y_dict.values())
        
        human_entrez_st_ids = set(entrez_to_kinase_st_dict.keys())
        human_entrez_y_ids = set(entrez_to_kinase_y_dict.keys())
        
        #which kinases are in both sets
        dual_specificity_kinases = st_kinases_uniprot & y_kinases_uniprot

        print('Dual specificity kinases')
        print(dual_specificity_kinases)
        
        if not os.path.exists(orthologs_dir):
            os.makedirs(orthologs_dir)
        
        """
        arguments = [
                        ('mouse', 10090, 'UP000000589'),
                        ('fly', 7227, 'UP000000803'),
                        ('worm', 6239, 'UP000001940'),
                        ('yeast', 4932, 'UP000002311'),
                        ('zebrafish', 7955, 'UP000000437')
                    ]
        """
        arguments = [
                        ('mouse', 10090),
                        ('fly', 7227),
                        ('worm', 6239),
                        ('yeast', 4932),
                        ('zebrafish', 7955)
                    ] 
        #for organism_name, taxon_id, proteome_id in arguments :
        #    output_file = os.path.join(orthologs_dir, f'{organism_name}_{str(taxon_id)}_{proteome_id}_orthologs.tsv')
        for organism_name, taxon_id in arguments :
            output_file = os.path.join(orthologs_dir, f'{organism_name}_{str(taxon_id)}_orthologs.tsv')
            if not os.path.exists(output_file):
                print(f'Building ortholog database for {organism_name}')
                #df_final = Utility.build_ortholog_database_for_organism(human_entrez_st_ids, human_entrez_y_ids, organism_name, taxon_id, proteome_id,threads=threads)            
                df_final = Utility.build_ortholog_database_for_organism(human_entrez_st_ids, human_entrez_y_ids, organism_name, taxon_id, threads=threads)
                print('Remove pseudogenes')
                df_final = Utility.remove_pseudogenes(df_final)
                df_final.to_csv(output_file, sep='\t', index=False)
            else :
                print(f'Ortholog database for {organism_name} already exists')

        entrez_to_human_kinase_dict = {**entrez_to_kinase_st_dict, **entrez_to_kinase_y_dict}

        #for organism_name, taxon_id, proteome_id in arguments :
        for organism_name, taxon_id in arguments :
            ortholog_file = os.path.join(orthologs_dir, f'{organism_name}_{str(taxon_id)}_orthologs.tsv')
            #ortholog_file = os.path.join(orthologs_dir, f'{organism_name}_{str(taxon_id)}_{proteome_id}_orthologs.tsv')
            output_file = os.path.join(orthologs_dir, f'{organism_name}_orthologs_refactored.tsv')
            #Utility.refactor_ortholog_file(ortholog_file, organism_name, entrez_to_human_kinase_dict, output_file)
            Utility.refactor_ortholog_file(ortholog_file, entrez_to_human_kinase_dict, output_file)
        
        
        for organism_name, _, in arguments :
            output_file = os.path.join(orthologs_dir, f'{organism_name}_orthologs_final.tsv')
            if os.path.exists(output_file):
                print(f'Final ortholog database for {organism_name} already exists')
                continue
            
            refactored_file = os.path.join(orthologs_dir, f'{organism_name}_orthologs_refactored.tsv')
            df = pd.read_csv(refactored_file, sep='\t')
            
            unique_geneid_types = df['gene_id_type'].unique()
            specific_geneid_types = unique_geneid_types[unique_geneid_types != 'GeneID']
            
            df_uniprot, uniprot_unmapped = Utility.add_supported_id_by_geneid(df, 'UniProtKB', 'GeneID')
            print(f'UniProt unmapped ({organism_name}): {uniprot_unmapped}')
            for geneid_type in specific_geneid_types :
                if len(uniprot_unmapped) > 0:
                    df_uniprot2, uniprot_unmapped = Utility.add_supported_id_by_geneid(df, 'UniProtKB', geneid_type, uniprot_unmapped)
                    if(len(df_uniprot2) > 0):
                        df_uniprot = pd.concat([df_uniprot, df_uniprot2])
            
            df_final = pd.concat([df, df_uniprot])
            
            Utility.build_final_orthologs_database(df_final, output_file)

            #print(df_final.head())
            #sort by organism, gene_id_type, kinase_name
            """
            df_final = df_final.sort_values(by=['organism', 'gene_id_type', 'kinase_name'])
            
            df_final['symbol'] = df_final['symbol'].astype(str)
            df_group = df_final.groupby(['kinase_type', 'organism', 'gene_id_type', 'kinase_name']).agg(set).reset_index()
            df_group['symbol'] = df_group['symbol'].apply(lambda x: list(x))
            #df_group['symbol'] = df_group['symbol'].apply(lambda x: sorted(x))
            df_group['symbol'] = df_group.apply(lambda x: f'({x["kinase_name"]}-like){"+".join(x["symbol"])}' if len(x['symbol']) > 1 else x['symbol'][0], axis=1)
            df_final = df_group.explode('gene_id')
            
            df_final.to_csv(output_file, sep='\t', index=False)
            """
    @staticmethod
    def build_organism_orthologs(organism_name : str, taxon_id : int, orthologs_dir : str, human_kinases_database_file : str, threads : int = 8) :
        if not os.path.exists(human_kinases_database_file):
            print('Missing human kinases database file')
            print('Please run Utility.DefaultConfiguration()')
            return
        
        kinase_to_uniprot_st_dict, entrez_to_kinase_st_dict  = Utility.load_human_kinases_database(human_kinases_database_file, 'ST')
        kinase_to_uniprot_y_dict, entrez_to_kinase_y_dict  = Utility.load_human_kinases_database(human_kinases_database_file, 'Y')
        
        st_kinases_uniprot = set(kinase_to_uniprot_st_dict.values())
        y_kinases_uniprot = set(kinase_to_uniprot_y_dict.values())
        
        human_entrez_st_ids = set(entrez_to_kinase_st_dict.keys())
        human_entrez_y_ids = set(entrez_to_kinase_y_dict.keys())
        
        #which kinases are in both sets
        dual_specificity_kinases = st_kinases_uniprot & y_kinases_uniprot

        print('Dual specificity kinases')
        print(dual_specificity_kinases)
        
        if not os.path.exists(orthologs_dir):
            os.makedirs(orthologs_dir)
    
        output_file = os.path.join(orthologs_dir, f'{organism_name}_{str(taxon_id)}_orthologs.tsv')
        if not os.path.exists(output_file):
            print(f'Building ortholog database for {organism_name}')
            df_final = Utility.build_ortholog_database_for_organism(human_entrez_st_ids, human_entrez_y_ids, organism_name, taxon_id, threads=threads)
            print('Remove pseudogenes')
            df_final = Utility.remove_pseudogenes(df_final)
            df_final.to_csv(output_file, sep='\t', index=False)
        else :
            print(f'Ortholog database for {organism_name} already exists')
        
        entrez_to_human_kinase_dict = {**entrez_to_kinase_st_dict, **entrez_to_kinase_y_dict}
        
        ortholog_file = os.path.join(orthologs_dir, f'{organism_name}_{str(taxon_id)}_orthologs.tsv')
        output_file = os.path.join(orthologs_dir, f'{organism_name}_orthologs_refactored.tsv')
        Utility.refactor_ortholog_file(ortholog_file, entrez_to_human_kinase_dict, output_file)
        
        output_file = os.path.join(orthologs_dir, f'{organism_name}_orthologs_final.tsv')
        if os.path.exists(output_file):
            print(f'Final ortholog database for {organism_name} already exists')
            return
                
        refactored_file = os.path.join(orthologs_dir, f'{organism_name}_orthologs_refactored.tsv')
        df = pd.read_csv(refactored_file, sep='\t')
        
        unique_geneid_types = df['gene_id_type'].unique()
        specific_geneid_types = unique_geneid_types[unique_geneid_types != 'GeneID']
        
        df_uniprot, uniprot_unmapped = Utility.add_supported_id_by_geneid(df, 'UniProtKB', 'GeneID')
        print(f'UniProt unmapped ({organism_name}): {uniprot_unmapped}')
        for geneid_type in specific_geneid_types :
            if len(uniprot_unmapped) > 0:
                df_uniprot2, uniprot_unmapped = Utility.add_supported_id_by_geneid(df, 'UniProtKB', geneid_type, uniprot_unmapped)
                if(len(df_uniprot2) > 0):
                    df_uniprot = pd.concat([df_uniprot, df_uniprot2])
        
        df_final = pd.concat([df, df_uniprot])
        
        Utility.build_final_orthologs_database(df_final, output_file)
        
def run_default_configuration() :
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--threads', type=int, default=8, help='Number of threads to use')
    argparser.add_argument('--data_dir', type=str, default='./data', help='Directory to store data files')
    argparser.add_argument('--orthologs_dir', type=str, default='./orthologs', help='Directory to store ortholog files')
    args = argparser.parse_args()
    Utility.DefaultConfiguration(threads=args.threads, data_dir=args.data_dir, orthologs_dir=args.orthologs_dir)

def run_add_organism() :
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--organism_name', type=str, help='Name of organism')
    argparser.add_argument('--taxon_id', type=int, help='DIOPT Taxon ID of organism')
    argparser.add_argument('--orthologs_dir', type=str, default='./orthologs', help='Directory to store ortholog files')
    argparser.add_argument('--human_kinases_database_file', type=str, default='./data/human_kinases_final.tsv', help='File containing human kinases database')
    argparser.add_argument('--threads', type=int, default=8, help='Number of threads to use')
    args = argparser.parse_args()
    
    Utility.build_organism_orthologs(args.organism_name, args.taxon_id, args.orthologs_dir, args.human_kinases_database_file, threads=args.threads)
    
if __name__ == '__main__' :
    print('Starting default configuration')
    run_default_configuration()
    print('Configuration complete')