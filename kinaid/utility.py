import requests
import time
import pandas as pd
import os
import numpy as np
from mpire import WorkerPool

class Utility :
    __MAPPING_API = 'https://rest.uniprot.org/idmapping/run/'
    __map_confidence = {'high': 3, 'moderate': 2, 'low': 1}


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
    def rearrange_matrices (matrices_file : str = './out/johnson_ST_matrices.xlsx',
                            sheet_name : str = 'ser_thr_all_norm_scaled_matrice', 
                        output_file : str = './out/ST-Kinases.xlsx', 
                        pos = ['-5', '-4', '-3', '-2', '-1', '1', '2', '3', '4']):
        ae_df = pd.read_excel(matrices_file, engine='openpyxl', sheet_name=sheet_name)
        #rename first column to Kinase
        ae_df.rename(columns={ae_df.columns[0]: 'Kinase'}, inplace=True)
        ae_df.set_index('Kinase', inplace=True)

        res = ['P','G','A','C','S','T','V','I','L','M','F','Y','W','H','K','R','Q','N','D','E','s','t','y']

        kinase_matrices = {}
        for k,row in ae_df.iterrows() :
            probs = row.to_numpy()
            prob_matrix = np.reshape(probs, (len(pos),len(res)))
            prob_matrix_t = prob_matrix.transpose()
            kdf = pd.DataFrame(prob_matrix_t, columns=pos, index=res)
            kdf.index.name = 'AA'
            kinase_matrices[k] = kdf

        with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
            any(df.to_excel(writer, sheet_name=k) for k, df in kinase_matrices.items())
        
    
    @staticmethod
    def make_map_ids_job(list_of_ids : list, from_db : str ='UniProtKB_AC-ID', to_db : str = 'GeneID') -> str:
        ids = ','.join(list_of_ids)
        url = Utility.__MAPPING_API
        params = {'from': from_db, 'to': to_db, 'ids': ids}
        response = requests.post(url, data=params)

        if not response.ok:
            print('Server responded:', response.status_code)
            print(response.text)
            return None
        else :
            job_id = response.text.split(':')[1].replace('\"', '').replace('}','')
            return job_id
    
    @staticmethod
    def get_map_ids_results(job_id : str,
                            url : str = 'https://rest.uniprot.org/idmapping/stream/',
                            exponential_backoff : bool = True,
                            max_tries : int = 5) -> dict:
        
        url = url + job_id

        for i in range(max_tries):
            response = requests.get(url)
            if response.ok:
                result =  response.json()
                id_mapping_df = pd.DataFrame.from_records(data=result['results'])
                if(id_mapping_df.empty):
                    return {}
                id_mapping_df.columns = ['from', 'to']
                id_mapping_group_df = id_mapping_df.groupby('from')['to'].apply(list)
                return id_mapping_group_df.to_dict()
                
            elif exponential_backoff:
                time.sleep(2**i)
            else:
                time.sleep(1)

        print('Server responded:', response.status_code)
        print(response.text)
        return None

    @staticmethod
    def get_orthologs(gene : str,
                    in_taxon : int = 9606,
                    out_taxon : int = 7227,
                    url : str = 'https://www.flyrnai.org/tools/diopt/web/diopt_api/v9/get_orthologs_from_entrez',
                    best_score : bool = True,
                    best_score_rev : bool = True,
                    species_specific : bool =True,
                    confidence = 'moderate',
                    exponential_backoff : bool = True,
                    max_tries : int = 5) -> list :
        id = str(gene)
        url = f'{url}/{in_taxon}/{id}/{out_taxon}/best_match'
        for i in range(max_tries) :
            response = None
            try: 
                response = requests.post(url)
            except requests.ConnectTimeout:
                print('Connection timeout')
            except requests.ConnectionError:
                print('Connection error, waiting for 60 seconds')
                time.sleep(60)
            if response == None or not response.ok:
                if response != None:
                    print('Server responded:', response.status_code)
                    print(response.text)
                if exponential_backoff:
                    time.sleep(2**i)
                else:
                    time.sleep(1)
            else :
                break

        confidence_filter = Utility.__map_confidence[confidence]
                        
        result =  response.json()
        if id not in result['results']:
            return [(-1, None, None, None, None, None, None)]
        else :
            id_block = result['results'][id]
            matches = []
            for match in id_block.keys() :
                best_match = int(match)
                result_best_score = id_block[match]['best_score']
                result_best_score_rev = id_block[match]['best_score_rev']
                result_confidence = Utility.__map_confidence[id_block[match]['confidence']]
                #result_symbol = id_block[match]['symbol']
                result_symbol = result['search_details']['gene_details'][0]['symbol']
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
                        species_specific_geneid = id_block[match]['species_specific_geneid']
                    #print('match')
                    matches.append((best_match, species_specific_geneid_type, species_specific_geneid, match_confidence, result_best_score == 'Yes', result_best_score_rev == 'Yes', result_symbol))
            if len(matches) == 0:
                #print('no matches')
                matches.append((-1, None, None, None, None, None, None))
            return matches

    @staticmethod
    def download_proteome(organism_name : str,
                        taxon_id : int,
                        proteome_id : str,
                        canonical_only : bool = False,
                        proteome_dir : str = 'proteomes'):
        if not os.path.exists(proteome_dir):
            os.makedirs(proteome_dir)
        
        canonical = 'con' if canonical_only else 'noncon'
        filename = os.path.join(proteome_dir, f'{organism_name}_{str(taxon_id)}_{proteome_id}_{canonical}.tsv')
        
        if os.path.exists(filename):
            return filename
        
        if not canonical_only :
            url = f'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name&format=tsv&query=%28%28proteome%3A{proteome_id}%29%29'
        else :
            url = f'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cprotein_name&format=tsv&query=%28%28proteome%3A{proteome_id}%29+AND+reviewed%3Dtrue%29'
        
        Utility.download_file(url, filename)
        return filename

    @staticmethod
    def get_entrez_ids_of_proteome(proteome_file : str) -> dict:
        #open the file (.tsv) with the proteome of the species as a dataframe
        proteome_df = pd.read_csv(proteome_file, sep='\t')

        #get uniprot ids of the proteome into a list
        uniprot_ids = proteome_df['Entry'].tolist()

        #get the entrez ids of the uniprot ids
        uniprot_to_entrez_job = Utility.make_map_ids_job(uniprot_ids, 'UniProtKB_AC-ID', 'GeneID')
        uniprot_to_entrez_dict = Utility.get_map_ids_results(uniprot_to_entrez_job, max_tries=10)

        #create a dict from entrez id to uniprot id
        entrez_to_uniprot_dict = {int(e): uniprot for uniprot, entrez_id_list in uniprot_to_entrez_dict.items() for e in entrez_id_list}

        return entrez_to_uniprot_dict

    @staticmethod
    def get_ortholog_in_human(tax_id:int, species_entrez_id:int, human_kinase_entrez_ids_st:set, human_kinase_entrez_ids_y:set, best_score:bool=True, best_score_rev:bool=False, confidence:str='moderate', species_specific:bool=True) :
        __human_tax_id = '9606'
        matches = Utility.get_orthologs(species_entrez_id, in_taxon=tax_id, out_taxon=__human_tax_id, best_score=best_score, best_score_rev=best_score_rev, confidence=confidence, species_specific=species_specific)
        matched_human_kinase_st = [('ST', species_entrez_id, ortholog, species_specific_geneid_type,
                                species_specific_geneid,
                                match_confidence,
                                    result_best_score,
                                    result_best_score_rev,
                                        result_symbol) for ortholog, 
                            species_specific_geneid_type,
                                species_specific_geneid,
                                match_confidence,
                                    result_best_score,
                                    result_best_score_rev,
                                        result_symbol in matches if ortholog in human_kinase_entrez_ids_st]
        matched_human_kinase_y = [('Y', species_entrez_id, ortholog, species_specific_geneid_type,
                                species_specific_geneid,
                                match_confidence,
                                    result_best_score,
                                    result_best_score_rev,
                                        result_symbol) for ortholog, 
                            species_specific_geneid_type,
                                species_specific_geneid,
                                match_confidence,
                                    result_best_score,
                                    result_best_score_rev,
                                        result_symbol in matches if ortholog in human_kinase_entrez_ids_y]
        
        matched_human_kinase = matched_human_kinase_st + matched_human_kinase_y
        
        if len(matched_human_kinase) == 0:
            return None
        else:
            return matched_human_kinase
            
    @staticmethod
    def job(shared_objects, entrez_id) :
        taxon_id = shared_objects['taxon_id']
        human_kinase_entrez_ids_st = shared_objects['human_kinase_entrez_ids_st']
        human_kinase_entrez_ids_y = shared_objects['human_kinase_entrez_ids_y']
        best_score = shared_objects['best_score']
        best_score_rev = shared_objects['best_score_rev']
        confidence = shared_objects['confidence']
        species_specific = shared_objects['species_specific']
        return Utility.get_ortholog_in_human(taxon_id, entrez_id, human_kinase_entrez_ids_st, human_kinase_entrez_ids_y, best_score, best_score_rev, confidence, species_specific)
    
    @staticmethod
    def build_ortholog_database_for_organism(
            human_kinase_entrez_ids_st : set,
            human_kinase_entrez_ids_y : set,
            organism_name : str,
            taxon_id : int,
            proteome_id : str,
            canonical_only : bool = False,
            best_score : bool = True,
            best_score_rev : bool = False,
            confidence : str = 'moderate',
            species_specific_geneid : bool = True,
            proteome_dir : str = 'proteomes',
            threads = 2) -> pd.DataFrame :
        Utility.download_proteome(organism_name, taxon_id, proteome_id, canonical_only, proteome_dir)

        filename = os.path.join(proteome_dir, f'{organism_name}_{str(taxon_id)}_{proteome_id}_{"con" if canonical_only else "noncon"}.tsv')

        entrez_id_dict = Utility.get_entrez_ids_of_proteome(filename)
        
        organism_entrez_ids = set(entrez_id_dict.keys())
        
        shared_objects = {'taxon_id' : taxon_id,
                          'human_kinase_entrez_ids_st' : human_kinase_entrez_ids_st,
                          'human_kinase_entrez_ids_y' : human_kinase_entrez_ids_y,
                          'best_score' : best_score,
                          'best_score_rev' : best_score_rev,
                          'confidence' : confidence,
                          'species_specific' : species_specific_geneid}
        
        with WorkerPool(threads, shared_objects=shared_objects, start_method='fork') as pool:
            results = pool.map_unordered(Utility.job, organism_entrez_ids, progress_bar=True)

        #remove None values from the list
        results = [x for x in results if x is not None]

        #flatten the list
        results = [item for sublist in results for item in sublist]
        
        print('Finished for species: ' + organism_name)
        
        result_df = pd.DataFrame(results, columns=['kinase_type', 'species_entrez_id', 'human_entrez_id', 'species_specific_geneid_type', 'species_specific_geneid', 'match_confidence', 'result_best_score', 'result_best_score_rev', 'result_symbol'])
        
        result_df['species_entrez_id'] = result_df['species_entrez_id'].astype(int) 
        result_df['human_entrez_id'] = result_df['human_entrez_id'].astype(int)
                
        return result_df
        
        #result_df.to_csv(output_file, sep='\t', index=False)

    @staticmethod
    def refactor_ortholog_file(ortholog_file : str,
                            organism : str,
                            entrez_to_human_kinase_dict : dict,
                            output_file : str,
                            match_confidence : str  = 'moderate',
                            result_best_score : bool = True,
                            result_best_score_rev : bool = True) :
        df = pd.read_csv(ortholog_file, sep='\t')
        df['kinase_name'] = df['human_entrez_id'].map(entrez_to_human_kinase_dict)
        df['organism'] = organism
        df['geneid_type'] = 'GeneID'
        df['gene_id'] = df['species_entrez_id']
        confidence = Utility.__map_confidence[match_confidence]
        df = df[df['match_confidence'] >= confidence]
        if result_best_score:
            df = df[df['result_best_score']]
        if result_best_score_rev:
            df = df[df['result_best_score_rev']]
        df.drop(columns=['species_specific_geneid', 'human_entrez_id', 'species_specific_geneid_type', 'species_entrez_id', 'match_confidence', 'result_best_score', 'result_best_score_rev'], inplace=True)

        df = df[['organism', 'kinase_name', 'kinase_type', 'geneid_type', 'gene_id', 'result_symbol']]
        df.to_csv(output_file, sep='\t', index=False)

        
if __name__ == '__main__' :
    data_dir = './data'
    threads = 8

    johnson_ST_matrices_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM4_ESM.xlsx'
    johnson_ST_matrices_original_file = os.path.join(data_dir,'johnson_ST_matrices.xlsx')

    if not os.path.exists(johnson_ST_matrices_original_file):
        Utility.download_file(johnson_ST_matrices_url, johnson_ST_matrices_original_file)
        
    johnson_ST_matrices_file = os.path.join(data_dir,'ST-Kinases.xlsx')
    Utility.rearrange_matrices(johnson_ST_matrices_original_file, sheet_name = 'ser_thr_all_norm_scaled_matrice', output_file=johnson_ST_matrices_file)

    densitometry_file = os.path.join(data_dir,'ST-Kinases_densitometry.xlsx')
    Utility.rearrange_matrices(johnson_ST_matrices_original_file, sheet_name = 'ser_thr_all_raw_matrices', output_file=densitometry_file)

    johnson_Y_matrices_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM4_ESM.xlsx'
    johnson_Y_matrices_original_file = os.path.join(data_dir,'johnson_Y_matrices.xlsx')

    if not os.path.exists(johnson_Y_matrices_original_file):
        Utility.download_file(johnson_Y_matrices_url, johnson_Y_matrices_original_file)

    johnson_Y_matrices_file = os.path.join(data_dir,'Y-Kinases.xlsx')
    Utility.rearrange_matrices(johnson_Y_matrices_original_file, sheet_name = 'tyrosine_all_norm_scaled_matric', pos = ['-5', '-4', '-3', '-2', '-1', '1', '2', '3', '4', '5'], output_file = johnson_Y_matrices_file)

    ST_matrix_to_uniprot_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM3_ESM.xlsx'
    ST_matrix_to_uniprot = os.path.join(data_dir,'ST-Kinases_to_Uniprot.xlsx')

    Y_matrix_to_uniprot_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM3_ESM.xlsx'

    ST_matrix_to_uniprot = os.path.join(data_dir,'ST-Kinases_to_Uniprot.xlsx')
    Y_matrix_to_uniprot = os.path.join(data_dir,'Y-Kinases_to_Uniprot.xlsx')

    if not os.path.exists(ST_matrix_to_uniprot):
        Utility.download_file(ST_matrix_to_uniprot_url, ST_matrix_to_uniprot)

    if not os.path.exists(Y_matrix_to_uniprot):
        Utility.download_file(Y_matrix_to_uniprot_url, Y_matrix_to_uniprot)

    ochoa_background_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05575-3/MediaObjects/41586_2022_5575_MOESM5_ESM.xlsx'


    ochoa_background_original_file = os.path.join(data_dir,'ochoa_background.xlsx')


    if not os.path.exists(ochoa_background_original_file):
        Utility.download_file(ochoa_background_url, ochoa_background_original_file)


    tyrosine_background_url = 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07407-y/MediaObjects/41586_2024_7407_MOESM5_ESM.xlsx'

    tyrosine_background_original_file = os.path.join(data_dir,'tyrosine_background.xlsx')

    if not os.path.exists(tyrosine_background_original_file):
        Utility.download_file(tyrosine_background_url, tyrosine_background_original_file)

    st_kinases_df = pd.read_excel(ST_matrix_to_uniprot, sheet_name='Table S1 Data')
    st_kinase_dict = dict(zip(st_kinases_df['Matrix_name'], st_kinases_df['Uniprot id']))
    
    y_kinases_df = pd.read_excel(Y_matrix_to_uniprot, sheet_name='Table_S1_Data')
    
    #remove rows that have entry in the SUBTYPE column
    y_kinases_df = y_kinases_df[~y_kinases_df['SUBTYPE'].isnull()]
    y_kinase_dict = dict(zip(y_kinases_df['MATRIX_NAME'], y_kinases_df['UNIPROT_ID']))
    
    st_kinases_uniprot = set(st_kinase_dict.values())
    y_kinases_uniprot = set(y_kinase_dict.values())

    #which kinases are in both sets
    dual_specificity_kinases = st_kinases_uniprot & y_kinases_uniprot

    print('Dual specificity kinases')
    print(dual_specificity_kinases)
    
    #human_kinase_dict = {**st_kinase_dict, **y_kinase_dict}

    
    #human_kinases_uniprot = st_kinases_uniprot | y_kinases_uniprot
    #human_uniprot_to_entrez_job = Utility.make_map_ids_job(human_kinases_uniprot, 'UniProtKB_AC-ID', 'GeneID')
    
    #human_kinases_uniprot_to_entrez_dict = Utility.get_map_ids_results(human_uniprot_to_entrez_job)
    
    #for uniprot, entrez_id_list in human_kinases_uniprot_to_entrez_dict.items():
    #    if len(entrez_id_list) != 1:
    #        print('Warning: entrez id list size not 1')
    #        print(uniprot)
    #        print(entrez_id_list)
            
    #uniprot_to_human_kinase = {v: k for k, v in human_kinase_dict.items()}
    #human_entrez_to_uniprot_dict = {int(e): uniprot for uniprot, entrez_id_list in human_kinases_uniprot_to_entrez_dict.items() for e in entrez_id_list}

    #human_entrez_ids = set(human_entrez_to_uniprot_dict.keys())

    human_uniprot_to_entrez_st_job = Utility.make_map_ids_job(st_kinases_uniprot, 'UniProtKB_AC-ID', 'GeneID')
    human_uniprot_to_entrez_y_job = Utility.make_map_ids_job(y_kinases_uniprot, 'UniProtKB_AC-ID', 'GeneID')
    
    human_kinases_uniprot_to_entrez_st_dict = Utility.get_map_ids_results(human_uniprot_to_entrez_st_job)
    human_kinases_uniprot_to_entrez_y_dict = Utility.get_map_ids_results(human_uniprot_to_entrez_y_job)
    
    uniprot_to_human_kinase_st_dict = {v: k for k, v in st_kinase_dict.items()}
    uniprot_to_human_kinase_y_dict = {v: k for k, v in y_kinase_dict.items()}
    
    human_kinases_entrez_to_uniprot_st_dict = {int(e): uniprot for uniprot, entrez_id_list in human_kinases_uniprot_to_entrez_st_dict.items() for e in entrez_id_list}
    human_kinases_entrez_to_uniprot_y_dict = {int(e): uniprot for uniprot, entrez_id_list in human_kinases_uniprot_to_entrez_y_dict.items() for e in entrez_id_list}

    human_entrez_st_ids = set(human_kinases_entrez_to_uniprot_st_dict.keys())
    human_entrez_y_ids = set(human_kinases_entrez_to_uniprot_y_dict.keys())
    
    orthologs_dir = 'orthologs'
    
    if not os.path.exists(orthologs_dir):
        os.makedirs(orthologs_dir)
    
    
    #Utility.build_ortholog_database_for_organism(human_entrez_ids, 'human', 9606, 'UP000005640')
    arguments = [
                    ('mouse', 10090, 'UP000000589'),
                    ('fly', 7227, 'UP000000803'),
                    ('worm', 6239, 'UP000001940'),
                    ('yeast', 4932, 'UP000002311'),
                    ('zebrafish', 7955, 'UP000000437')
                ]
    
    for organism_name, taxon_id, proteome_id in arguments :
        output_file = os.path.join(orthologs_dir, f'{organism_name}_{str(taxon_id)}_{proteome_id}_orthologs.tsv')
        if not os.path.exists(output_file):
            print(f'Building ortholog database for {organism_name}')
            df_final =Utility.build_ortholog_database_for_organism(human_entrez_st_ids, human_entrez_y_ids, organism_name, taxon_id, proteome_id,threads=threads)            
            df_final.to_csv(output_file, sep='\t', index=False)
        else :
            print(f'Ortholog database for {organism_name} already exists')

    uniprot_to_human_kinase_st_dict = {u : k for k,u in st_kinase_dict.items()}
    uniprot_to_human_kinase_y_dict = {u : k for k,u in y_kinase_dict.items()}

    entrez_to_human_kinase_st_dict = {e : uniprot_to_human_kinase_st_dict[u] for e,u in human_kinases_entrez_to_uniprot_st_dict.items()}
    entrez_to_human_kinase_y_dict = {e : uniprot_to_human_kinase_y_dict[u] for e,u in human_kinases_entrez_to_uniprot_y_dict.items()}

    entrez_to_human_kinase_dict = {**entrez_to_human_kinase_st_dict, **entrez_to_human_kinase_y_dict}

    for organism_name, taxon_id, proteome_id in arguments :
        ortholog_file = os.path.join(orthologs_dir, f'{organism_name}_{str(taxon_id)}_{proteome_id}_orthologs.tsv')
        output_file = os.path.join(orthologs_dir, f'{organism_name}_orthologs_refactored.tsv')
        Utility.refactor_ortholog_file(ortholog_file, organism_name, entrez_to_human_kinase_dict, output_file)
