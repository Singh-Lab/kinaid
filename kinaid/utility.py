import requests
import time
import pandas as pd
import os

class Utility :
    MAPPING_API = 'https://rest.uniprot.org/idmapping/run/'

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
        ids = ','.join(list_of_ids)
        url = Utility.MAPPING_API
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
                    confidence = 'moderate') -> list:
        
        id = gene
        url = f'{url}/{in_taxon}/{id}/{out_taxon}/best_match'
        response = requests.post(url)

        map_confidence = {'high': 3, 'moderate': 2, 'low': 1}
        confidence_filter = map_confidence[confidence]
                        
        if not response.ok:
            print('Server responded:', response.status_code)
            print(response.text)
            return None
        else :
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
                    result_confidence = map_confidence[id_block[match]['confidence']]
                    result_symbol = id_block[match]['symbol']
                    if best_score and result_best_score == "No":
                        best_match = -1
                    if best_score_rev and result_best_score_rev == "No":
                        best_match = -1
                    match_confidence = result_confidence
                    if match_confidence < confidence_filter:
                        best_match = -1
                    if best_match > 0:
                        species_specific_geneid_type = None
                        species_specific_geneid = None
                        if species_specific:
                            species_specific_geneid_type = id_block[match]['species_specific_geneid_type']
                            species_specific_geneid = id_block[match]['species_specific_geneid']
                        matches.append((best_match, species_specific_geneid_type, species_specific_geneid, match_confidence, result_best_score == 'Yes', result_best_score_rev == 'Yes', result_symbol))
                if len(matches) == 0:
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
    def build_ortholog_database_for_organism(
            human_kinase_entrez_ids : set,
            organism_name : str,
            taxon_id : int,
            proteome_id : str,
            canonical_only : bool = False,
            best_score : bool = True,
            best_score_rev : bool = False,
            confidence : str = 'moderate',
            species_specific_geneid : bool = True,
            proteome_dir : str = 'proteomes') -> pd.DataFrame :
        Utility.download_proteome(organism_name, taxon_id, proteome_id, canonical_only, proteome_dir)

        filename = os.path.join(proteome_dir, f'{organism_name}_{str(taxon_id)}_{proteome_id}_{"con" if canonical_only else "noncon"}.tsv')

        entrez_id_dict = Utility.get_entrez_ids_of_proteome(filename)
        print(len(entrez_id_dict))

    def get_ortholog_in_human(tax_id:int, species_entrez_id:int, human_kinase_entrez_ids:set, best_score:bool=True, best_score_rev:bool=False, confidence:str='moderate', species_specific:bool=True) :
        __human_tax_id = '9606'
        matches = Utility.get_orthologs(species_entrez_id, in_taxon=tax_id, out_taxon=__human_tax_id, best_score=best_score, best_score_rev=best_score_rev, confidence=confidence, species_specific=species_specific)
        
        matched_human_kinase = [(species_entrez_id, ortholog,species_specific_geneid_type,
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
                                        result_symbol in matches if ortholog in human_kinase_entrez_ids]
        if len(matched_human_kinase) == 0:
            return None
        else:
            return matched_human_kinase
        
if __name__ == '__main__' :

    Utility.build_ortholog_database_for_organism({}, 'human', 9606, 'UP000005640')
    Utility.build_ortholog_database_for_organism({}, 'mouse', 10090, 'UP000000589')
    Utility.build_ortholog_database_for_organism({}, 'fly', 7227, 'UP000000803')
    Utility.build_ortholog_database_for_organism({}, 'worm', 6239, 'UP000001940')
    Utility.build_ortholog_database_for_organism({}, 'yeast', 4932, 'UP000002311')
    Utility.build_ortholog_database_for_organism({}, 'zebrafish', 7955, 'UP000000437')