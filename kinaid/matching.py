import pandas as pd
import math
from bisect import bisect_left
from typing import List
from tqdm.notebook import tqdm_notebook
import tqdm

class PWM_Matrices : 

    @staticmethod
    def score_peptide_with_matrix(sequence : str, pwm : dict, cols : list, favorability : dict = None, log_score : bool = False) -> float:
        if (len(sequence) - 1) != len(cols) :
            raise ValueError(f'Invalid sequence length. Expected {len(cols)} but got {len(sequence)}')
        flanks = sequence[:len(sequence) // 2] + sequence[len(sequence) // 2 + 1:]
        phosphorylated_site = sequence[len(sequence) // 2]
        if phosphorylated_site not in Scoring.__valid_phosphorylation_sites__ :
            raise ValueError('Invalid phosphorylated site')
        try :
            scores = [pwm[cols[i]][aa] for i,aa in enumerate(flanks) if (aa != '_' and aa != '-' and aa != 'X' and aa != 'x')]
        except KeyError as e:
            raise ValueError(f'Invalid amino acid: {e}')
        
        if favorability is not None :
            if phosphorylated_site == 'S' or phosphorylated_site == 'T' :
                scores.append(favorability[phosphorylated_site])
        
        if log_score :
            return sum([math.log(s) for s in scores])
        else :
            return math.prod(scores)
        
    def __init__(self, pwm_xlsx : str, ignore_suffix : str = '', debug : bool = False) :
        pwm_excel = pd.ExcelFile(pwm_xlsx, engine = 'openpyxl')
        if ignore_suffix != '' :
            self._kinase_names = [k for k in pwm_excel.sheet_names if not k.endswith(ignore_suffix)]
        else :
            self._kinase_names = pwm_excel.sheet_names
            
        pwm_dfs = {k:pwm_excel.parse(k) for k in self._kinase_names}
        
        assert(len(pwm_dfs[self._kinase_names[0]]) == l for l in (len(d) for d in pwm_dfs.values()))
        assert(all(d.columns[0] == 'AA' for d in pwm_dfs.values()))
        assert(all(pwm_dfs[self._kinase_names[0]]['AA'].equals(d['AA']) for d in pwm_dfs.values()))
        
        self._aa = list(pwm_dfs[self._kinase_names[0]]['AA'])
        
        any(d.set_index('AA', inplace=True) for d in pwm_dfs.values())
        self._pwm_dicts = {k:d.to_dict() for k,d in pwm_dfs.items()}
        self._pwm_dicts = {k:{int(c):aa for c,aa in cl.items()} for k,cl in self._pwm_dicts.items()}
        self._pwm_cols = {k:list(d.keys()) for k,d in self._pwm_dicts.items()}
        self._pwm_cols = {k:[int(c) for c in cl] for k,cl in self._pwm_cols.items()}
        any(c.sort() for c in self._pwm_cols.values())
        
        assert all(self._pwm_cols[self._kinase_names[0]] == c for c in self._pwm_cols.values())
        self._range = (self._pwm_cols[self._kinase_names[0]][0], self._pwm_cols[self._kinase_names[0]][-1])
        
        if debug :
            print(self._kinase_names)
            print(len(self._kinase_names))
            print(self._range)
            
        self._has_favorability = False
    
    def add_densitometry(self, densitometry_xlsx : str) :
        densitometry_excel = pd.ExcelFile(densitometry_xlsx, engine = 'openpyxl')
        densitometry_dfs = {k:densitometry_excel.parse(k) for k in self._kinase_names}
        
        assert(len(densitometry_dfs[self._kinase_names[0]]) == l for l in (len(d) for d in densitometry_dfs.values()))
        assert(all(d.columns[0] == 'AA' for d in densitometry_dfs.values()))
        assert(all(densitometry_dfs[self._kinase_names[0]]['AA'].equals(d['AA']) for d in densitometry_dfs.values()))
        
        
        any(d.set_index('AA', inplace=True) for d in densitometry_dfs.values())
        densitometry_dicts = {k:d.to_dict() for k,d in densitometry_dfs.items()}
        densitometry_dicts = {k:{int(c):aa for c,aa in cl.items()} for k,cl in densitometry_dicts.items()}
        densitometry_cols = {k:list(d.keys()) for k,d in densitometry_dicts.items()}
        densitometry_cols = {k:[int(c) for c in cl] for k,cl in densitometry_cols.items()}
        any(c.sort() for c in densitometry_cols.values())
        
        assert all(densitometry_cols[self._kinase_names[0]] == c for c in densitometry_cols.values())

        self._favorability_dicts = {}
  
        for k in self._kinase_names :
            col_map = densitometry_cols[k]
            S_S = sum([densitometry_dicts[k][c]['S'] for c in col_map])
            S_T = sum([densitometry_dicts[k][c]['T'] for c in col_map])
            S_ctrl = 0.75 * S_S - 0.25 * S_T
            T_ctrl = 0.75 * S_T - 0.25 * S_S
            S_0 = S_ctrl / (max(S_ctrl, T_ctrl))
            T_0 = T_ctrl / (max(S_ctrl, T_ctrl))
            self._favorability_dicts[k] = {'S' : S_0, 'T' : T_0}
        
        self._has_favorability = True

    def get_matrix_and_cols(self, kinase_name : str) :
        return self._pwm_dicts[kinase_name], self._pwm_cols[kinase_name]

    def get_favorability(self, kinase_name : str) :
        if not self._has_favorability :
            return None
        if kinase_name not in self._kinase_names :
            raise ValueError(f'Invalid kinase name: {kinase_name}')
        return self._favorability_dicts[kinase_name]
    
    def score_peptide(self, sequence : str, kinase : str, log_score : bool = False) -> float:
        return PWM_Matrices.score_peptide_with_matrix(sequence, *self.get_matrix_and_cols(kinase), self.get_favorability(kinase), log_score)
    
    def get_kinase_names(self) :
        return self._kinase_names
    
class Scoring:
    __valid_phosphorylation_sites__ = {'S', 'T', 'Y', 's', 't', 'y'}
    
    @staticmethod
    def get_sequence_format(sequence : str) -> str:
        # * format - * is right of the phosphorylated site and only exists once
        if '*' in sequence :
            split_sequence = sequence.split('*')
            if len(split_sequence) > 2 :
                raise ValueError('Invalid sequence format')
            if split_sequence[0][-1] not in Scoring.__valid_phosphorylation_sites__ :
                raise ValueError('Invalid sequence format')
            return 'star'
        # center format - the center position is the phosphorylated site
        elif sequence[len(sequence) // 2] in Scoring.__valid_phosphorylation_sites__:
            return 'center'
        else:
            raise ValueError('Invalid sequence format')
    

    def __init__(self, pwm_matrices : PWM_Matrices) :
        self._pwm_matrices = pwm_matrices
        self._kinase_names = pwm_matrices._kinase_names
        self._has_st_favorability = pwm_matrices._has_favorability
        self._range = pwm_matrices._range

    def clean_sequence(self, sequence : str, mode = None) -> str:
        if mode is None:
            mode = Scoring.get_sequence_format(sequence)
        
        split_sequence = ''
        phosphorylated_site = ''
        if mode == 'star':
            split_sequence = sequence.split('*')
            if len(split_sequence) > 2:
                raise ValueError('Invalid sequence format: too many stars')
            phosphorylated_site = split_sequence[0][-1]
            split_sequence[0] = split_sequence[0][:-1]

        elif mode == 'center':
            split_sequence = [sequence[:len(sequence) // 2], sequence[len(sequence) // 2 + 1:]]
            phosphorylated_site = sequence[len(sequence) // 2]



        if phosphorylated_site not in Scoring.__valid_phosphorylation_sites__:
            raise ValueError(f'Invalid phosphorylated site: {phosphorylated_site}')

        phosphorylated_site = phosphorylated_site.upper()

        left_side_len = 0 - self._range[0]
        right_side_len = self._range[1]

        left_side = '_' * (left_side_len - len(split_sequence[0])) + split_sequence[0] if len(split_sequence[0]) < left_side_len else split_sequence[0][-left_side_len:]

        right_side = split_sequence[1] + '_' * (right_side_len - len(split_sequence[1])) if len(split_sequence[1]) < right_side_len else split_sequence[1][:right_side_len]

        return left_side + phosphorylated_site + right_side

    
    def score_peptide(self, sequence : str, kinase: str, mode : str  = None, log_score : bool = False) -> float:
        if not mode == 'as_is' :
            cleaned_sequence = self.clean_sequence(sequence, mode)
        else :
            cleaned_sequence = sequence

        phosphorylation_site = cleaned_sequence[len(cleaned_sequence) // 2]

        if phosphorylation_site == 'Y' or phosphorylation_site == 'S' or phosphorylation_site == 'T' :
            return self._pwm_matrices.score_peptide(cleaned_sequence, kinase, log_score)
        else :
            raise ValueError('Invalid phosphorylated site')
        
    def get_kinase_names(self) :
        return self._kinase_names

class PeptideBackground :
    @staticmethod
    def background_factory(background_file : str,
                           sheet_name : str,
                           scoring : Scoring,
                           output_file : str,
                           progress : str = None) :
    
        #read the original background file
        background_df = pd.read_excel(background_file, engine='calamine', sheet_name=sheet_name, usecols=['SITE_+/-7_AA'])


        background_df['clean_seq'] = background_df['SITE_+/-7_AA'].apply(scoring.clean_sequence)

        background_df = background_df.drop(columns=['SITE_+/-7_AA'])
        
        kinases = scoring.get_kinase_names()

        kinases_loop = {None:kinases, 'notebook':tqdm_notebook(kinases), 'terminal':tqdm(kinases)}[progress]
        
        scores_df = pd.concat([
            pd.Series(
                background_df['clean_seq'].apply(
                    lambda x: scoring.score_peptide(
                        x,
                        kinase=k,
                        mode='as_is',
                        log_score=False
                    )
                ),
                index = background_df.index,
                name = k
            ) 
            for k in kinases_loop
        ], axis =1)
        scores_df = pd.concat([background_df, scores_df], axis=1)
        scores_df.to_csv(output_file, sep='\t', index=True, index_label='sequence')
        
        return PeptideBackground(output_file)
    
    def __init__(self, background_file : str) :
        df = pd.read_csv(background_file, sep = '\t', usecols=lambda x: x not in {'sequence', 'clean_seq'})
        self._background_size = len(df)
        kinase_names = list(df.columns)
        self._kinase_names = kinase_names
        self._background = {k:sorted(df[k].to_list()) for k in kinase_names}

    def get_percentile(self, score : float, kinase : str, low_score_skip : bool = False) :
        if low_score_skip and ((type(score) == float  or type(score) == int) and (score <= 0)) :
            return 0
        i = bisect_left(self._background[kinase], score)
        percentile = i/self._background_size * 100
        return percentile

class Match :
    def __init__(self, scoring : Scoring, background : PeptideBackground) :
        if not scoring._kinase_names == background._kinase_names :
            raise ValueError('Kinase names do not match between scoring and background')
        self._scoring = scoring
        self._background = background
        self._kinase_names = scoring._kinase_names

    def get_percentile(self, sequence : str, kinase : str, mode : str = None, low_score_skip : bool = False) :
        score = self._scoring.score_peptide(sequence, kinase, mode)
        return self._background.get_percentile(score, kinase, low_score_skip)
    
    def get_percentiles_for_kinase(self, sequences : List[str], kinase : str, mode : str = None, low_score_skip : bool = False) :
        return [self.get_percentile(s, kinase, mode, low_score_skip) for s in sequences]
    
    def get_percentiles(self, sequences : List[str], mode : str = None, low_score_skip : bool = False) :
        return {k:self.get_percentiles_for_kinase(sequences, k, mode, low_score_skip) for k in self._kinase_names}
    
    def get_percentiles_df(self, sequences : List[str], mode : str = None, low_score_skip : bool = False) :
        percentiles = self.get_percentiles(sequences, mode, low_score_skip)
        return pd.DataFrame(percentiles, index = sequences)
    
    def get_matched_kinases(self, sequence : str, threshold = 90.0, mode : str = None, low_score_skip : bool = False) :
        return [k for k in self._kinase_names if self.get_percentile(sequence, k, mode, low_score_skip) > threshold]
    
class MatchWithMapping :
    def __init__(self, matching : Match, mapping : dict) :
        self._matching = matching
        self._mapping = mapping
        self._mapped_kinase_names = list(mapping.values())
        if not all(k in matching._kinase_names for k in self._mapped_kinase_names) :
            raise ValueError('Invalid kinase names in mapping')
    
    def get_percentiles_for_kinase(self, sequences : List[str], kinase_symbol : str, mode : str = None, low_score_skip : bool = False) :
        return [self._matching.get_percentile(s, self._mapping[kinase_symbol], mode, low_score_skip) for s in sequences]
    
    def get_percentiles_for_all_kinases(self, sequences : List[str], mode : str = None, low_score_skip : bool = False) :
        return {sk:self.get_percentiles_for_kinase(sequences, hk, mode, low_score_skip) for sk,hk in self._mapping}