import pandas as pd
class PWM_Matrices : 
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
        #densitometry_range = (densitometry_cols[self._kinase_names[0]][0], densitometry_cols[self._kinase_names[0]][-1])

        favorability_dict = {}
        #S_S = math.prod([matrix[c]['S'] for c in col_map])
        #S_T = math.prod([matrix[c]['T'] for c in col_map])
        for k in self._kinase_names :
            col_map = densitometry_cols[k]
            S_S = sum([densitometry_dicts[k]['S'][c] for c in col_map])
            S_T = sum([densitometry_dicts[k]['T'][c] for c in col_map])
            S_ctrl = 0.75 * S_S - 0.25 * S_T
            T_ctrl = 0.75 * S_T - 0.25 * S_S
            S_0 = S_ctrl / (max(S_ctrl, T_ctrl))
            T_0 = T_ctrl / (max(S_ctrl, T_ctrl))
            favorability_dict[k] = {'S' : S_0, 'T' : T_0}
        
        self._has_favorability = True
        
class Scoring:
    __valid__phosphorylation_sites__ = {'S', 'T', 'Y', 's', 't', 'y'}
    
    @staticmethod
    def get_sequence_format(sequence : str) -> str:
        # * format - * is right of the phosphorylated site and only exists once
        if '*' in sequence :
            split_sequence = sequence.split('*')
            if len(split_sequence) > 2 :
                raise ValueError('Invalid sequence format')
            if split_sequence[0][-1] not in Scoring.__valid__phosphorylation_sites__ :
                raise ValueError('Invalid sequence format')
            return 'star'
        # center format - the center position is the phosphorylated site
        elif sequence[len(sequence) // 2] in Scoring.__valid__phosphorylation_sites__:
            return 'center'
        else:
            raise ValueError('Invalid sequence format')
    
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
            if phosphorylated_site.islower():
                phosphorylated_site = phosphorylated_site.upper()

        elif mode == 'center':
            split_sequence = [sequence[:len(sequence) // 2], sequence[len(sequence) // 2 + 1:]]
            phosphorylated_site = sequence[len(sequence) // 2]

        if phosphorylated_site == 'Y':
            left_side_len = 0 - self._y_range[0]
            right_side_len = self._y_range[1]
        elif phosphorylated_site == 'S' or phosphorylated_site == 'T':
            left_side_len = 0 - self._st_range[0]
            right_side_len = self._st_range[1]
        else:
            raise ValueError(f'Invalid phosphorylated site: {phosphorylated_site}')

        left_side = '_' * (left_side_len - len(split_sequence[0])) + split_sequence[0] if len(split_sequence[0]) < left_side_len else split_sequence[0][-left_side_len:]

        right_side = split_sequence[1] + '_' * (right_side_len - len(split_sequence[1])) if len(split_sequence[1]) < right_side_len else split_sequence[1][:right_side_len]

        #print(split_sequence)
        #print('left side len:', left_side_len)
        #print('right side len:', right_side_len)
        #print('phosphorylated site:', phosphorylated_site)
        #print('left side:', left_side)
        #print('right side:', right_side)


        return left_side + phosphorylated_site + right_side

    def __init__(self, st_pwm_matrices : PWM_Matrices, y_pwm_matrices : PWM_Matrices) :
        self._st_pwm_matrices = st_pwm_matrices
        self._y_pwm_matrices = y_pwm_matrices
        self._st_kinase_names = st_pwm_matrices._kinase_names
        self._y_kinase_names = y_pwm_matrices._kinase_names
        self._has_st_favorability = st_pwm_matrices._has_favorability
        self._st_range = st_pwm_matrices._range
        self._y_range = y_pwm_matrices._range
   