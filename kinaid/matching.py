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
            
        self._has_densitometry = False
    
    def add_densitometry(self, densitometry_xlsx : str) :
        densitometry_excel = pd.ExcelFile(densitometry_xlsx, engine = 'openpyxl')
        densitometry_dfs = {k:densitometry_excel.parse(k) for k in self._kinase_names}
        
        assert(len(densitometry_dfs[self._kinase_names[0]]) == l for l in (len(d) for d in densitometry_dfs.values()))
        assert(all(d.columns[0] == 'AA' for d in densitometry_dfs.values()))
        assert(all(densitometry_dfs[self._kinase_names[0]]['AA'].equals(d['AA']) for d in densitometry_dfs.values()))
        
        
        any(d.set_index('AA', inplace=True) for d in densitometry_dfs.values())
        self._densitometry_dicts = {k:d.to_dict() for k,d in densitometry_dfs.items()}
        self._densitometry_dicts = {k:{int(c):aa for c,aa in cl.items()} for k,cl in self._densitometry_dicts.items()}
        self._densitometry_cols = {k:list(d.keys()) for k,d in self._densitometry_dicts.items()}
        self._densitometry_cols = {k:[int(c) for c in cl] for k,cl in self._densitometry_cols.items()}
        any(c.sort() for c in self._densitometry_cols.values())
        
        assert all(self._densitometry_cols[self._kinase_names[0]] == c for c in self._densitometry_cols.values())
        self._densitometry_range = (self._densitometry_cols[self._kinase_names[0]][0], self._densitometry_cols[self._kinase_names[0]][-1])
        
        self._has_densitometry = True
        
class Scoring: 
    def __init__(self, st_pwm_matrices : PWM_Matrices, y_pwm_matrices : PWM_Matrices) :
        self._st_pwm_matrices = st_pwm_matrices
        self._y_pwm_matrices = y_pwm_matrices
        self._st_kinase_names = st_pwm_matrices._kinase_names
        self._y_kinase_names = y_pwm_matrices._kinase_names
        self._has_st_densitometry = st_pwm_matrices._has_densitometry