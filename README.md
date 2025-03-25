## KINAID
<img src="https://kinaid.princeton.edu/assets/logo3.png" width="500">


# Public webserver version
https://kinaid.princeton.edu

# Install library
  ```
  python3 -m pip install kinaid
  ```

### install resources
can take up to 3GB HDD memory.

## Requires the following libraries

- `biopython`
- `dash`
- `dash-cytoscape`
- `mpire`
- `numpy`
- `openpyxl`
- `pandas`
- `pyxlsb`
- `scipy`
- `statsmodels`

  ```
  cd <DIRECTORY TO INSTALL KINAID RESOURCES>
  kinaid-install
  ```

### add additional organism
Choose organism from DIOPT
https://www.flyrnai.org/diopt

```
usage: kinaid-add [-h] [--organism_name ORGANISM_NAME] [--taxon_id TAXON_ID] [--orthologs_dir ORTHOLOGS_DIR]
                  [--human_kinases_database_file HUMAN_KINASES_DATABASE_FILE] [--threads THREADS]

options:
  -h, --help            show this help message and exit
  --organism_name ORGANISM_NAME
                        Name of organism
  --taxon_id TAXON_ID   DIOPT Taxon ID of organism
  --orthologs_dir ORTHOLOGS_DIR
                        Directory to store ortholog files
  --human_kinases_database_file HUMAN_KINASES_DATABASE_FILE
                        File containing human kinases database
  --threads THREADS     Number of threads to use
```

#### Example
```
kinaid-add --organism_name clawed_frog --taxon_id 8364 --threads 4
```


# Private webserver version

- clone this repository
- change to kinaid root directory

### install resources
can take up to 3GB HDD memory.

  ```
  python3 -m kinaid.utility --threads 4
  ```
### run dashboard
  ```
  usage: dashboard.py [-h] [--port PORT] [--host HOST] [--no_debug]

  options:
  -h, --help   show this help message and exit
  --port PORT
  --host HOST
  --no_debug
  ```

  ```
  python3 dashboard.py
  ```

## Examples

### sandbox.ipynb

basic examples of loading phosphoproteomics data and running library version of kinaid

### yeast_tests.ipynb

notebook to generate yeast experiment data for Supplementary Table 4. The data must be loaded by the dashboard.

### offline_tests.ipynb

notebook to generate and run experiments for Supplementary Table 3. This experiment evaluated the performance of the library and not the dashboard.
