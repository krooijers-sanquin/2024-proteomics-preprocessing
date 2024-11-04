## Contents

The goal is to identify the minimal (parsimonious) set of proteins that is responsible to generate a set of peptides, observed across a set of samples.

A significant part of this script is based on [this StackOverflow answer](https://stackoverflow.com/a/50027975).
It has a good explanation and visualization of the bipartite graph. In our case, proteins are nodes on the left, peptides are nodes on the right. 

A lot of the "algorithmic heavy lifting" is done by the [NetworkX max_flow_min_cost function](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.flow.max_flow_min_cost.html). 
A twist in this script, is that we use an additional set of weights, such that curated Uniprot entries (from "SwissProt"), as well as entries for which there is strong evidence of protein existence, are favored over others.

See `Steps to reproduce` below, to see how to get this running. 


## Scripts

- `scripts/merge_peptide_tables.py`
  - The main script which assigns peptides to genes
- `scripts/convert_uniprot_xml_to_duckdb.py`
  - Script to convert Uniprot XML entries to a database which is used for protein sequence lookup (amongst others)
  - *NB: internally, we create larger databases which also combine data from MSigDB and HUGO, termed the OmniDB*
- `scripts/extract_uniprot_entries_subset.py`
  - Script to extract a small subset of entries from the Uniprot XML dump, to speed up testing


## Steps to reproduce

### Setup python/poetry environment

*(I assume familiarity with python, poetry and pyenv)*

- `pyenv install 3.12.7`
- `pyenv local 3.12.7`
- `pyenv shell 3.12.7`
- `poetry init -n`
- `poetry env use $(pyenv which python)`
- `poetry install --no-root`

### Obtain Uniprot proteome dump

See:
- <https://www.uniprot.org/proteomes/UP000005640>
- <https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.xml.gz>

(NB: I named the obtained file `KR20231004.proteomes.UP000005640.xml.gz`)

### Subset Uniprot XML

*(just for testing purposes)*

- `poetry run python scripts/extract_uniprot_entries_subset.py KR20231004.proteomes.UP000005640.xml.gz tests/data/uniprot_5640_subset.xml`
- `gzip tests/data/uniprot_5640_subset.xml`

### Create reference database (DuckDB file format)

- `poetry run python scripts/convert_uniprot_xml_to_duckdb.py tests/data/uniprot_5640_subset.xml.gz tests/data/uniprot_5640_subset_reference.db`

### Obtain per-peptide gene assignments

- `poetry run python scripts/merge_peptide_tables.py -o peptide_assignments.tsv -d ./tests/data/uniprot_5640_subset_reference.db ./tests/data/peptides_1.tsv ./tests/data/peptides_2.tsv`


## See also

- https://github.com/TomSmithCGAT/CamProt/tree/master/camprot/scripts
