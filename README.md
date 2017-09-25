# ConfPipeline

Scripts to run conformer generation pipeline

Example usage:

- source env.sh
- export TURBODIR=/path/to/turbomole
- source $TURBODIR/Config_turbo_env

```
import yaml
from pipeline import run_though_pipeline as conf_gen
from rdkit import Chem

config = yaml.load(open('pipeline.yml').read())
mol = Chem.MolFromSmiles('CCCCCC')
mol = conf_gen(mol, config)

print(mol.politzer[0])
```
