
from pathlib import Path

import numpy as np
import pandas as pd


DATA_DIR = Path(__file__).absolute().parent

CENTROMERE_B37_PATH = DATA_DIR.joinpath('centromere_b37.tsv')
CENTROMERE_B38_PATH = DATA_DIR.joinpath('centromere_b38.tsv')
REFFLAT_B37_PATH = DATA_DIR.joinpath('refFlat_b37.tsv')
REFFLAT_B38_PATH = DATA_DIR.joinpath('refFlat_b38.tsv')
<<<<<<< HEAD
RECOMB_B37_PATH = DATA_DIR.joinpath('recomb_b37.tsv')
RECOMB_B38_PATH = DATA_DIR.joinpath('recomb_b38.tsv')
=======
RECOMB_B37_PATH = DATA_DIR.joinpath('recomb_rate_b37.tsv')
RECOMB_B38_PATH = DATA_DIR.joinpath('recomb_rate_b38.tsv')
>>>>>>> master



CENTROMERE_B37 = pd.read_csv(
            CENTROMERE_B37_PATH,
            sep = '\t',
            header = 0,
            dtype = {
                'chrom': np.int64,
                'start': np.int64,
                'end': np.int64
            }
)

CENTROMERE_B38 = pd.read_csv(
            CENTROMERE_B37_PATH,
            sep = '\t',
            header = 0,
            dtype = {
                'chrom': np.int64,
                'start': np.int64,
                'end': np.int64
            }
)