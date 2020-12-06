import argparse
from tempfile import NamedTemporaryFile
import sys
from subprocess import run
import os
from pathlib import Path


src_path = Path(os.path.abspath(os.path.dirname(__file__)))

REFERENCE_ALN = src_path / "data" / "reference_alignment.fasta"
REFERENCE_FEATURES = src_path / "data" / "reference_feature_columns.json"
HMM = src_path / "data" / "DNA_pol_B.hmm"
MODEL = src_path / "data" / "catboost_model"
DEFAULT_THRES = 5

parser = argparse.ArgumentParser(description="",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--query',
                    type=str,
                    nargs='?',
                    required=True,
                    metavar='<>.fasta',
                    help='query protein sequences in fasta format')
parser.add_argument('--threshold',
                    type=float,
                    nargs='?',
                    default=DEFAULT_THRES,
                    metavar='x',
                    help='HMM e-value threshold for prefiltering')
parser.add_argument('--output',
                    type=str,
                    nargs='?',
                    required=True,
                    metavar='<>.csv',
                    help='output csv file with predictions ranked by proba')

args = parser.parse_args(sys.argv[1:])
query_file = args.query
threshold = args.threshold
output_file = args.output


filter_hmm_path = src_path / "scripts" / "filter_hmm.py"
feature_extraction_path = src_path / "scripts" / "feature_extraction.py"
catboost_predict_path = src_path / "scripts" / "catboost_predict.py"

with NamedTemporaryFile() as filtered, NamedTemporaryFile() as features:
    print(f'Filtering proteins with HMM at E-value <= {threshold} threshold...', file=sys.stderr)
    run(f"python {filter_hmm_path} --query {query_file} --out {filtered.name} --hmm {HMM} --thres {threshold}", shell=True)
    print(f'Generating features from filtered proteins...', file=sys.stderr)
    run(f"python {feature_extraction_path} --reference_aln {REFERENCE_ALN} --sequences {filtered.name} --reference_features {REFERENCE_FEATURES} --output {features.name}", shell=True)
    print(f'Running predictions for filtered proteins...', file=sys.stderr)
    run(f"python {catboost_predict_path} --feature_file {features.name} --model_file {MODEL} --output {output_file}", shell=True)
    print('Finished.', file=sys.stderr)