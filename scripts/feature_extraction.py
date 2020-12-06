import json
from sys import argv
import argparse
from tempfile import NamedTemporaryFile
from subprocess import run, DEVNULL
from Bio import AlignIO, SeqIO
from tqdm import tqdm


parser = argparse.ArgumentParser(description="Feature extraction from given fasta sequences.")
parser.add_argument('--reference_aln',
                    required=True,
                    nargs='?',
                    type=str,
                    metavar='<>.fasta',
                    help='reference alignment file')
parser.add_argument('--sequences',
                    required=True,
                    nargs='?',
                    type=str,
                    metavar='<>.fasta',
                    help='sample proteins to generate features for')
parser.add_argument('--reference_features',
                    required=True,
                    nargs='?',
                    type=str,
                    metavar='<>.json',
                    help='reference alignment feature columns')
parser.add_argument('--output',
                    required=True,
                    nargs='?',
                    type=str,
                    metavar='<>.json',
                    help='output json file with extracted features for each protein')

args = parser.parse_args(argv[1:])
reference_alignment_file = args.reference_aln
sample_sequences_file = args.sequences
reference_feature_columns_file = args.reference_features
output_file = args.output

reference_alignment = AlignIO.read(reference_alignment_file, 'fasta')
with open(reference_feature_columns_file, 'r') as infile:
    reference_features = json.load(infile)
sample_sequences = list(SeqIO.parse(sample_sequences_file, 'fasta'))


def create_sample_alignment(reference_file, sample_sequence):
    with NamedTemporaryFile('r+') as temp_seq:
        with NamedTemporaryFile('r+') as temp_aln:
            SeqIO.write(sample_sequence, temp_seq.name, 'fasta')
            run(f'muscle -profile -in1 {reference_file} -in2 {temp_seq.name} -out {temp_aln.name}', shell=True, stderr=DEVNULL)
            sample_alignment = AlignIO.read(temp_aln, 'fasta')
    return sample_alignment


def get_alignment_map_index(reference_alignment, sample_alignment):
    reference_length = reference_alignment.get_alignment_length()
    sample_length = sample_alignment.get_alignment_length()
    if reference_length == sample_length:
        alignment_map_index = {i: i
                               for i in range(reference_alignment.get_alignment_length())}
    else:
        alignment_map_index = dict()
        sample_pos = 0
        for reference_pos in range(reference_length):
            reference_column = reference_alignment[:, reference_pos]
            sample_column = sample_alignment[:-1, sample_pos]
            while sample_pos < sample_length and reference_column != sample_column:
                sample_pos += 1
                sample_column = sample_alignment[:-1, sample_pos]
            alignment_map_index[reference_pos] = sample_pos
    return alignment_map_index, all([i in alignment_map_index.keys()
                                     for i in range(reference_length)])


def generate_features(sample_alignment, alignment_map, feature_columns):
    sample_sequence = sample_alignment[-1]
    sample_feature_columns = dict()
    for name, columns in feature_columns.items():
        if isinstance(columns[0], int):
            sample_columns = [sample_sequence[alignment_map[column]]
                              for column in columns]
        elif isinstance(columns[0], list):
            sample_columns = [[sample_sequence[alignment_map[pair[0]]],
                               sample_sequence[alignment_map[pair[1]]]]
                              for pair in columns]
        else:
            raise ValueError(f'Not list and not int in {name} feature')
        sample_feature_columns[name] = sample_columns
    return sample_feature_columns


sequences_features = dict()
for sequence in tqdm(sample_sequences):
    sample_alignment = create_sample_alignment(reference_alignment_file,
                                               sequence)
    alignment_map, correct = get_alignment_map_index(reference_alignment,
                                                     sample_alignment)
    if not correct:
        print(f'{sequence.name}: some alignment columns are incorrect.')
    sample_features = generate_features(sample_alignment,
                                        alignment_map,
                                        reference_features)
    sequences_features[sequence.name] = sample_features


with open(output_file, 'w') as outfile:
    json.dump(sequences_features, outfile)
