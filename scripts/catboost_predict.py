# this module takes as input file json from feature_extraction.py
# output files contains a ranked list of predictied phi29 activity
# with scores and labels where 0 stands for lack of activity, 1 stands for existing activity.
import pandas as pd
import json
import catboost as cat
from sys import argv
import argparse


parser = argparse.ArgumentParser(description="Predict phi29 activity based on generated features.")
parser.add_argument('--feature_file',
                    type=str,
                    required=True,
                    nargs='?',
                    metavar='<>.json',
                    help='input json file from feature_extraction.py')
parser.add_argument('--model_file',
                    type=str,
                    required=True,
                    nargs='?',
                    metavar="MODEL",
                    help='catboost model')
parser.add_argument('--output',
                    type=str,
                    required=True,
                    nargs='?',
                    metavar="<>.csv",
                    help='output csv file')

args = parser.parse_args(argv[1:])
feature_file = args.feature_file
model_file = args.model_file
output_file = args.output


def read_json(file_name):
    with open(file_name, "r") as read_file:
        data = json.load(read_file)
    return(pd.DataFrame(data).T)


def make_binary_disulf(ao_lists):
    result = ["1" if pair == ["C", "C"] else "0" for pair in ao_lists]
    return(result)


# read df and selected necessary features
df = read_json(feature_file)
df = df[['catalytic',
         'dntp_binding',
         'primer_binding',
         'exonuclease',
         'Mg',
         'replication_activity',
         'conservative',
         'ssbonds',
         'pockets']]
df.ssbonds = df.ssbonds.apply(make_binary_disulf)
new_colnames = []
counter = 0
for i, j in zip(df.columns[:], df.iloc[0, :]):
    for k in range(len(j)):
        new_colnames.append(str(counter) + '_' + i + '_' + str(k))
        counter += 1

df_label = pd.DataFrame()
for i in df.columns:
    df_label = pd.concat([df_label,
                          pd.DataFrame(df[i].tolist(),
                                       index=df.index)],
                         axis=1)

df_label.columns = new_colnames

selected = ['5_dntp_binding_2',
            '101_pockets_2',
            '14_primer_binding_4',
            '43_replication_activity_3',
            '15_primer_binding_5',
            '8_dntp_binding_5',
            '90_conservative_29',
            '66_conservative_5',
            '70_conservative_9',
            '39_Mg_5',
            '27_exonuclease_9',
            '93_conservative_32',
            '67_conservative_6',
            '28_exonuclease_10',
            '37_Mg_3']
df_label = df_label[selected]

clf_cat = cat.CatBoostClassifier()  # parameters are not required.
clf_cat.load_model(model_file)
y_pred = clf_cat.predict(df_label)
y_pred_proba = clf_cat.predict_proba(df_label)[:, -1]

labelsdf = pd.DataFrame()
labelsdf['name'] = df_label.index
labelsdf['proba'] = y_pred_proba
labelsdf = labelsdf.sort_values(by='proba')
labelsdf = labelsdf.sort_values(by='proba', ascending=False)
labelsdf.to_csv(output_file, sep='\t', index=False)
