# commands/mirrorsubclone.py

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score as adjustedRandIndex
from sklearn.metrics import adjusted_mutual_info_score as AMI
import os


def add_mirrorsubclone_subparser(subparsers):
    """
    Add the mirrorsubclone subcommand to the CLI tool.

    This subcommand is used to calculate the RMSE,ACC of Mirror-subclone CNA.
    """
    parser = subparsers.add_parser(
        "mirrorsubclone",
        help="Calculate the RMSE,ACC of Mirror-subclone CNA."
    )
    parser.add_argument(
        "--tool-cna-files", 
        type=str, 
        nargs = "+",
        required=True, 
        metavar="",
        help="Paths to the CNA files generated by different tools."
    )
    parser.add_argument(
        "--tool-names", 
        type=str, 
        nargs = "+",
        required=True, 
        metavar="",
        help="List of tool names corresponding to the classification files"
    )
    parser.add_argument(
        "--changes-file", 
        type=str, 
        required=True, 
        metavar="",
        help="Path to the changes file that contains information about changes in CN status."
    )
    parser.add_argument(
        "--output-dir", 
        type=str, 
        default="./output", 
        metavar="",
        help="Path to the file where the results will be saved, default: './output' "
    )
    return parser

def process_mirror_change(change_file_path):

    change_ground_truth = pd.read_csv(change_file_path)

    change_ground_truth['Segment'] = (
        change_ground_truth['Chromosome'] + ":" +
        change_ground_truth['Start'].astype(str) + "-" +
        change_ground_truth['End'].astype(str)
    )

    return change_ground_truth


def calculate_predictions(change_data, predict_data):
    predict_sorted = predict_data.merge(change_data, left_on='region', right_on='Segment', how='inner')

    if predict_sorted.empty:
        print("Predict sorted is empty.")
        return pd.DataFrame(columns= ['Segment','Clone1', 'Clone2', 'Clone1_CNA','Clone2_CNA', 'Clone1_predict_CNA', 'Clone2_predict_CNA'])

    def calculate_mode_for_columns(row, columns_prefix, df):
        columns = [col for col in df.columns if col.startswith(row[columns_prefix])]
        if columns:
            values = df.loc[row.name, columns]
            mode_value = values.mode()
            return mode_value.iloc[0] if not mode_value.empty else None
        else:
            return None

    predict_sorted['Clone1_predict_CNA'] = predict_sorted.apply(
        calculate_mode_for_columns, axis=1, columns_prefix='Clone1', df=predict_sorted
    )
    
    predict_sorted['Clone2_predict_CNA'] = predict_sorted.apply(
        calculate_mode_for_columns, axis=1, columns_prefix='Clone2', df=predict_sorted
    )

    columns_to_keep = ['Segment','Clone1', 'Clone2', 'Clone1_CNA','Clone2_CNA', 'Clone1_predict_CNA', 'Clone2_predict_CNA']
    predict_sorted_cleaned = predict_sorted[columns_to_keep]
    
    return predict_sorted_cleaned


##python src/main.py mirrorsubclone --tool-cna-files data/chisel_5x_combined_cnv.csv --tool-names chisel_5x --changes-file profile/mirrored_subclonal_cnas.csv

def run(args):
    # Validate input
    if len(args.tool_cna_files) != len(args.tool_names):
        raise ValueError("The number of classification files must match the number of tool names.")
    
    os.makedirs(args.output_dir, exist_ok=True)

    change_data = process_mirror_change(args.changes_file)

    results = []
    for tool_path, tool_name in zip(args.tool_cna_files, args.tool_names):

        predict_combined = pd.read_csv(tool_path)

        combined_result = calculate_predictions(change_data, predict_combined)

        combined_result['Clone1_maternal_CNA'] = combined_result['Clone1_CNA'].str.split('|').str[0]
        combined_result['Clone1_paternal_CNA'] = combined_result['Clone1_CNA'].str.split('|').str[1]
        combined_result['Clone2_maternal_CNA'] = combined_result['Clone2_CNA'].str.split('|').str[0]
        combined_result['Clone2_paternal_CNA'] = combined_result['Clone2_CNA'].str.split('|').str[1]

        combined_result['Clone1_predict_maternal_CNA'] = combined_result['Clone1_predict_CNA'].str.split('|').str[0]
        combined_result['Clone1_predict_paternal_CNA'] = combined_result['Clone1_predict_CNA'].str.split('|').str[1]
        combined_result['Clone2_predict_maternal_CNA'] = combined_result['Clone2_predict_CNA'].str.split('|').str[0]
        combined_result['Clone2_predict_paternal_CNA'] = combined_result['Clone2_predict_CNA'].str.split('|').str[1]

        combined_result['acc_result'] = (
            (combined_result['Clone1_predict_maternal_CNA'].astype(str) == combined_result['Clone1_maternal_CNA'].astype(str) )  &
            (combined_result['Clone1_predict_paternal_CNA'].astype(str) ==combined_result['Clone1_paternal_CNA'].astype(str) ) &
            (combined_result['Clone2_predict_maternal_CNA'].astype(str) == combined_result['Clone2_maternal_CNA'].astype(str) )  &
            (combined_result['Clone2_predict_paternal_CNA'].astype(str) ==combined_result['Clone2_paternal_CNA'].astype(str) ) 
        )

        combined_result['Clone1_maternal_error'] = (
            (combined_result['Clone1_predict_maternal_CNA'].astype(float) - combined_result['Clone1_maternal_CNA'].astype(float))**2
        )

        combined_result['Clone1_paternal_error'] = (
            (combined_result['Clone1_predict_paternal_CNA'].astype(float) - combined_result['Clone1_paternal_CNA'].astype(float))**2
        )

        combined_result['Clone2_maternal_error'] = (
            (combined_result['Clone2_predict_maternal_CNA'].astype(float) - combined_result['Clone2_maternal_CNA'].astype(float))**2
        )

        combined_result['Clone2_paternal_error'] = (
            (combined_result['Clone2_predict_paternal_CNA'].astype(float) - combined_result['Clone2_paternal_CNA'].astype(float))**2
        )

        combined_result['RMSE_result'] = (
            np.sqrt(combined_result['Clone1_maternal_error']) +
            np.sqrt(combined_result['Clone1_paternal_error']) +
            np.sqrt(combined_result['Clone2_maternal_error']) +
            np.sqrt(combined_result['Clone2_paternal_error'])
        )

        rmse = combined_result['RMSE_result'].mean()

        acc = combined_result['acc_result'].mean()

        results.append({
            "Tool": tool_name,
            "RMSE": rmse,
            "ACC": acc
        })
    
    result_df = pd.DataFrame(results)
    output_file = os.path.join(args.output_dir, "mirror_subclone_result.csv")
    result_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
