# commands/cnstates.py

import pandas as pd
import numpy as np
import os
from sklearn.metrics import auc, precision_recall_curve, roc_auc_score, average_precision_score, accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, cohen_kappa_score, brier_score_loss

import warnings


warnings.filterwarnings("ignore")

def add_cnstates_subparser(subparsers):
    """
    Add the cnstates subcommand to the CLI tool.

    This subcommand is used to calculate CN state detection metrics at the bin level.
    """
    parser = subparsers.add_parser(
        "cnstates",
        help="Calculate CN state detection metrics at the bin level."
    )
    parser.add_argument(
        "--tool-maternal-cna-files", 
        type=str, 
        nargs = "+",
        metavar="",
        help="Paths to the maternal CNA files generated by different tools. , separated by spaces"
    )
    parser.add_argument(
        "--tool-paternal-cna-files", 
        type=str, 
        nargs = "+",
        metavar="",
        help="Paths to the maternal CNA files generated by different tools. , separated by spaces"
    )
    parser.add_argument(
        "--tool-names", 
        type=str, 
        nargs = "+",
        metavar="",
        help="List of tool names corresponding to the CNA files."
    )
    parser.add_argument(
        "--profile-maternal-cna-files", 
        type=str, 
        metavar="",
        help="Path to the maternal ground truth CNA file."
    )
    parser.add_argument(
        "--profile-paternal-cna-files", 
        type=str, 
        metavar="",
        help="Path to the paternal ground truth CNA file."
    )
    parser.add_argument(
        "--output-dir", 
        type=str, 
        default="./output", 
        metavar="",
        help="Path to the file where the results will be saved, default: './output' "
    )
    return parser

def categorize_and_save_cnv(ground_truth_file, predict_file, save_dir,tool, condition,type = "maternal"):
    """
    Categorize and save CNVs based on a specific condition.

    Parameters:
        ground_truth_file (str): Path to the ground truth file.
        predict_file (str): Path to the prediction file.
        save_dir (str): Directory to save categorized results.     
        tool (str): Tool name.
        condition (str): Condition for categorization (e.g., '>=2').
        cna_type (str): Type of CNA ('maternal' or 'paternal').
    """
    ground_truth_df = pd.read_csv(ground_truth_file)
    predict_df = pd.read_csv(predict_file)
    
    def categorize_cn(value):
        if condition == ">=2":
            return int(value >= 2)
        
        if condition == ">=9":
            return int(value >= 9)
        

        if condition.startswith("="):
            try:
                threshold = float(condition[1:])  
                return int(value == threshold)
            except ValueError:
                raise ValueError(f"Invalid threshold in condition: {condition}")
        raise ValueError(f"Unsupported condition: {condition}")
    
    
    ground_truth_categorized = ground_truth_df.copy()
    predict_categorized = predict_df.copy()

    for col in ground_truth_df.columns[1:]:  
        ground_truth_categorized[col] = ground_truth_df[col].apply(categorize_cn)
        predict_categorized[col] = predict_df[col].apply(categorize_cn)

    os.makedirs(save_dir, exist_ok=True)
    ground_truth_categorized.to_csv(os.path.join(save_dir, f'ground_truth_{type}.csv'), index=False)
    predict_categorized.to_csv(os.path.join(save_dir, f'{tool}_predict_{type}.csv'), index=False)


def align_dataframes(ground_truth_df, predict_df):
    """
    Align two dataframes by their common columns and index.

    Parameters:
        ground_truth_df (DataFrame): Ground truth dataframe.
        predict_df (DataFrame): Predicted values dataframe.

    Returns:
        Tuple[DataFrame, DataFrame]: Aligned ground truth and predicted dataframes.
    """
    common_columns = ground_truth_df.columns.intersection(predict_df.columns)
    ground_truth_df = ground_truth_df[common_columns]
    predict_df = predict_df[common_columns]
    
    common_index = ground_truth_df.index.intersection(predict_df.index)
    ground_truth_df = ground_truth_df.loc[common_index]
    predict_df = predict_df.loc[common_index]
    
    return ground_truth_df, predict_df


def calculate_metrics(ground_truth_df, predict_df):
    """
    Calculate evaluation metrics for CNV detection.

    Parameters:
        ground_truth_df (DataFrame): Ground truth values.
        predict_df (DataFrame): Predicted values.

    Returns:
        dict: A dictionary containing evaluation metrics.
    """
    y_true = ground_truth_df.values.flatten()  
    y_pred = predict_df.values.flatten()  

    unique_y_true = set(y_true)
    unique_y_pred = set(y_pred)
    
    metrics = {}
    
    if len(set(y_true)) > 1 and len(set(y_pred)) > 1:
        metrics['AUROC'] = roc_auc_score(y_true, y_pred)
        precision, recall, _ = precision_recall_curve(y_true, y_pred)
        metrics['AUPRC'] = auc(recall, precision)
    else:
        metrics['AUROC'] = None
        metrics['AUPRC'] = None

    metrics['Accuracy'] = accuracy_score(y_true, y_pred)
    metrics['Precision'] = precision_score(y_true, y_pred, zero_division=0)
    metrics['Recall'] = recall_score(y_true, y_pred, zero_division=0)
    metrics['F1 Score'] = f1_score(y_true, y_pred, zero_division=0)
    metrics['Kappa'] = cohen_kappa_score(y_true, y_pred)
    metrics['Brier Score'] = brier_score_loss(y_true, y_pred)

    unique_labels = sorted(unique_y_true.union(unique_y_pred))

    # Check if y_true and y_pred contain more than one class before calculating confusion matrix
    if len(set(y_true)) > 1 and len(set(y_pred)) > 1:
        # Calculate confusion matrix
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred,labels=unique_labels).ravel()
        metrics['Specificity'] = tn / (tn + fp) if (tn + fp) != 0 else None
        metrics['PPV'] = metrics['Precision']  # PPV is the same as Precision in binary classification
        metrics['NPV'] = tn / (tn + fn) if (tn + fn) != 0 else None
    else:
        # If we cannot calculate confusion matrix due to lack of classes
        metrics['Specificity'] = None
        metrics['PPV'] = None
        metrics['NPV'] = None

    return metrics


def process_files_in_folder(folder_path,tool):
    results = {}
    for file_type in ['maternal', 'paternal']:
        ground_truth_file = os.path.join(folder_path, f'ground_truth_{file_type}.csv')
        predict_file = os.path.join(folder_path, f'{tool}_predict_{file_type}.csv')
        
        if os.path.exists(ground_truth_file) and os.path.exists(predict_file):
            ground_truth_df = pd.read_csv(ground_truth_file)
            predict_df = pd.read_csv(predict_file)

            ground_truth_df, predict_df = align_dataframes(ground_truth_df, predict_df)

            results[file_type] = {}
            
            for clone in predict_df.columns[1:].str.split('_').str[0].unique():
                ground_truth_subclone_df = ground_truth_df.loc[:, ground_truth_df.columns.str.split('_').str[0] == clone]
                predict_subclone_df = predict_df.loc[:, predict_df.columns.str.split('_').str[0] == clone]

                metrics = calculate_metrics(ground_truth_subclone_df, predict_subclone_df)

                results[file_type][clone] = metrics

    return results

def save_metrics_to_csv(results, folder_name, folder_path):
    for file_type,clone_metrics in results.items():
        all_metrics = []
        for clone, metrics in clone_metrics.items():
            metrics_with_info = {
                'Clone': clone
            }
            metrics_with_info.update(metrics)
            all_metrics.append(metrics_with_info)

        metrics_df = pd.DataFrame(all_metrics)
        save_path = os.path.join(folder_path, f'metrics_{folder_name}_{file_type}.csv')
        metrics_df.to_csv(save_path, index=False)



# python src/main.py cnstates 
# --tool-maternal-cna-files data/chisel_5x_maternal_cnv.csv 
# --tool-paternal-cna-files data/chisel_5x_paternal_cnv.csv 
# --tool-name chisel_5x 
# --profile-maternal-cna-files data/ground_truth_maternal_cnv.csv 
# --profile-paternal-cna-files data/ground_truth_paternal_cnv.csv

def run(args):
     # Validate input
    if (len(args.tool_maternal_cna_files) != len(args.tool_names) ) and (len(args.tool_maternal_cna_files) != len(args.tool_paternal_cna_files)):
        raise ValueError("The number of classification files must match the number of tool names.")
    
    profile_maternal_path = args.profile_maternal_cna_files
    profile_paternal_path = args.profile_paternal_cna_files
    
    conditions_and_dirs = [
        ('>=2', 'CN_Gain'),
        ('=1', 'CN_Neutral'),
        ('=0', 'CN_Loss'),
        ('=2', 'CN_equal_2'),
        ('=3', 'CN_equal_3'),
        ('=4', 'CN_equal_4'),
        ('=5', 'CN_equal_5'),
        ('=6', 'CN_equal_6'),
        ('=7', 'CN_equal_7'),
        ('=8', 'CN_equal_8'),
        ('>=9', 'CN_over_9'),
    ]
    
    os.makedirs(args.output_dir, exist_ok=True)

    all_metrics = []
    for tool_maternal_path,tool_paternal_path, tool_name in zip(args.tool_maternal_cna_files, args.tool_paternal_cna_files,args.tool_names):
        for condition, save_dir in conditions_and_dirs:
            categorize_and_save_cnv(
                ground_truth_file=profile_maternal_path,  
                predict_file=tool_maternal_path,  
                save_dir=f"{args.output_dir}/{tool_name}/{save_dir}",
                tool=tool_name,
                condition=condition,
                type='maternal'
            )

        for condition, save_dir in conditions_and_dirs:
            categorize_and_save_cnv(
                ground_truth_file= profile_paternal_path,  
                predict_file=tool_paternal_path,  
                save_dir=f"{args.output_dir}/{tool_name}/{save_dir}",
                tool=tool_name,
                condition=condition,
                type='paternal'
            )  

        for condition, save_dir in conditions_and_dirs:
            results = process_files_in_folder(f"{args.output_dir}/{tool_name}/{save_dir}",tool_name)
            for file_type, clone_metrics in results.items():
                for clone, metrics in clone_metrics.items():
                    metrics_with_info = {
                        'Type': save_dir,
                        'Haplotype': file_type,
                        'Clone': clone,
                        'tool':tool_name
                    }
                    metrics_with_info.update(metrics)  
                    
                    all_metrics.append(metrics_with_info)
    
    result_df = pd.DataFrame(all_metrics)
    output_file = os.path.join(args.output_dir, "cnstates_results.csv")
    result_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")