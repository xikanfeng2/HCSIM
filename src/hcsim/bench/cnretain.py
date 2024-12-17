# commands/cnretain.py

import ast
import pandas as pd
from sklearn.metrics import adjusted_rand_score as adjustedRandIndex
from sklearn.metrics import adjusted_mutual_info_score as AMI
import os
import numpy as np
from Bio import Phylo

def add_cnretain_subparser(subparsers):
    """
    Add the cnretain subcommand to the CLI tool.

    This subcommand is used to calculate Evolutionary CN Stability.
    """
    parser = subparsers.add_parser(
        "cnretain",
        help="Calculate the ACC of Evolutionary CN Stability."
    )
    parser.add_argument(
        "--tool-maternal-cna-files", 
        type=str, 
        nargs = "+",
        required=True, 
        metavar="",
        help="Paths to the maternal CNA files generated by different tools, separated by spaces"
    )
    parser.add_argument(
        "--tool-paternal-cna-files", 
        type=str, 
        nargs = "+",
        required=True, 
        metavar="",
        help="Paths to the paternal CNA files generated by different tools, separated by spaces"
    )
    parser.add_argument(
        "--tool-names", 
        type=str, 
        nargs = "+",
        required=True, 
        metavar="",
        help="List of tool names corresponding to the maternal and paternal CNA files provided."
    )
    parser.add_argument(
        "--changes-file", 
        type=str, 
        required=True, 
        metavar="",
        help="Path to the changes file that contains information about changes in CN status."
    )
    parser.add_argument(
        "--tree-file", 
        type=str, 
        required=True, 
        metavar="",
        help="Path to the tree file that provides the hierarchical clonal structure."
    )
    parser.add_argument(
        "--output-dir", 
        type=str, 
        default="./output", 
        metavar="",
        help="Path to the file where the results will be saved, default: './output' "
    )
    return parser


def is_descendant_or_self(tree, ancestor_name, descendant_name):
    """
    Check if a node is a descendant of or identical to an ancestor node in the tree.

    Parameters:
    tree: Tree object to search.
    ancestor_name (str): Name of the ancestor node.
    descendant_name (str): Name of the descendant node.

    Returns:
    bool: True if the node is a descendant of or identical to the ancestor, otherwise False.
    """
    ancestor = tree.find_any(name=ancestor_name)
    descendant = tree.find_any(name=descendant_name)
    
    if ancestor is None or descendant is None:
        return False  
    if ancestor == descendant :
        return True

    for clade in ancestor.find_clades():
        if clade == descendant:
            return True
            
    return False

def get_descendant(tree,ancestor_name):
    """
    Retrieve the names of all descendant nodes of a given ancestor in the tree.

    Parameters:
    tree: Tree object to search.
    ancestor_name (str): Name of the ancestor node.

    Returns:
    list: List of names of descendant nodes.
    """
    ancestor = tree.find_any(name=ancestor_name)
    descendant = ancestor.find_clades()
            
    descendant_names = [clade.name for clade in descendant if clade.name is not None]
    
    return descendant_names

def add_check_list(tree,change_data):
    """
    Add a 'check_child_list' column to the change data for nodes requiring verification.

    Parameters:
    tree: Tree object containing hierarchical relationships.
    change_data (DataFrame): DataFrame with change data.

    Returns:
    DataFrame: Updated DataFrame with 'check_child_list' column.
    """

    if 'check_child_list' not in change_data.columns:
        change_data['check_child_list'] = None  
    change_data = change_data.astype({'check_child_list': 'object'})  
    for index,row in change_data.iterrows():
        start = row['Child'] 
        segment = row['Segment'] 
        haplotype = row['Haplotype']
        descendant = get_descendant(tree,start)
        sub_clone = change_data[(change_data['Segment']==segment) & (change_data['Haplotype'] == haplotype)]['Child'].unique()

        del_clone_list = []
        for clone in sub_clone:
            if (clone == start):
                continue
            elif is_descendant_or_self(tree, start,clone) :
                sub_descendant = list(get_descendant(tree, clone))
                del_clone_list += sub_descendant


        update_descendant = list(set(descendant) - set(del_clone_list))
        change_data.at[index, 'check_child_list'] = update_descendant

    return change_data


def calculate_mode_for_columns(row, columns_prefix, df):
    """
    Calculate the mode for a subset of columns matching a specific prefix.

    Parameters:
    row (Series): A row from the DataFrame.
    columns_prefix (str): Prefix to filter relevant columns.
    df (DataFrame): DataFrame containing the data.

    Returns:
    Any: The most frequent value (mode) or None if no values exist.
    """
    columns = [col for col in df.columns if col.startswith(f"{columns_prefix}_")]
    if columns:
        values = row[columns]
        mode_value = values.mode()
        return mode_value.iloc[0] if not mode_value.empty else None
    else:
        return None

def process_haplotype_result(row, predict_data, haplotype, change_data, index):
    """
    Process a single row of change data for a specific haplotype.

    Parameters:
    row (Series): Row of the change_data DataFrame.
    prediction_data (DataFrame): DataFrame containing prediction results.
    haplotype (str): Haplotype being processed ('maternal' or 'paternal').
    change_data (DataFrame): DataFrame containing change data.
    index (int): Index of the current row in change_data.
    """
    change_result = int(row['Change'].split('->')[1])
    region = row['Segment']

    tool_row = predict_data[predict_data['region'] == region]

    if tool_row.empty:
        return

    tool_row = tool_row.iloc[0]
    change_data.loc[index, 'result'] = 1
    clone_list = row['check_child_list']
    if isinstance(clone_list, str):
        clone_list = ast.literal_eval(clone_list)

    for clone in clone_list:
        clone_result = calculate_mode_for_columns(tool_row,clone,predict_data)
        if clone_result != change_result:
            change_data.loc[index, 'result'] = 0
    

def process_change_data(change_data, maternal_result, paternal_result):
    for index, row in change_data.iterrows():
        if row['Haplotype'] == 'maternal':
            process_haplotype_result(row, maternal_result, 'maternal', change_data, index)
        elif row['Haplotype'] == 'paternal':
            process_haplotype_result(row, paternal_result, 'paternal', change_data, index)


# python src/main.py cnretain 
# --tool-maternal-cna-files data/chisel_5x_maternal_cnv.csv 
# --tool-paternal-cna-files data/chisel_5x_paternal_cnv.csv 
# --tool-names chisel_5x --changes-file profile/changes.csv 
# --tree-file profile/tree.newick

def run(args):
    # Validate input
    if (len(args.tool_maternal_cna_files) != len(args.tool_names) ) and (len(args.tool_maternal_cna_files) != len(args.tool_paternal_cna_files)):
        raise ValueError("The number of classification files must match the number of tool names.")
    
    
    os.makedirs(args.output_dir, exist_ok=True)

    tree = Phylo.read(args.tree_file, "newick")
    change_data = pd.read_csv(args.changes_file)
    change_data = add_check_list(tree,change_data)
    print(change_data.head(),"----")

    results = []
    for tool_maternal_path,tool_paternal_path, tool_name in zip(args.tool_maternal_cna_files, args.tool_paternal_cna_files,args.tool_names):

        predict_maternal = pd.read_csv(tool_maternal_path)
        predict_paternal = pd.read_csv(tool_paternal_path)

        process_change_data(change_data, predict_maternal, predict_paternal)

        for type in change_data['Type'].unique().tolist():

            acc = change_data[change_data['Type']==type]['result'].mean()

            results.append({
                "Tool": tool_name,
                "Type": type,
                "ACC": acc
            })
    
    result_df = pd.DataFrame(results)
    output_file = os.path.join(args.output_dir, "evolution_cn_stability_acc.csv")
    result_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
