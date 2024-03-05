# Import Libraries

import pandas as pd
import numpy as np
import seaborn as sns
from typing import Optional
from tqdm import tqdm
from time import sleep




def distance_matrix(string1: str,
                    string2: str,
                    show_heatmap: Optional[bool]=False):
    """Calculates the distance matrix between two strings

    Args:
        string1 (str): any character string
        string2 (str): any character string
        show_heatmap (Optional[bool], optional): Optionally creates a visual heatmap \ 
        of the distance matrix. Defaults to False.

    Returns:
        pd.DataFrame. The distance matrix between any two strings as a DataFrame. \ 
        Optionally can return a heatmap.
    """
    
    # Makes empty matrix
    matrix = np.empty((len(string1),len(string2)))
    matrix[:]=np.nan
    matrix_df = pd.DataFrame(matrix)

    # Compares each element in string1 and string2
    # Code puts a 0 where element is the same and 0 otherwise
    for i in enumerate(string1):
        for j in enumerate(string2):
            matrix_df.iloc[i[0],j[0]]=int(i[1]==j[1])
        
    matrix_df = matrix_df.astype(int) #Makes sure matrix is an integer
    
    # Makes the rows and columns of the dataframe string1 and string2, respectively
    matrix_df.columns = list(string2)
    matrix_df.index = list(string1)
    
    # If heatmap is set to True, draws a heatmap according to the heatmap function
    if show_heatmap:
        distance_matrix_heatmap(matrix_df)
    
    return matrix_df



def distance_matrix_heatmap(distance_matrix_df: pd.DataFrame):
    """Visualization of distance matrix as a heatmap

    Args:
        distance_matrix_df (pd.DataFrame): A dataframe from the output of distance_matrix function.
    """
  
    # Create the heatmap
    g = sns.heatmap(
        distance_matrix_df,
        annot = True,
        square=True,
        cbar_kws={'fraction' : 0.01},
        cmap='OrRd',
        linewidth=1)
    
    # Customize the heatmap appearance 
    g.tick_params(left=False, bottom=False)
    g.set_xticklabels(g.get_xticklabels(), rotation=0, horizontalalignment='right')
    g.set_yticklabels(g.get_yticklabels(), rotation=0, horizontalalignment='right')
    None # prevent the list of label objects showing up annoyingly in the output
    return



def percent_similarity(string1: str,
                 string2: str,
                 view_diags: Optional[bool]=False,
                 diag_type: Optional[str]='all'):
    """Finds the axis of greatest similarity between two strings based on the distance matrix 
    and returns the fractional similarity along the main diagonal or finds the maximum 
    similarity along any diagonal.

    Args:
        string1 (str): Any character string
        string2 (str): Any character string
        view_diags (Optional[bool], optional): If True, returns all possible string alignments. Defaults to False.
        diag_type (Optional[str], optional): Specifies the type of diagonal to consider. 
            Choose from:
            - 'main': Considers only the main diagonal of the distance matrix.
            - 'all': Computes similarity across all possible alignments. Defaults to 'all'.

    Returns:
        Tuple[float, Optional[list]]: The fractional similarity along the main diagonal or the maximum similarity 
            along any diagonal, and optionally, a list of diagonals if `view_diags` is True.
    """

    # Create a distance matrix
    # Put the distance matrix into a square padded by NaNs
    A = distance_matrix(string1,string2,show_heatmap=False)
    dim = max(A.shape)
    matrix1 = np.empty((dim,dim))
    matrix1[:]=np.nan
    matrix1[0:A.shape[0],0:A.shape[1]] = A
    
    #Take the sum of the diagonals
    # Main returns percent similarity along the main diagonal of the square
    # The main diagonal is equivalent to aligning the stings left-to-right
    if diag_type=='main':
        diag_sum = np.nansum(np.diagonal(matrix1,offset=0))
        diag_len = len(np.diagonal(matrix1,offset=0))
        sim_max = 100*(np.array(diag_sum)/np.array(diag_len))
    # All returns percent similarity along the max diagonal
    # The max of all diagonals is equivalent to the max alignment of the strings
    elif diag_type=='all':
        diag_sum = [np.nansum(np.diagonal(matrix1,offset=i)) for i in range(-len(matrix1)+1,len(matrix1),1)]
        #diag_len = [len(np.diagonal(matrix1,offset=i)) for i in range(-len(matrix1)+1,len(matrix1),1)]
        diag_len = len(np.diagonal(matrix1,offset=0))
        sim_max = max(100*(np.array(diag_sum)/np.array(diag_len)))
    
    # When true, all diagonals of the square are returned 
    if view_diags:
        diags_view = [np.diagonal(matrix1,offset=i) for i in range(-len(matrix1)+1,len(matrix1),1)]
    else:
        diags_view = None
        
    return sim_max, diags_view


def run_all_CDR3s(data: pd.DataFrame,
                  Num_top_CDR3: Optional[int]=10,
                  Observations_col: Optional[str]='umi_number',
                  CDR3_col_name: Optional[str]='pep'):
    """Runs the similarity functions on a list of TCR CDR3 strings and returns a similarity matrix of the top 10 CDR3s with all CDR3s in the list

    Args:
        data (pd.DataFrame): DataFrame of with a column of strings and a column of frequency of observation.
        Num_top_CDR3 (Optional[int], optional): Number of top CDR3s to compare to the whole. Higher number takes more time. Defaults to 10.
        Observations_col (Optional[str], optional): Name of column with string frequencies. Defaults to 'umi_number'.
        CDR3_col_name (Optional[str], optional): Name of column that contains CDR3 strings. Defaults to 'pep'.

    Returns:
        pd.DataFrame: A matrix of the top most frequent CDR3s and their similarity to the entire list of CDR3s in the data.
    """
    #Prepare the data
    data_sorted = data.sort_values(by=Observations_col, ascending=False)
    data_sorted_top = data_sorted.iloc[:Num_top_CDR3]
    
    dim_all = len(data_sorted)
    dim_top = len(data_sorted_top)
    CDR3_list = list(data_sorted[CDR3_col_name])
    CDR3_list_umi = data_sorted[Observations_col]
    umi_sum = sum(CDR3_list_umi)
    CDR3_list_umi.index = CDR3_list
    umi_percent = 100*(CDR3_list_umi/umi_sum)
    CDR3_list_top = list(data_sorted_top[CDR3_col_name])
    
    matrix2 = np.empty((dim_all, dim_top))
    matrix2[:]=np.nan
    matrix2_df = pd.DataFrame(matrix2)
    matrix2_df.columns = CDR3_list_top
    matrix2_df.index = CDR3_list
    
    print("matrix is made")

    # Create a progress bar for loop
    for i in tqdm(range(dim_top), desc="Processing", unit="item"):
        
        if i %10 == 0:
            print(f'{i}_loop')
        
        for j in range(dim_all):
            matrix2_df.iloc[j,i],_ = percent_similarity(CDR3_list_top[i], CDR3_list[j])
    
    # Tack information columns onto the dataframe
    matrix2_df.insert(0,'umi_number',CDR3_list_umi)
    matrix2_df.insert(0,'umi_percent',umi_percent)

    return matrix2_df


def filtersimilarity(cdr3map: pd.DataFrame,
                     cdr3: [str],
                     umi_cols: [list] =['umi_percent','umi_number'],
                     truncate: Optional[bool]=True,
                     sim_threshold: [int]=70):
    """Filters CDR3s from the cdr3map for similarity, allowing the visibility of the top most similar CDR3s

    Args:
        cdr3map (pd.DataFrame): cdr3 similarity map given by run_all_CDR3s
        cdr3 (str]): _description_
        umi_cols (list], optional): _description_. Defaults to ['umi_percent','umi_number'].
        truncate (Optional[bool], optional): _description_. Defaults to True.
        sim_threshold (int], optional): _description_. Defaults to 70.

    Returns:
        _type_: _description_
    """
    
    cdr3map_final = cdr3map
    umi_cols.append(cdr3)
    matrix = cdr3map_final.loc[:,umi_cols]
    sorted_df = matrix.sort_values(by=cdr3,ascending=False)
    filtered_df = sorted_df[sorted_df[cdr3] > 70]
    
    if truncate:
        return filtered_df
    else:
        return sorted_df