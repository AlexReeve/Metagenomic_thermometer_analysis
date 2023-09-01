import pandas as pd
from itertools import combinations
import csv



def all_combinations(items_list):
    """
    creates a list of all combinations of amino acids 
    i.e. if the input was [a,b,c}, all_combinations would return [[a],[b],[c],[a,b],[b,c],[a,c],[a,b,c]].
    
    Arguments:
        items_list(list): List of all amino acids being investigated 
    
    Returns:
        result(list): list of all possible combinations of sets of amino acids.
    """
    result = []
    for i in range(1, len(items_list) + 1):
        for combo in combinations(items_list, i):
            result.append(combo)
    return result 


def listToString(s):
    """
    creates column names as a string from the input list.
    
    Arguments:
        s(list): list of amino acids being investigated
    
    Returns:
        str1(string): a string made up of all items from input list
    """
 
    # initialize an empty string
    str1 = ""
 
    # traverse in the string
    for ele in s:
        str1 += ele
 
    # return string
    return str1

def add_acids(data, combination_list):
    '''
    Finds the ten strongest postively and ten strongest negatively correlated sets of amino acids with a set of genomes. Searches through every possible comination of amino acid sets. Saves the best results as CSV files "best_positive_corrs.csv" and "best_negative_corrs.csv"
    Prints the current amino acid set bing tested, the current best positive and negative correlations, and the percentage completeness of the program.
    
    Arguments:
        data(dataframe): A dataframe containing the amino acid fraction of each genome. Must contain the folowing columns: ["OGT", "I_fraction", "V_fraction", "Y_fraction", "W_fraction", "R_fraction", "E_fraction", "L_fraction", "A_fraction", "N_fraction", "D_fraction", "C_fraction", "Q_fraction", "G_fraction", "H_fraction", "K_fraction", "M_fraction", "F_fraction", "P_fraction", "S_fraction", "T_fraction"]
        
        combination_list(list): a list of all combinations of amino acid sets.
    
    Returns:
        Nothing.
    
    '''
    # create our output dataframe
    best_corr= 0
    best_no = 0
    best_corr_dict = {}
    worst_corr = 0
    worst_no = 0
    worst_corr_dict = {}
    
    for combination in combination_list:
        # tell me how far through we are
        print("------------------------------------------------------")
        print(combination)
        index = combination_list.index(combination)
        print(f"{index +1} out of {len(combination_list)}")
        print(f"{((index+1)/(len(combination_list)))*100}% complete")
        print(f"best positive corr is {best_corr}")
        print(f"best positive corr number is{best_no}")
        print(f"best negative corr is {worst_corr}")
        print(f"best negative number is{worst_no}")
        #create column
        string = listToString(combination)
        
        add_dataframe = pd.DataFrame()
        add_dataframe["OGT"] = data["OGT"]
       
        ## add columns by acids
        
        add_dataframe[string]=0
        for acid in combination:
            add_dataframe[string] = add_dataframe[string] + data[acid]
        
        #corralate it
        corr = add_dataframe.corr()
        
        corr_number = corr.loc["OGT",string]
        print(f"corr number is {corr_number}")
        if corr_number >= best_no:
            best_no = corr_number
            best_corr = combination
        if corr_number <= worst_no:
            worst_no = corr_number
            worst_corr = combination
        ## making top ten correlations
        if index<10:
            best_corr_dict[combination] = corr_number
            worst_corr_dict[combination] = corr_number

        elif index>= 10:
            # find the maximum and minimum values in the respective dictionaries
            min_key = min(best_corr_dict, key=best_corr_dict.get)
            min_value = min(best_corr_dict.values())
            max_key = max(worst_corr_dict, key=worst_corr_dict.get)
            max_value = max(worst_corr_dict.values())
            ## see if corr value is larger or smaller then those values
            #if so, replace those values
            if corr_number >= min_value:
                del best_corr_dict[min_key]
                best_corr_dict[combination] = corr_number
            if corr_number <= max_value:
                del worst_corr_dict[max_key]
                worst_corr_dict[combination] = corr_number


    #Open a csv file 
    with open("best_positive_corrs.csv", "w", newline="") as fp:
   
        writer = csv.DictWriter(fp, fieldnames=best_corr_dict.keys())

        #Write the header 
        writer.writeheader()

        #Write the data 
        writer.writerow(best_corr_dict)
        print('Done writing best positive dictionary to a csv file')
    ###########
     # Open a file 
    with open("best_negative_corrs.csv", "w", newline="") as fp:
    
        writer = csv.DictWriter(fp, fieldnames=worst_corr_dict.keys())

        #Write the header
        writer.writeheader()

        #Write the data 
        writer.writerow(worst_corr_dict)
        print('Done writing best negative dictionary to a csv file')

            
        



def best_sequence_finder(data):
    '''
    Finds the ten strongest postively and ten strongest negatively correlated sets of amino acids with a set of genomes. Searches through every possible comination of amino acid sets. Saves the best results as CSV files (see add_acids function)
    
    Arguments:
        data(dataframe): A dataframe containing the amino acid fraction of each genome. Must contain the folowing columns: ["OGT", "I_fraction", "V_fraction", "Y_fraction", "W_fraction", "R_fraction", "E_fraction", "L_fraction", "A_fraction", "N_fraction", "D_fraction", "C_fraction", "Q_fraction", "G_fraction", "H_fraction", "K_fraction", "M_fraction", "F_fraction", "P_fraction", "S_fraction", "T_fraction"]
    
    Returns:
        Nothing.
    '''
    list_of_potential_amino_acids = ["I_fraction", "V_fraction", "Y_fraction", "W_fraction", "R_fraction", "E_fraction", "L_fraction", "A_fraction", "N_fraction", "D_fraction", "C_fraction", "Q_fraction", "G_fraction", "H_fraction", "K_fraction", "M_fraction", "F_fraction", "P_fraction", "S_fraction", "T_fraction"]
    all_acid_combinations = all_combinations(list_of_potential_amino_acids)
    add_acids(data, all_acid_combinations)
    


# get database from input as a csv file
# eg database.csv
input_database = sys.argv[1]
#import the genome list
data = pd.read_csv(f'{input_database}')
# clean up dataset
data =  data[["OGT", "I_fraction", "V_fraction", "Y_fraction", "W_fraction", "R_fraction", "E_fraction", "L_fraction", "A_fraction", "N_fraction", "D_fraction", "C_fraction", "Q_fraction", "G_fraction", "H_fraction", "K_fraction", "M_fraction", "F_fraction", "P_fraction", "S_fraction", "T_fraction"]]

best_sequence_finder(data)
