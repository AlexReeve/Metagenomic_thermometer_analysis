import pandas as pd
import numpy as np
import time 
import sys

        
def find_taxids_in_both(data, genomes_data):
    '''
    Creates a list of all unique taxonomic identification numbers found in both given datasets  
    
    Arguments:
        data(dataframe): The dataset comtaining the OGTs of microorganisms. "Growth temperatures for 21,498 microorganisms" (Engqvist, 2018)
        
        genomes_data(dataframe): The dataset containing the FTP path to the genome.
        
    Returns:
        unique_taxids_in_both_arrays(array): An array containing all unique taxinomic identification numbers that are present in both datasets.

    '''
    ##find unique taxids in genomic dataframe
    unique_taxid_ids = genomes_data["species_taxid"].unique()
    #set up list
    taxids_in_both =[]
    #find if taxid is in both dataframes
    for ogt_taxid in data["taxid"]:
        for genomes_taxid in unique_taxid_ids :
            # if in both dataframes, add to list
            if genomes_taxid == ogt_taxid:
                taxids_in_both.append(ogt_taxid)
                
    #make it an array, and make sure the array contains only one of each taxanomic identification number
    taxids_in_both_array = np.array(taxids_in_both)
    unique_taxids_in_both_arrays = list(dict.fromkeys(taxids_in_both))
    unique_taxids_in_both_arrays = np.array(unique_taxids_in_both_arrays)
    return(unique_taxids_in_both_arrays)

def make_new_dataframe(genomes_data, unique_taxids_in_both_arrays):
    '''
    Creates a dataframe of all genomes with known optimal growth temperatures.
    
    Arguments: 
        genomes_data(dataframe): The NCBI summary list of available genomes.
        
        unique_taxids_in_both_arrays(array): An array of taxonomic identification numbers. This is an array of all taxids present in the genome dataframe and the known OGT dataframe.
        
    Returns:
        working_data(dataframe): A dataframe of all genomes available with known optimal growth temperatures (note- dataframe does not contain known optimal growth temperatures)
    '''
    #set up dataframe
    genomes_with_known_ogt_dataframe = pd.DataFrame()
    # go through taxids, and if genome has taxid, put genome in output dataframe
    for taxid_list_item in unique_taxids_in_both_arrays:
        for index, row in genomes_data.iterrows():
            if taxid_list_item == row["species_taxid"]:
                genomes_with_known_ogt_dataframe = pd.concat([genomes_with_known_ogt_dataframe, row], axis = 1)
    # format dataframe
    working_data = genomes_with_known_ogt_dataframe.T
    # tell me how many species we have
    print(f"number of species in data is {working_data['species_taxid'].nunique()}")
    return(working_data)

def marry_main(genomes_data, ogt_data):
    '''
    For each genome, takes a given optimal growth temperature (found using marry_find) and puts it in the column "OGT"
    
    Arguments:
        genomes_data(dataframe): The dataframe of genomes with known optimal growth temperatures.
        
        ogt_data(dataframe): The dataset comtaining the OGTs of microorganisms. "Growth temperatures for 21,498 microorganisms" (Engqvist, 2018)
        
    Returns:
        output(dataframe): A dataframe of genomes with their optimal growth temperatures.
    '''
    output = genomes_data
    #for each genome
    for index1, genomes_row in genomes_data.iterrows():
        #get the OGT
        found_ogt = marry_find(genomes_row, ogt_data)
        #put ogt in the right column
        output.loc[index1, "OGT"] = found_ogt
        
    return output

def marry_find(genomes_row, ogt_data):
    '''
    For a given genome, finds the optimal growth temperature from the dataset comtaining the OGTs of microorganisms. "Growth temperatures for 21,498 microorganisms" (Engqvist, 2018). Identifies OGT using taxonomic identification number.
    
    Arguments:
        genomes_row(dataframe): the row of the genome dataframe. 
        
        ogt_data(dataframe): a dataset containing the ogts of microorganisms by taxonomic identification number
        
    Returns:
        output(int): The optimal growth temperature of the input genome
    '''
    #look through ogt data for relevent taxid
    for index2, ogt_item in ogt_data.iterrows():
        if genomes_row["species_taxid"] == ogt_item["taxid"]:
            #get the ogt
            output = ogt_item["temperature"]
            return output
        

def make_dataset(data, genomes_data):
    '''
    Creates a dataframe of available genomes with known optimal growth temperatures from two dataframes, one of available genomes, and one of known OGTs.
    
    Arguments:
        data(dataframe): The dataframe of microorganisms with a known optimal growth temperature. "Growth temperatures for 21,498 microorganisms" (Engqvist, 2018)
        
        genomes_data(dataframe): The list of whole genomes available on NCBI
    '''
    # find taxids in both
    unique_taxids_in_both_arrays = find_taxids_in_both(data, genomes_data)
    # shorten the genome dataset to only include genomes we have the OGT of
    working_data = make_new_dataframe(genomes_data, unique_taxids_in_both_arrays)
    #Find the OGTs of the shortened genome list
    genome_ogt = marry_main(working_data, data)
    
    return(genome_ogt)

        
temp_data = pd.read_csv('temperature_data.csv')

input_name = sys.argv[1]

output_name = 'marry_find_output_dataset.csv'
    

genomes_data = pd.read_csv(f'{input_name}', encoding='latin-1')

output = make_dataset(temp_data, genomes_data)

output.to_csv(output_name)