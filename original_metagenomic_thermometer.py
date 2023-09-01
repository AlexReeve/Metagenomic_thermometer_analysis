import pandas as pd
from Bio import SeqIO
import wget
import gzip
import shutil
import os
from ftplib import FTP
import time
import sys

def thermometer(genome):
    '''
    Calculates the predicted optimal growth temperature of any given genome.
    Counts the IVYRWEL protein fraction and uses that to predict OGT.
    Calculation is based upon Zeldovich, Berezovsky, and Shakhnovich (2007).
    
    Arguments:
        genome(protein fasta file): The genome used to predict OGT. Must be a protein fasta file. 
    
    Returns:
        predicted_ogt_2dp(int):The predicted optimal growth temperature to two decimal places.
    '''
    #format data
    fasta_batch = SeqIO.parse(genome, "fasta")
    # loop through sequences counting IVYRWEL
    IVYRWEL_count = 0
    other_count = 0
    for record in fasta_batch:
        for item in record.seq:
            if item == "I":
                IVYRWEL_count = IVYRWEL_count +1
            elif item == "V":
                IVYRWEL_count = IVYRWEL_count +1
            elif item == "Y":
                IVYRWEL_count = IVYRWEL_count +1
            elif item == "W":
                IVYRWEL_count = IVYRWEL_count +1
            elif item == "R":
                IVYRWEL_count = IVYRWEL_count +1
            elif item == "E":
                IVYRWEL_count = IVYRWEL_count +1
            elif item == "L":
                IVYRWEL_count = IVYRWEL_count +1
            else:
                other_count = other_count +1
    #return values           
    length_count = IVYRWEL_count + other_count
    IVYWREL_fraction = (IVYRWEL_count/length_count)
    predicted_ogt = (937 * IVYWREL_fraction) - 335
    predicted_ogt_2dp = "{:.2f}".format(predicted_ogt)
    
    return predicted_ogt_2dp

def file_checker(genome_address, genome_file):
    '''
    Assesses the presence of the appropriate protein fasta file in the NCBI genome directory.
    
    Arguments:
        genome_address(string): the directory of the genome.
        genome_file(string): the file name you wish to analyse.
    Returns:
        file_present(boolean): If true, file is present. If false, file is not present.
    '''
    #set up boolean
    file_present = False
    #find ftp
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    # move to the right directory
    genome_directory = genome_address[28::] + "/"
    ftp.cwd(genome_directory) 
    #find the file
    if genome_file in ftp.nlst():
        file_present = True
    return file_present

    
def protein_file_generate(genome_address):
    '''
    Creates the protein.faa.gz NCBI ftp file path based upon the file path provided.
    if no file is detected by file_checker, it returns "No File"
    
    Arguments: 
        genome_address(string): The original ftp path to the genome directory
    Returns:
        actual_genome_address(string): The new ftp path to the appropriate protein file within the original directory
    '''
    reversed_genome_file = ""
    #find the file name
    for item in genome_address[::-1]:
        if item == "/":
            reversed_genome_file = reversed_genome_file + item
            break
        else:
            reversed_genome_file = reversed_genome_file + item
    #join the parts
    part_1_genome_file = ''.join(reversed(reversed_genome_file))
    part_2_genome_file = "_protein.faa.gz"
    actual_genome_address = genome_address + str(part_1_genome_file) + part_2_genome_file
    genome_file = str(part_1_genome_file) + part_2_genome_file
    short_genome_file = genome_file[1::]
    file_present = file_checker(genome_address, short_genome_file)
    if file_present == False:
        actual_genome_address = "No File"
    
    return actual_genome_address




def download_and_unzip_genome(actual_genome_address):
    '''
    Downloads the genome from the address provided and saves the file to the current directory. Unzips the file and saves to the current directory.
    Saves the files as "zipped_genome.txt" and "unzipped_genome.txt"
    
    Arguments:
        actual_genome_address(string): The file path to the file to be downloaded and unzipped.
    '''
    #dowload the protein genome
    zipped_genome = wget.download(actual_genome_address, out = "zipped_genome.txt")
    
    ## unzip

    with gzip.open(zipped_genome, 'rb') as f_in:
        with open('unzipped_genome.txt', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    

def create_predicted_vs_actual(data):
    '''
    Creates a dataset containing the predicted optimal growth temperature of microorganisms. Uses the dataset entered to find and download the relevent genome, and analyses the genome to predict OGT. Creates a new column in the dataset for the predicted OGT. Deletes the downloaded genome for each organism after calculating OGT.
    if no files are found, the predicted OGT is "N/A".
    
    Arguments:
        data(dataframe): The dataset containing the ftp path to the NCBI genome database, and the actual optimal growth temperature of each organism.
    Returns:
        output(dataframe): returns the original datafram with an additional column, containing the predicted OGT of each organism based upon their IVYWREL concentration.
    '''
    data["predicted_OGT"] = "empty"
    output = data
    ##errors
    if len(data["ftp_path"])<1:
        sys.exit(" Dataset has no ftp path column for genome download")

    
    for index, row in data.iterrows():
        #visualise what is happening
        print("--------------------------------")
        org_name = data.loc[index,"organism_name"]
        print(org_name)
        #get genome address
        genome_address = data.loc[index, "ftp_path"]
        actual_genome_address = protein_file_generate(genome_address)
        # is the file present in the directory?
        if actual_genome_address == "No File":
            predicted_ogt = "N/A"
            print("---could not be downloaded---")
        else:
            #download genome
            download_and_unzip_genome(actual_genome_address)
            #tell the computer its a fasta file
            unzipped_genome = 'unzipped_genome.txt'

            #predict ogt
            predicted_ogt = thermometer(unzipped_genome)
            
            ### remove the downloaded files
            #sometimes the files are still open, need time for files to close
            time.sleep(0.5)
            os.remove("zipped_genome.txt")
            os.remove("unzipped_genome.txt")
            print(f" predicted OGT is {predicted_ogt}")
        #put predicted ogt into dataset
        output.loc[index, "predicted_OGT"] = predicted_ogt
        
        
    return output
        





# get database from input as a csv file
# eg database.csv
input_name = sys.argv[1]
# get output name from input


output_name = 'metagenomic_thermometer_results.csv'
    
#import the genome list
data = pd.read_csv(f'{input_name}')

# do the program
new_data = create_predicted_vs_actual(data)

# save the results

new_data.to_csv(output_name)