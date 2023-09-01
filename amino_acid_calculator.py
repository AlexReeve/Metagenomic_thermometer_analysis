import pandas as pd
from Bio import SeqIO
import wget
import gzip
import shutil
import os
from ftplib import FTP
import time
import sys

def amino_acid_count(genome):
    '''
    Counts all amino acids in genome and returns their fraction of the total genome.
    
    
    Arguments:
        genome(protein fasta file): The genome. Must be a protein fasta file. 
    
    Returns:
        I_fraction(int):The fraction of amino acid I present in the genome, up to five decimal places.
        V_fraction(int):The fraction of amino acid V present in the genome, up to five decimal places.
        Y_fraction(int):The fraction of amino acid Y present in the genome, up to five decimal places.
        W_fraction(int):The fraction of amino acid W present in the genome, up to five decimal places.
        R_fraction(int):The fraction of amino acid R present in the genome, up to five decimal places.
        E_fraction(int):The fraction of amino acid E present in the genome, up to five decimal places.
        L_fraction(int):The fraction of amino acid L present in the genome, up to five decimal places.
        A_fraction(int):The fraction of amino acid A present in the genome, up to five decimal places.
        N_fraction(int):The fraction of amino acid N present in the genome, up to five decimal places.
        D_fraction(int):The fraction of amino acid D present in the genome, up to five decimal places.
        C_fraction(int):The fraction of amino acid C present in the genome, up to five decimal places.
        Q_fraction(int):The fraction of amino acid Q present in the genome, up to five decimal places.
        G_fraction(int):The fraction of amino acid G present in the genome, up to five decimal places.
        H_fraction(int):The fraction of amino acid H present in the genome, up to five decimal places.
        K_fraction(int):The fraction of amino acid K present in the genome, up to five decimal places.
        M_fraction(int):The fraction of amino acid M present in the genome, up to five decimal places.
        F_fraction(int):The fraction of amino acid F present in the genome, up to five decimal places.
        P_fraction(int):The fraction of amino acid P present in the genome, up to five decimal places.
        S_fraction(int):The fraction of amino acid S present in the genome, up to five decimal places.
        T_fraction(int):The fraction of amino acid T present in the genome, up to five decimal places.
        B_fraction(int):The fraction of amino acid B present in the genome, up to five decimal places.
        Z_fraction(int):The fraction of amino acid Z present in the genome, up to five decimal places.
        O_fraction(int):The fraction of amino acid O present in the genome, up to five decimal places.
        U_fraction(int):The fraction of amino acid U present in the genome, up to five decimal places.
        X_fraction(int):The fraction of X present in the genome, up to five decimal places.
        J_fraction(int):The fraction of amino acid J present in the genome, up to five decimal places.
        other_fraction(int):The fraction of other amino acids present in the genome, up to five decimal places.
    '''
    #format data
    fasta_batch = SeqIO.parse(genome, "fasta")
    # loop through sequences 
    I_count = 0
    V_count = 0
    Y_count = 0
    W_count = 0
    R_count = 0
    E_count = 0
    L_count = 0
    A_count = 0
    N_count = 0
    D_count = 0
    C_count = 0
    Q_count = 0
    G_count = 0
    H_count = 0
    K_count = 0
    M_count = 0
    F_count = 0
    P_count = 0
    S_count = 0
    T_count = 0
    B_count = 0
    Z_count = 0
    O_count = 0
    U_count = 0
    X_count = 0
    J_count = 0
    other_count = 0
    sequence_count = 0
    for record in fasta_batch:
        for item in record.seq:
            sequence_count = sequence_count+1
            if item == "I":
                I_count = I_count +1
            elif item == "V":
                V_count = V_count +1
            elif item == "Y":
                Y_count = Y_count +1
            elif item == "W":
                W_count = W_count +1
            elif item == "R":
                R_count = R_count +1
            elif item == "E":
                E_count = E_count +1
            elif item == "L":
                L_count = L_count +1
            elif item == "A":
                A_count = A_count +1
            elif item == "N":
                N_count = N_count +1
            elif item == "D":
                D_count = D_count +1
            elif item == "C":
                C_count = C_count +1
            elif item == "Q":
                Q_count = Q_count +1
            elif item == "G":
                G_count = G_count +1
            elif item == "H":
                H_count = H_count +1
            elif item == "K":
                K_count = K_count +1
            elif item == "M":
                M_count = M_count +1
            elif item == "F":
                F_count = F_count +1
            elif item == "P":
                P_count = P_count +1
            elif item == "S":
                S_count = S_count +1
            elif item == "T":
                T_count = T_count +1
            elif item == "B":
                B_count = B_count +1
            elif item == "Z":
                Z_count = Z_count +1
            elif item == "O":
                O_count = O_count +1
            elif item == "U":
                U_count = U_count +1
            elif item == "X":
                X_count = X_count +1
            elif item == "J":
                J_count = J_count +1
            else:
                other_count = other_count +1
            
    #return values           
    
    
    I_fraction = "{:.5f}".format(I_count/sequence_count)
    V_fraction = "{:.5f}".format(V_count/sequence_count)
    Y_fraction = "{:.5f}".format(Y_count/sequence_count)
    W_fraction = "{:.5f}".format(W_count/sequence_count)
    R_fraction = "{:.5f}".format(R_count/sequence_count)
    E_fraction = "{:.5f}".format(E_count/sequence_count)
    L_fraction = "{:.5f}".format(L_count/sequence_count)
    A_fraction = "{:.5f}".format(A_count/sequence_count)
    N_fraction = "{:.5f}".format(N_count/sequence_count)
    D_fraction = "{:.5f}".format(D_count/sequence_count)
    C_fraction = "{:.5f}".format(C_count/sequence_count)
    Q_fraction = "{:.5f}".format(Q_count/sequence_count)
    G_fraction = "{:.5f}".format(G_count/sequence_count)
    H_fraction = "{:.5f}".format(H_count/sequence_count)
    K_fraction = "{:.5f}".format(K_count/sequence_count)
    M_fraction = "{:.5f}".format(M_count/sequence_count)
    F_fraction = "{:.5f}".format(F_count/sequence_count)
    P_fraction = "{:.5f}".format(P_count/sequence_count)
    S_fraction = "{:.5f}".format(S_count/sequence_count)
    T_fraction = "{:.5f}".format(T_count/sequence_count)
    B_fraction = "{:.5f}".format(B_count/sequence_count)
    Z_fraction = "{:.5f}".format(Z_count/sequence_count)
    O_fraction = "{:.5f}".format(O_count/sequence_count)
    U_fraction = "{:.5f}".format(U_count/sequence_count)
    X_fraction = "{:.5f}".format(X_count/sequence_count)
    J_fraction = "{:.5f}".format(J_count/sequence_count)
    other_fraction = "{:.5f}".format(other_count/sequence_count)
    
    
    return I_fraction, V_fraction, Y_fraction, W_fraction, R_fraction, E_fraction, L_fraction, A_fraction, N_fraction, D_fraction, C_fraction, Q_fraction, G_fraction, H_fraction, K_fraction, M_fraction, F_fraction, P_fraction, S_fraction, T_fraction, B_fraction, Z_fraction, O_fraction, U_fraction, X_fraction, J_fraction, other_fraction

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
    

def create_dataset(data):
    '''
    Creates a dataset containing the fractions of amino acids present in the genome of microorganisms. Uses the dataset entered to find and download the relevent genome. Creates a new column in the dataset for each amino acid. Deletes the downloaded genome for each organism after calculating OGT.
    if no files are found, the amino acid fraction is "N/A".
    
    Arguments:
        data(dataframe): The dataset containing the ftp path to the NCBI genome database, and the actual optimal growth temperature of each organism.
    Returns:
        output(dataframe): returns the original datafram with an additional 28 columns, containing the amino acids (and other potential letters present) of each organismn.
    '''
    data["I_fraction"]= "empty"
    data["V_fraction"]= "empty"
    data["Y_fraction"]= "empty"
    data["W_fraction"]= "empty"
    data["R_fraction"]= "empty"
    data["E_fraction"]= "empty"
    data["L_fraction"]= "empty"
    data["A_fraction"]= "empty"
    data["N_fraction"]= "empty"
    data["D_fraction"]= "empty"
    data["C_fraction"]= "empty"
    data["Q_fraction"]= "empty"
    data["G_fraction"]= "empty"
    data["H_fraction"]= "empty"
    data["L_fraction"]= "empty"
    data["K_fraction"]= "empty"
    data["M_fraction"]= "empty"
    data["F_fraction"]= "empty"
    data["P_fraction"]= "empty"
    data["S_fraction"]= "empty"
    data["T_fraction"]= "empty"
    data["B_fraction"]= "empty"
    data["Z_fraction"]= "empty"
    data["O_fraction"]= "empty"
    data["U_fraction"]= "empty"
    data["X_fraction"]= "empty"
    data["J_fraction"]= "empty"
    data["other_fraction"]= "empty"
    output = data
    count = 0
    ## error message for no ftp path column
    if len(data["ftp_path"])<1:
        sys.exit(" Dataset has no ftp path column for genome download")
        
        
    ######
    for index, row in data.iterrows():
        ## visualise for me
        print("--------------------------------")
        org_name = data.loc[index,"organism_name"]

        #get genome address
        genome_address = data.loc[index, "ftp_path"]
        actual_genome_address = protein_file_generate(genome_address)
        count = count +1
        # is the file present in the directory?
        if actual_genome_address == "No File":
            I_fraction= "N/A"
            V_fraction= "N/A"
            Y_fraction= "N/A"
            W_fraction= "N/A"
            R_fraction= "N/A"
            E_fraction= "N/A"
            L_fraction= "N/A"
            A_fraction= "N/A"
            N_fraction= "N/A"
            D_fraction= "N/A"
            C_fraction= "N/A"
            Q_fraction= "N/A"
            G_fraction= "N/A"
            H_fraction= "N/A"
            L_fraction= "N/A"
            K_fraction= "N/A"
            M_fraction= "N/A"
            F_fraction= "N/A"
            P_fraction= "N/A"
            S_fraction= "N/A"
            T_fraction= "N/A"
            B_fraction= "N/A"
            Z_fraction= "N/A"
            O_fraction= "N/A"
            U_fraction= "N/A"
            X_fraction= "N/A"
            J_fraction= "N/A"
            other_fraction= "N/A"
            print("---could not be downloaded---")
            
        else:
            #download genome
            download_and_unzip_genome(actual_genome_address)
            #tell the computer its a fasta file
            unzipped_genome = 'unzipped_genome.txt'

            #get amino acid count
            I_fraction, V_fraction, Y_fraction, W_fraction, R_fraction, E_fraction, L_fraction, A_fraction, N_fraction, D_fraction, C_fraction, Q_fraction, G_fraction, H_fraction, K_fraction, M_fraction, F_fraction, P_fraction, S_fraction, T_fraction, B_fraction, Z_fraction, O_fraction, U_fraction, X_fraction, J_fraction, other_fraction = amino_acid_count(unzipped_genome)
            
            ### remove the downloaded files
            #sometimes the files are still open, need time for files to close
            time.sleep(0.5)
            os.remove("zipped_genome.txt")
            os.remove("unzipped_genome.txt")
            
        #put amino acids into dataset
        output.loc[index, "I_fraction"] = I_fraction
        output.loc[index, "V_fraction"] = V_fraction
        output.loc[index, "Y_fraction"] = Y_fraction
        output.loc[index, "W_fraction"] = W_fraction
        output.loc[index, "R_fraction"] = R_fraction
        output.loc[index, "E_fraction"] = E_fraction
        output.loc[index, "L_fraction"] = L_fraction
        output.loc[index, "A_fraction"] = A_fraction
        output.loc[index, "N_fraction"] = N_fraction
        output.loc[index, "D_fraction"] = D_fraction
        output.loc[index, "C_fraction"] = C_fraction
        output.loc[index, "Q_fraction"] = Q_fraction
        output.loc[index, "G_fraction"] = G_fraction
        output.loc[index, "H_fraction"] = H_fraction
        output.loc[index, "K_fraction"] = K_fraction
        output.loc[index, "M_fraction"] = M_fraction
        output.loc[index, "F_fraction"] = F_fraction
        output.loc[index, "P_fraction"] = P_fraction
        output.loc[index, "S_fraction"] = S_fraction
        output.loc[index, "T_fraction"] = T_fraction
        output.loc[index, "B_fraction"] = B_fraction
        output.loc[index, "Z_fraction"] = Z_fraction
        output.loc[index, "O_fraction"] = O_fraction
        output.loc[index, "U_fraction"] = U_fraction
        output.loc[index, "X_fraction"] = X_fraction
        output.loc[index, "J_fraction"] = J_fraction
        output.loc[index, "other_fraction"] = other_fraction
        print(f" xxxxxxxx Succesfully Downloaded and Calculated {org_name} xxxxxxxx")
  
       
    return output
        


# get database from input as a csv file
# eg database.csv
input_name = sys.argv[1]

output_name = 'amino_acid_count.csv'
    
#import the genome list
data = pd.read_csv(f'{input_name}')

# do the program
new_data = create_dataset(data)

# save the results

new_data.to_csv(output_name)