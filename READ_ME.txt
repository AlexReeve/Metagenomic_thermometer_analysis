ReadMe File

I created this repositoty as part of my Bioinformatics masters at the University of Bristol.
The following is an analysis of the metagenomic thermometer of Kurakawa et al (2023). I also analyse Zeldovich et al's(2007) IVYWREL predictor of OGT.
In this file are all datasets used in the study, and all datasets needed for the metagenomic_thermometer_analysis_workbook.ipynb.

Before using any of the python files below, open the metagenomic_thermometer_analysis_workbook.ipynb. It will take you throughout all stages of analysis, explain each dataset, and tell you how to use each python program.

-------------------------------------------------------------------------------------------------------------------------------

This file conatains the following items:

amino_acid_and_met_thermometer.py - python program.From a file of genome summaries, uses the ftp path to download the genomes. Counts the IVYWREL fraction and predicts the OGT. Also counts every amino acid fraction. Program is a combination of the 'original_metagenomic_thermometer.py' program and 'amino_acid_calculator.py' in this file. Combined for conveniance.
                                    Use: python amino_acid_and_met_thermometer.py input_example_dataset.csv

amino_acid_calculator.py - python program.From a file of genome summaries, uses the ftp path to download the genomes. Calculates every amino acid fraction in each genome.
                                    Use: python amino_acid_calculator.py input_example_dataset.csv

archaea_complete_genomes.csv - csv dataset. The list of all complete archaean genomes available on the NCBI genome database.

archaea_dataset_with_ogt_and_genomes.csv - csv dataset. Contains all complete archaean genomes available on the NCBI genome database with a known optimal growth temperature.

bact_dataset_with_ogt_and_genomes.csv - csv dataset. Contains all complete bacterial genomes available on the NCBI genome database with a known optimal growth temperature.

bacteria_complete_genomes.csv - csv dataset. The list of all complete bacterial genomes available on the NCBI genome database.

find_best_correlation.py - python program. From a file of genome summaries with known optimal growth temperatures and amino acid genome fractions, finds the best amino acid set for predicting OGT. See workbook.
			Use: python find_best_correlation.py input_example_dataset.csv

full_archaea_amino_acids.csv - csv dataset. Contains all available complete archaean genomes from NCBI with known OGTS. Contains all amino acid fractions for each genome.

full_archaea_predicted_ogt.csv - csv dataset. Contains all available archaean genomes from NCBI with known ogts. Contains predicted OGTS, calculated using 'original_metagenomic_thermometer', and the Zeldovich et al (2007) equation.

full_bacteria_amino_acids.csv - csv dataset. Contains all available complete bacterial genomes from NCBI with known OGTS. Contains all amino acid fractions for each genome.

full_bacteria_predicted_ogt.csv - csv dataset. Contains all available bacterial genomes from NCBI with known ogts. Contains predicted OGTS, calculated using 'original_metagenomic_thermometer', and the Zeldovich et al (2007) equation.

full_prokaryote_amino_acids.csv - csv dataset. Contains all available complete prokaryote genomes from NCBI with known OGTS. Contains all amino acid fractions for each genome.

full_prokaryote_predicted_ogt.csv - csv dataset. Contains all available prokaryote genomes from NCBI with known ogts. Contains predicted OGTS, calculated using 'original_metagenomic_thermometer', and the Zeldovich et al (2007) equation.

marry_find.py - python program. Using 'temperature_data.csv'(Engqvist, 2018), takes a dataset of genomes and produces a dataset of genomes with known OGTS.
		Use: python marry_find.py input_database_of_genomes.csv

metagenomic_thermometer_analysis_workbook.ipynb - Jupyter workbook. Read this workbook before attempting to use other datasets and python files. Contains all steps in the metagenomic thermometor analysis.

original_metagenomic_thermometer.py - python file - From an input dataset of genomes, downloads each genome using the FTP path. Calulates the predicted OGT for each genome using the IVYWREL fraction and the equation from Zeldovich et al(2007).

representative_archaea_amino_acids.csv - csv file. Contains the 'representative' archaean dataset used in the report. Contains all amino acid fractions for each genome.

representative_archaea_predicted_ogt.csv- csv file. Contains the 'representative' archaean dataset used in the report. Contains predicted OGTS, calculated using 'original_metagenomic_thermometer', and the Zeldovich et al (2007) equation.

representative_bacteria_amino_acids.csv - csv file. Contains the 'representative' bacterial dataset used in the report. Contains all amino acid fractions for each genome.

representative_bacteria_predicted_ogt.csv -csv file. Contains the 'representative' bacterial dataset used in the report. Contains predicted OGTS, calculated using 'original_metagenomic_thermometer', and the Zeldovich et al (2007) equation.

representative_prokaryote_amino_acids.csv - csv file. Contains the 'representative' prokaryotic dataset used in the report. Contains all amino acid fractions for each genome.

representative_prokaryote_predicted_ogt.csv - csv file. Contains the 'representative' prokaryotic dataset used in the report. Contains predicted OGTS, calculated using 'original_metagenomic_thermometer', and the Zeldovich et al (2007) equation.

temperature_data.csv - csv file. Growth temperatures for 21,498 microorganisms. Engqvist, 2018.
