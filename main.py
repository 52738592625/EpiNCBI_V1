#Install packages necessary to run the code
import subprocess
subprocess.run(
    f"python -m pip install pandas",
    shell=True,
    check=True
)
import os
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

####################################
#User selects a pathogen to download
print("Type the NCBI Taxonomy number of the pathogen you'd like to download, and press enter")
print("Escherichia coli = 562, Salmonella enterica = 590, etc.")
choice = input("")

####################################
#Download step. Filters out sequencing data. Cutoff date is Jan 1, 2024
print("Downloading. Please wait a while.")
subprocess.run(
    f"datasets download genome taxon {choice} --include none --released-before 2024-01-01",
    shell=True,
    check=True
)
print("Download successful")
####################################
#Data Formatting Step. Further filters out unused data fields


print("Formatting the data. Please wait a while.")
subprocess.run(
    f"dataformat tsv genome --package ncbi_dataset.zip --fields accession,assminfo-biosample-attribute-name,assminfo-biosample-attribute-value,assminfo-submitter > {choice}_data.tsv",
    shell=True,
    check=True
)
print("Data formatting successful")
if os.path.exists('ncbi_dataset.zip'):
    os.remove('ncbi_dataset.zip')
####################################
#Data Pivoting Step


print("Cleaning the data. Please wait a while.")
df = pd.read_csv(f'{choice}_data.tsv', sep='\t')
df.columns = ['genome', 'attribute', 'attribute_value', 'submitter']
df = df.pivot_table(index=['genome','submitter'], columns='attribute', values='attribute_value', aggfunc='first').reset_index()
df.columns.name = None
col_interest = ['genome', 'geo_loc_name','collection_date','strain','host', 'serotype','isolate','isolation_source','source_type', 'submitter']
df = df[col_interest]

####################################
#Data Cleaning Step

#Create a dictionary from the lookup table from the top US cities and their states
df_city_state_lookup = pd.read_csv('./lookup_tables/top_us_cities.csv', dtype=str)
df_city_state_lookup.head(5)
#create a dictionary with the states as keys and a list of the cities in that state as values
city_state_dict = {}
for index, row in df_city_state_lookup.iterrows():
    state = row['State']
    city = row['City']
    if state not in city_state_dict:
        city_state_dict[state] = []
    city_state_dict[state].append(city)
city_state_dict
#lowercase the keys and values in the dictionary
for state in city_state_dict:
    city_state_dict[state] = [city.lower() for city in city_state_dict[state]]
    state = state.lower()

#Create a column called is_usa_based full of NaN values initally
df['geo_loc_name'] = df['geo_loc_name'].fillna('')
df['geo_loc_name_is_missing'] = df['geo_loc_name'] == ''
df['is_usa_based'] = np.nan
df['usa_state'] = '-'

#read in the lookup table for US states
df_state_lookup = pd.read_csv('./lookup_tables/us_states.csv', dtype=str)
df_state_lookup.head(5)
state_abbreviation_dict = {}
for index, row in df_state_lookup.iterrows():
    state = row['State'].upper()
    abbreviation = row['Abbreviation'].upper()
    state_abbreviation_dict[abbreviation.upper()] = state

countries_df = pd.read_csv('./lookup_tables/countries_of_the_world.csv', dtype=str)
countries_list = countries_df['Country'].str.upper().tolist()

#set is_usa_based to 1 if geo_loc_name contains 'USA' or 'UNITED STATES'
df.loc[df['geo_loc_name'].str.upper().str.contains('USA|UNITED STATES'), 'is_usa_based'] = 1
for state in state_abbreviation_dict.values():
    df.loc[df['geo_loc_name'].str.upper().str.contains(state), 'is_usa_based'] = 1
#if geo_loc_name contains a string part that is in countries_list set is_usa_based to 0
for country in countries_list:
    df.loc[df['geo_loc_name'].str.upper().str.contains(country), 'is_usa_based'] = 0
#If the geo_loc_name is in the list of countries, set is_usa_based to 0
df.loc[df['geo_loc_name'].str.upper().isin(countries_list), 'is_usa_based'] = 0
#If the geo_loc_name.upper() == 'MISSING' or 'UNKNOWN', set usa_based to NaN
df.loc[df['geo_loc_name'].str.upper().isin(['MISSING', 'UNKNOWN', 'NOT DETERMINED', 'NOT APPLICABLE', 'NOT COLLECTED', 'NA', 'NOT PROVIDED', '', 'NOT AVAILABLE', 'RESTRICTED ACCESS']), 'is_usa_based'] = 99
#all others set is_usa_based to 98
df.loc[df['is_usa_based'].isna(), 'is_usa_based'] = 98

#if geo_loc_name == 'USA' or 'UNITED STATES', set usa_state to 'Unclassified'
df.loc[df['geo_loc_name'] == 'USA', 'usa_state'] = 'Unknown'
#where is_usa_based == 1, if geo_loc_name contains a string part that is in the city_state_dict, set usa_state to the state
for state, cities in city_state_dict.items():
    for city in cities:
        df.loc[df['geo_loc_name'].str.upper().str.contains(" " + city.upper()), 'usa_state'] = state
        df.loc[df['geo_loc_name'].str.upper().str.contains(":" + city.upper()), 'usa_state'] = state
#where is_usa_based == 1, if geo_loc_name contains a string part that is in the state_abbreviations_dict, set usa_state to the state
for abbreviation, state in state_abbreviation_dict.items():
    df.loc[df['geo_loc_name'].str.upper().str.contains(" " + abbreviation), 'usa_state'] = state
    df.loc[df['geo_loc_name'].str.upper().str.contains(" " + state), 'usa_state'] = state
    df.loc[df['geo_loc_name'].str.upper().str.contains(":" + abbreviation), 'usa_state'] = state
    df.loc[df['geo_loc_name'].str.upper().str.contains(":" + state), 'usa_state'] = state
df.loc[(df['is_usa_based'] == 1) & (df['usa_state'] == '-'), 'usa_state']  = 'No state'
df.loc[df['is_usa_based'].isin({0, 98, 99}), 'usa_state'] = 'N/A'
#remove all rows where is_usa_based is 0, 98, or 99, leaving only USA based rows
df = df[df['is_usa_based'] == 1]

df.set_index('genome', inplace=True)
df['collection_date'] = df['collection_date'].fillna('missing')
df['collection_year'] = df['collection_date'].replace(
    ['Missing', 'Unknown', 'unknown', 'none', 'not applicable', 'not collected', 'not provided', 'Not collected', 'MISSING', 'not available', 'Not Collected', 'unspecified', 'Not available'],
    'missing'
)
years = [str(x) for x in range(1900, 2026)]

def extract_year(date_str):
    if date_str == 'missing':
        return 'missing'
    r = date_str[-4:]
    l = date_str[:4]
    if r in years:
        return r
    elif l in years:
        return l
    else:
        return 'missing'

df['collection_year'] = df['collection_date'].apply(extract_year)
df['month'] = df['collection_date'].apply(lambda x: x[5:7] if len(x) > 5 and x[4] == '-' else None)
df['month'] = pd.to_numeric(df['month'], errors='coerce')


##############################################
#Saving the data output

df.to_csv(f'./ncbi_{choice}_data.csv', index=True)

if os.path.exists(f'{choice}_data.tsv'):
    os.remove(f'{choice}_data.tsv')
if os.path.exists(f'{choice}_unstacked.csv'):
    os.remove(f'{choice}_unstacked.csv')
print("Done")

