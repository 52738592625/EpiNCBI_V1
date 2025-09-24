This readme documents the procedure for running the code to generate the results shown in "Leveraging NCBI Genomic Metadata for Epidemiological Insights:
Example of Enterobacterales".

The code that carries out the data pipeline in Figure 1 is detailed in main.py. To run the code, make sure you have Python downloaded. You can download Python at https://www.python.org/downloads/
Follow the default installation instructions provided for your machine.

Downloading the GitHub repo can be done via cloning the repo through Git, or by direct download of a zip file of the repo, and extracting the contents somewhere on the machine.

The code that carries out the data pipeline in Figure 1 is detailed in main.py. There are many ways
to run this code. For example, after installing Python on a Windows PC, and extracting a zip download of the repo,
navigate to the newly extracted folder
containing the files. click the search bar at the top of file explorer, type "cmd", and press
enter to open the command prompt. From the command prompt, type: "python main.py" and press enter.
The code will automatically install any required Python packages to carry out the functions in the code.

The program will then prompt the user to enter the NCBI Taxonomy ID of the pathogen they'd like to download, and to
press enter. IDs can be looked up for individual pathogens on https://www.ncbi.nlm.nih.gov/datasets/

E. Coli's ID is 562, and Salmonella's is 590.

The code will then begin the download, and conclude with cleaning and saving the output data to the same directory
main.py resides. Due to the limited bandwidth of NCBI's databases, it may take several minutes for the data to
download.

Pre-downloaded copies of E. Coli and Salmonella NCBI data, ncbi_562_data.csv and ncbi_590_data.csv respectively, as well as a copy of the NORS data used for analysis, NORS_20250702.csv, are provided in the root folder as well.
