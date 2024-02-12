
# Import necessary libraries
import pandas as pd
import re

# Function to clean the data
def clean_data(file_path):
    # Load the data
    data = pd.read_csv(file_path)

    # Remove brackets and words in brackets from the 'EG.PrecursorId' column
    data['EG.PrecursorId'] = data['EG.PrecursorId'].apply(lambda x: re.sub(r'\[.*?\]', '', x))

    # Remove '.#' from the 'EG.PrecursorId' column
    data['EG.PrecursorId'] = data['EG.PrecursorId'].apply(lambda x: re.sub(r'\.\d$', '', x))

    # Keep only the columns before 'EG.PrecursorId' and the 'EG.PrecursorId' column itself
    data = data.loc[:, :'EG.PrecursorId']

    # Save the cleaned data as a new CSV file
    data.to_csv('cleaned_data.csv', index=False)

# Call the function with the path to the file to be cleaned
clean_data('path_to_file.csv')
