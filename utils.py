
import os
import pandas as pd

from itertools import product

'''
MAKE SURE THIS FILE IS NEXT TO YOUR NOTEBOOK, AND YOUR DATASET FOLDER IS THERE TOO <3

Usage:

from utils import get_genome_examples

# this will grab the first 10_000 from each class
train_test_df = get_genome_examples(directory_name = 'trainingdata', num_samples = 10_000, delay_start = 0)

# this will grab 200_000 examples and IGNORE the 10_000 examples grabbed above 
validation_df = get_genome_examples(directory_name = 'trainingdata', num_samples = 200_000, delay_start = 10_000)


**************************************************
Args:
    - directory_name - the name of the directory of downloaded genome sequences; the default is 'trainingdata'
    - num_samples    - the number of samples to take from each class. If you use 10_000, it will take 10_000 from each class,
                       so 60_000 total for all 6 classes
    - delay_start    - the number of samples to ignore. If you create a train/test split, you can use delay_start to ignore the
                       train/test split and start grabbing samples afterwards 

'''

def get_genome_examples(directory_name: str = 'trainingdata',
                        num_samples: int = 60_000, delay_start: int = 0) -> pd.DataFrame:
    
    # ensures the data path is in ./trainingdata
    DATA_PATH = os.path.join(os.getcwd(), directory_name)

    X = pd.DataFrame([], dtype=str)

    # loop through virus names
    for virus_name in os.listdir(DATA_PATH):

        virus_series = pd.DataFrame([], dtype=str)
        virus_delay = delay_start
        
        # loop through hasta files
        for fasta in os.listdir(os.path.join(DATA_PATH, virus_name)): 
            file_name = os.path.join(DATA_PATH, virus_name, fasta)  

            # determine how many samples are needed
            needed_samples = num_samples - len(virus_series) + virus_delay
    
            # read hasta files and fetch the correct number of samples of genome
            df = pd.read_csv(file_name, header=None, names=['genome'], 
                             nrows=needed_samples, skiprows=lambda i: not i % 2).squeeze("columns").to_frame()
            df['label'] = virus_name  # assign the labels
    
            # throw away the samples that are used for training/testing purposes
            if virus_delay > 0:
                if len(df) > virus_delay:
                    df = df.iloc[virus_delay:]  # get rid of the examples
                    virus_delay = 0        
                # of not enough samples in hasta, just go to next file
                else:
                    virus_delay -= len(df)
                    continue
            
            virus_series = pd.concat([virus_series, df])

            # check if we have enough samples from this virus class
            if len(virus_series) >= num_samples:
                break  

        X = pd.concat([virus_series, X])        
    
    print(f"Created set with {len(X)} samples")
    
    return X


def get_vocab(n_gram: int = 5):
    vocab = dict()
    for p in product('ACTG', repeat=n_gram):
        vocab[' '.join(p)] = len(vocab)
    return vocab