import pandas as pd
import pyarrow

HOMEPATH="/home/fridald4/projects/def-gsarah/fridald4/renal_genetics_project/MR_analysis/"

df = pd.read_parquet(HOMEPATH + 'Nerve_Tibial.v8.EUR.allpairs.chr3.parquet')
df.to_csv(HOMEPATH + 'Nerve_Tibial.v8.EUR.allpairs.chr3.csv')
