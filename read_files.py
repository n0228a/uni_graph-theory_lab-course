import glob
import pandas as pd
import random

for file in glob.glob("data/*.xlsx"):
    data = pd.read_excel(file)
    print(f"File: {file}, Rows: {len(data)}")
    # Print the first 20 rows
    print(data.head(20))