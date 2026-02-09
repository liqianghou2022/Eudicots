import pandas as pd

# Read CSV file
df = pd.read_csv('37-ACEK-202.id.csv')

# Transpose and reset index
df_transposed = df.T.reset_index()

# Rename columns
df_transposed.columns = ['new_column'] + list(range(1, len(df_transposed.columns)))

# Save transposed data
df_transposed.to_csv('37-ACEK-202_transposed.csv', index=False)