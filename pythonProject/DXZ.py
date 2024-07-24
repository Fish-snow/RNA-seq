import pandas as pd
df=pd.read_excel("/Users/yuxueqian/Desktop/test.xlsx")
df=df.drop(df.index[[6,7,8]])
