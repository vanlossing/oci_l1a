#!/usr/bin/env python

from pprint import pprint

import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

file = "/Users/cfield/Downloads/intervals.csv"

df = pd.read_csv(file)
print(df.columns)
df.t_id_start = pd.to_datetime(df.t_id_start, unit='s')
df.t_start = pd.to_datetime(df.t_start, unit='us')
df.t_end = pd.to_datetime(df.t_end, unit='us')
df["meta"]  = [eval(x) for x in df.meta]

pprint(df, compact=False)

print(type(df.loc[0, "meta"]))
print(type(df.loc[24, "meta"]))


np.ma.compress_cols
