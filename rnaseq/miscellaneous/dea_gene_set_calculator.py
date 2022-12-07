#!/usr/bin/env python3
'''
input: 
    a_vs_b.xls
    c_vs_b.xls
output:
    co-up.xls
    co-down.xls

当前只支持读取gfold比较结果。
'''

import sys
import pandas as pd

def read_gfold(path):
    return pd.read_csv(path, sep='\t', index_col=0)

def get_gfold_up(data, gfold=2, rpkm=10):
    c1 = data["GFOLD(0.01)"] > gfold
    c2 = (data['1stRPKM'] + data['2ndRPKM'])/2 > rpkm
    r = data[c1 & c2]
    return set(r.index)

def get_gfold_down(data, gfold=-2, rpkm=10):
    c1 = data["GFOLD(0.01)"] < gfold
    c2 = (data['1stRPKM'] + data['2ndRPKM'])/2 > rpkm
    r = data[c1 & c2]
    return set(r.index)

# 1. read
# 2. get up or down regulated gene sets. 
# 3. set operation

# table1 = "gfold/AC_last_vs_MS2_last.anno.xls"
# table2 = "gfold/AC_vs_DYR.anno.xls"
table1 = sys.argv[1]
table2 = sys.argv[2]

table1 = read_gfold(table1)
table2 = read_gfold(table2)

table1_up = get_gfold_up(table1)
table2_up = get_gfold_up(table2)
table1_down = get_gfold_down(table1)
table2_down = get_gfold_down(table2)

co_up = table1_up & table2_up
co_down = table1_down & table2_down
print("Num of up regulated genes in table1:", len(table1_up))
print("Num of up regulated genes in table2:", len(table2_up))
print("Num of up regulated genes in both tables:", len(co_up))
print()
print("Num of down regulated genes in table1:", len(table1_down))
print("Num of down regulated genes in table2:", len(table2_down))
print("Num of down regulated genes in both tables:", len(co_down))

print("co-up list: ")
for i in co_up:
    print(i)

print("co-down list: ")
for i in co_down:
    print(i)
