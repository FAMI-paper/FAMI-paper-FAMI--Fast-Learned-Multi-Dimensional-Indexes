# FAMI
## About this repo
- `include`: codes of FAMI
- `libmorton`: open source codes of Z-order  https://github.com/Forceflow/libmorton
- `libspatialindex`: open source codes of R*-tree  https://github.com/libspatialindex/libspatialindex
- `stx-btree-0.9`: open source codes of B+ tree  https://github.com/learnedsystems/SOSD/tree/master/competitors/stx-btree-0.9

## Requirements
- g++

## BUILD
g++ main.cpp -o main  -I. -l spatialindex -lpthread -O2
## RUN
./main <data_file> <data_num> <query_num> <dimension> <exp> <index_type> <num_radix_bits> <cell_size> <rs_group_size> <query_type> <range_size> <update_ratio>

example: ./main testdata.txt 1000000 1000 2 7 grid+r-tree 18 100 200 1 400 1
+ <data_file>: Name of data file
+ <data_num>: Number of data used, be careful not to exceed the total number in data file
+ <query_num>: Number of querie
+ <dimension>: Dimension of Data
+ <exp>：Used to transfer float data to integer data, integer data = round(float data * 10^exp)
+ <index_type>: Index types, including: grid+r-tree, avggrid+r-tree, r*-tree, avggrid+b-tree
+ <num_radix_bits>: The size of all radix tables, including one on the x-axis and n on the y-axis (n is the number of segments divided on the x-axis), each radix_table is equal in size 
+ <cell_size>: The number of data points in each grid 
+ <rs_group_size>: The number of data points in each group on the y-axis 
+ <query_type>: Query type: 0 for point query, 1 for range query, 2 for update query 
+ <range_size>: Range query size, the size in each dimension is 1/<range_size> of the entire range 
+ <update_ratio>: The percentage of updated data points

+ ## DATASET
+ TPCH: https://drive.google.com/file/d/1Qdq5acHdwKmOSptPl0JqB-G_F44Eynpq/view?usp=sharing
