# npcd
Neighbourhood Proximity based Community Detection (NPCD)
To to able to use this code, first compile it using the C++ compiler on a Unix machine as follows:
g++ -std=c++11 npcd.cpp -o npcd
Then an executable file named npcd will be created. This can be executed, for weighted graphs, as follows:
./npcd network_file -rh [option] -ov [option] -w
For unweighted graphs, you need to simply run it as follows:
./npcd network_file -rh [option] -ov [option]

Here, -rh, and -ov are the algorithm parameters whose values are to be specified at the place of option. 
If you use this code for research purpose, please cite the following article:

Kumar, P., & Dohare, R. (2019). A neighborhood proximity based algorithm for overlapping community structure detection in weighted networks. Frontiers of computer science, 13(6), 1353-1355.
