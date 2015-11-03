KBF (README under construction)
------------------------------
C++ implementation of the k-mer Bloom filters described in: Pellow, Filippova, and Kingsford, "Improving Bloom filter performance on sequence data using k-mer Bloom filters"

#### Dependencies:

Boost: 

libbf: 

###### Compilation:

g++ -std=c++11 -O3 -o kbf main.cpp -I <libbf include path> -I <Boost include path> -L <libbf lib path> -L <Boost lib path> -lbf
