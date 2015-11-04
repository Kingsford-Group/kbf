k-mer Bloom filters
------------------------------
# (REAME under construction)

C++ implementation of the k-mer Bloom filters described in: Pellow, Filippova, and Kingsford, "Improving Bloom filter performance on sequence data using k-mer Bloom filters"

The C++ source files are in the directory cpp-src

Implementations of one-sided, two-sided, sparse (strict contains function) and sparse (relaxed contains function) kbfs are in the files KBF1.hpp, KBF2.hpp, KBFSparse.hpp and KBFSparseRelaxed.hpp

#### Dependencies:

Boost: http://www.boost.org/doc/libs/1_59_0/more/getting_started/unix-variants.html

libbf: https://github.com/mavam/libbf/ 

###### Compilation:
```
g++ -std=c++11 -O3 -o kbf main.cpp -I [libbf include path] -I [Boost include path] -L [libbf lib path] -L [Boost lib path] -lbf
```

