#k-mer Bloom filters
------------------------------

Sequence k-mer content is often used to compare sequences, enabling significant 
performance improvements in metagenomic species identification, estimation of transcript abundances, 
and alignment-free comparison of sequencing data. k-mer sets often reach hundreds 
of millions of elements making traditional data structures impractical for k-mer set storage.
Probabilistic Bloom filters and their variants are often used instead. 
Bloom filters reduce the memory footprint and allow for fast set containment queries.
Since k-mers are derived from sequencing reads, the information about k-mer overlap 
can be used to reduce the false positive rate up to two orders of magnitude 
with little or no additional memory and with set containment queries that are 
1.3 - 1.6 times slower. 
Alternatively, we can leverage k-mer overlap information to store k-mer sets in about half 
the space while maintaining the original false positive rate. 

More details are available at:

``` 
Pellow, Filippova, and Kingsford. "Improving Bloom filter performance on sequence data using k-mer Bloom filters" To appear in RECOMb 2016.
```

--------

#### Dependencies

Boost: http://www.boost.org/doc/libs/1_59_0/more/getting_started/unix-variants.html

libbf: https://github.com/mavam/libbf/ 

#### Compilation

```
g++ -std=c++11 -O3 -o kbf main.cpp -I [libbf include path] -I [Boost include path] -L [libbf lib path] -L [Boost lib path] -lbf
```

#### *k*BF Variants

The C++ source files are in the directory `cpp-src`

* KBF1.hpp -- one-sided Bloom filter that improves false postive rate three fold without using any additional storage

* KBF2.hpp -- two-sided Bloom filter that improve FPR by an order of magnitude while using very little additional memory

* KBFSparse.hpp -- sparse Bloom filter with a strict `contains` function that uses 1/2 of space to store the same set of kmers and guarantees the same FPR as a classic Bloom filter

* KBFSparseRelaxed.hpp -- same as above, but `contains` is relaxed

#### Example usage

```
string kmer = "ACGTACGTACGTACGTACGT";
kmer_t binary_kmer = mer_string_to_binary(kmer, 20);
KBF1 kbf(20 /* kmer length */, read_kmers /* collection of kmers to store */);
bool in_set = kbf.contains(binary_kmer);
cerr << "Kmer " << kmer << ": " << in_set << endl;
```
