/*************************************************************
* Copyright (c) David Pellow, Darya Filippova, Carl Kingsford
*************************************************************/

#ifndef LIB_BASE_BLOOM_FILTER
#define LIB_BASE_BLOOM_FILTER

#include <string>
#include <bf.h>

using namespace std;

// kmer type
typedef uint64_t kmer_t;

class BaseBloomFilter {
protected:
	bf::basic_bloom_filter bf_;	// generic bloom filter for holding the entries
	const int k; // kmer length
	size_t num_inserted = 0; // number of kmers inserted
	const size_t bf_size; // size of the bf
	const int h; // number of hash functions

public:

	BaseBloomFilter(const int k, const size_t num_elems=1024*1024*32, const size_t size_factor=10,const int num_hashes = 2)
	: bf_(bf::make_hasher(num_hashes), num_elems*size_factor), k(k), bf_size(num_elems*size_factor), h(num_hashes)
	{};

	virtual ~BaseBloomFilter(){ //print out stats in destructor
		cerr << "# kmers inserted: " << num_inserted << endl;
		cerr << "BF size: " << bf_size << endl;
		cerr << "# hashes: " << h << endl;
	}

	void add(const kmer_t & kmer) {
		bf_.add(kmer);
		num_inserted++;
	}

	virtual bool contains(const kmer_t & kmer) {
		// cerr << "calling classic contains" << endl;
		return bf_.lookup(kmer);
	}

	virtual void populate(const unordered_set<kmer_t> & kmer_set){
		for(auto km : kmer_set){
			add(km);
		}
	}

};

#endif
