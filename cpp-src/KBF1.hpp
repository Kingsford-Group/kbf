/*************************************************************
* Copyright (c) David Pellow, Darya Filippova, Carl Kingsford
*************************************************************/

#ifndef LIB_KBF_1
#define LIB_KBF_1

#include <bf.h>
#include "BaseBloomFilter.hpp"
#include "KBFUtil.hpp"


class KBF1 : public BaseBloomFilter {

private:

  const unsigned extend_len; // length to extend up to and check for kmers

  size_t extended_check = 0; // number of extended checks for all kmers (each side counted)

public:

  KBF1(const int k, const size_t num_elems = 1024 * 1024, const unsigned extend_len = 1)
    : BaseBloomFilter(k, num_elems),
      extend_len(extend_len) {}

  //An alternative constructor that populates with kmer set
  KBF1(const int k, const unordered_set<kmer_t> & kmer_set, const unsigned extend_len=1)
    : BaseBloomFilter(k, kmer_set.size()),
      extend_len(extend_len){
    populate(kmer_set);
  }

  ~KBF1() {
    cerr << "Extended checks: " << extended_check  << endl;
  }

  bool contains(const kmer_t & A) {
    auto in_bf = bf_.lookup(A);
    // if A is not in BF -- return right away
    if (!in_bf) {
      return in_bf;
    }
    extended_check++;
    auto left_in_bf = searchLeft(A, extend_len, k, &bf_);
    if (in_bf && left_in_bf) {
      return true;
    }
    // here, even if kmer A is in the BF, may return false b.c. its extension is not in BF
    extended_check++;
    auto right_in_bf = searchRight(A, extend_len, k, &bf_);
    if (in_bf && right_in_bf) {
      return true;
    }
    return false;
  }

};

#endif
