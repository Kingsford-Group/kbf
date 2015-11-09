/*************************************************************
* Copyright (c) David Pellow, Darya Filippova, Carl Kingsford
*************************************************************/

#ifndef LIB_KBF_2
#define LIB_KBF_2

#include <bf.h>
#include "BaseBloomFilter.hpp"
#include "KBFUtil.hpp"

class KBF2 : public BaseBloomFilter {
private:
  const unsigned extend_len;
  unordered_set<kmer_t> edge_kmers; //set of kmers on sequence edges - will fail BF contains()
  size_t extended_check = 0; // # of times checks extensions
  size_t edge_check = 0; // # of times checks edge set
  size_t edge_pass = 0; // # of times returns true from edge check
public:
  // constructor 1:
  // takes in kmer sets - kmers and potential edges (read edges)
  KBF2(const int k, unordered_set<kmer_t> & kmer_set, unordered_set<kmer_t> & edge_set, const unsigned extend_len = 1, const size_t size_factor = 10) : BaseBloomFilter(k, kmer_set.size(), size_factor), extend_len(extend_len){
    populate(kmer_set);
    getEdges(edge_set);
    cerr << "Edges: " << edge_kmers.size() << endl;
    extended_check = 0;
    edge_check = 0;
    edge_pass = 0;
  }

  ~KBF2(){
    cerr << "# of extended checks: " << extended_check << endl;
    cerr << "# of edge checks: " << edge_check << endl;
    cerr << "# that pass: " << edge_pass << endl;
  }

  // any read edge kmer that is not in the filter needs
  // to be stored in edge_kmers
  void getEdges(unordered_set<kmer_t>& edges){
    for(auto km : edges){
      if(!contains(km))
        edge_kmers.insert(km);
    }
  }

  bool contains(const kmer_t & kmer){
    if(!bf_.lookup(kmer)){
      return false;
    }
    extended_check++;
    bool containsLeft = searchLeft(kmer,extend_len,k,&bf_);
    bool containsRight = searchRight(kmer,extend_len,k,&bf_);
    if(containsLeft && containsRight){
      return true;
    }
    if(containsLeft || containsRight){
      //assumes read length > k+extend_len
      edge_check++;
      if(edge_kmers.find(kmer)!=edge_kmers.end()){
        edge_pass++;
        return true;
      }
    }
    return false;
  }

};

#endif
