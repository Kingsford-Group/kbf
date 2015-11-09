/*************************************************************
* Copyright (c) David Pellow, Darya Filippova, Carl Kingsford
*************************************************************/

#ifndef LIB_KBF_SR
#define LIB_KBF_SR

#include <bf.h>
#include "BaseBloomFilter.hpp"
#include "KBFUtil.hpp"

class KBFSparseRelaxed : public BaseBloomFilter {
private:
  unordered_set<kmer_t> edge_kmers;
  const unsigned skip_len; // distance between stored kmers
  size_t extended_check = 0; // found in bf - check both sides
  size_t sparse_check = 0; // not found in bf - check both sides
  size_t edge_check = 0; // checks of the edge set
  size_t edge_pass = 0; // edge checks that pass

  // constructor actions
  void constructKBF(unordered_set<kmer_t> & kmer_set, unordered_set<kmer_t> & edge_set){
    populate(kmer_set);
    getEdges(edge_set);
    cerr << "Edges: " << edge_kmers.size() << endl;
    sparse_check = 0;
    extended_check = 0;
    edge_check = 0;
    edge_pass = 0;
  }

public:

  KBFSparseRelaxed(const int k, unordered_set<kmer_t> & kmer_set, unordered_set<kmer_t> & edge_set, const unsigned skip_len = 1, const size_t size_factor = 10) : BaseBloomFilter(k, kmer_set.size(), size_factor), skip_len(skip_len){
    constructKBF(kmer_set, edge_set);
  }

  //Constructor 2 - takes size explicitly
  KBFSparseRelaxed(const size_t num_elems, const int k, unordered_set<kmer_t> & kmer_set, unordered_set<kmer_t> & edge_set, const unsigned skip_len = 1, const size_t size_factor = 10) : BaseBloomFilter(k,num_elems, size_factor), skip_len(skip_len){
    constructKBF(kmer_set, edge_set);
  }

  ~KBFSparseRelaxed(){
    cerr << "# contained in BF: " << extended_check << endl;
    cerr << "# of sparse checks (not in BF): " << sparse_check << endl;
    cerr << "# edge set checks: " << edge_check << endl;
    cerr << "# that pass: " << edge_pass << endl;
  }

  void getEdges(unordered_set<kmer_t> & edges) {
    for(auto km : edges){
      if(!contains(km))
	edge_kmers.insert(km);
    }
  }

  // Relaxed version of the contains() fnc from SeqBFSparse
  // just require the sum of the skip lengths on either side to be <= skip_len -1 (not exactly equal to it)
  bool contains(const kmer_t & kmer){
    if(!bf_.lookup(kmer)){
      // did not find kmer, check if the neighbouring ones are there
      bool left = false;
      bool right = false;
      for(unsigned i = 0; i < skip_len; i++){
        sparse_check++;
        left |= skipAndSearchLeft(kmer, i+1, k, &bf_);
        right |= skipAndSearchRight(kmer,skip_len-i,k,&bf_);
        if(left && right)
          return true;
      }
      if(left || right){
        edge_check++;
        if(edge_kmers.find(kmer)!=edge_kmers.end()){
          edge_pass++;
          return true;
        }
      }
      return false;
    }
    else{ // query is in the BF - check neighbours
      bool left = false;
      bool right = false;
      for(unsigned i = 0; i <= skip_len; i++){
        extended_check++;
        left |= skipAndSearchLeft(kmer, i+1, k, &bf_);
        right |= skipAndSearchRight(kmer,i+1,k,&bf_);
        if(left && right)
          return true;
      }
      if(left || right){
        edge_check++;
        if(edge_kmers.find(kmer)!=edge_kmers.end()){
          edge_pass++;
          return true;
        }
      }
      return false;
    }
  }

};

#endif
