/*************************************************************
* Copyright (c) David Pellow, Darya Filippova, Carl Kingsford
*************************************************************/

#ifndef LIB_KBF_S
#define LIB_KBF_S

#include <bf.h>
#include "BaseBloomFilter.hpp"
#include "KBFUtil.hpp"

class KBFSparse : public BaseBloomFilter{
private:
  unordered_set<kmer_t> edge_kmers;
  const unsigned skip_len; //distance between stored kmers
  size_t extended_check = 0; // # found in bf, checked both sides
  size_t sparse_check = 0; // # not found in BF - check both sides
  size_t edge_check = 0; // # of checks in edge set
  size_t edge_pass = 0; // # of edge checks that pass

  // constructor action - avoid code duplication
  void constructKBF(unordered_set<kmer_t> & kmer_set,unordered_set<kmer_t> & edge_set){
    populate(kmer_set);
    getEdges(edge_set);
    cerr << "Edges: " << edge_kmers.size() << endl;
    sparse_check = 0;
    extended_check = 0;
    edge_check = 0;
    edge_pass = 0;
  }

public:

  KBFSparse(const int k, unordered_set<kmer_t> & kmer_set, unordered_set<kmer_t> & edge_set, const unsigned skip_len = 1, const size_t size_factor = 10) : BaseBloomFilter(k, kmer_set.size(), size_factor), skip_len(skip_len){
    constructKBF(kmer_set, edge_set);
  }

  //Constructor 2 - takes size explicitly
  KBFSparse(const size_t num_elems, const int k, unordered_set<kmer_t> & kmer_set, unordered_set<kmer_t> & edge_set, const unsigned skip_len = 1, const size_t size_factor = 10) : BaseBloomFilter(k, num_elems, size_factor), skip_len(skip_len){
    constructKBF(kmer_set,edge_set);
  }

  ~KBFSparse() {
    cerr << "# contained in BF: " << extended_check << endl;
    cerr << "# of sparse checks (not in BF): " << sparse_check << endl;
    cerr << "# edge set checks: " << edge_check << endl;
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
    if(bf_.lookup(kmer)){ // query is in the BF - skip and check neighbours
      extended_check++;
      bool left = skipAndSearchLeft(kmer, skip_len+1, k, &bf_);
      bool right = skipAndSearchRight(kmer, skip_len+1, k, &bf_);
      if (left && right)
        return true;
      if(left || right){
        edge_check++;
        if(edge_kmers.find(kmer)!=edge_kmers.end()){
          edge_pass++;
          return true;
        }
      }
    }

    // did not find kmer, check if the neighbouring ones are there
    // if neighbours exist the distance between the two edges should match
    // for example, if a sequence is made up of kmers:
    // A = a0 a1 a2 a3 a4 a5 and skip_len is 4, the kmers a0 and a5 will be saved in the filter
    // if a2 is being queried then we need to skip 1 when we search left and 2 when searching right
    // the sum of skip lengths on each side is restricted to be skip_len - 1 for a true match
    bool found_one = false; // could be an edge kmer if only find match on one side
    for(unsigned i = 0; i < skip_len; i++){
      sparse_check++;
      bool left = skipAndSearchLeft(kmer, i+1, k, &bf_);
      bool right = skipAndSearchRight(kmer,skip_len-i,k,&bf_);
      if(left && right)
        return true;
      found_one |= (left||right);
    }
    if(found_one){
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
