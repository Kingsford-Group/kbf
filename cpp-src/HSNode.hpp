/*************************************************************
* Copyright (c) David Pellow, Darya Filippova, Carl Kingsford
*************************************************************/

#ifndef LIB_HS_NODE
#define LIB_HS_NODE

#include <unordered_set>
#include <boost/heap/fibonacci_heap.hpp>

typedef uint64_t kmer_t;

class HSNode {
private:
  const kmer_t kmer;
  unsigned degree = 0;
  unordered_set<unordered_set<kmer_t>*> neighbours;

public:

  HSNode(kmer_t kmer) : kmer(kmer){
  }

  kmer_t get_kmer(){
    return kmer;
  }

  unordered_set<unordered_set<kmer_t>*> get_neighbours(){
    return neighbours;
  }

  unsigned get_degree(){
    return degree;
  }

  void add_neighbour(unordered_set<kmer_t>* N){
    neighbours.insert(N);
    degree++;
  }

  void remove_neighbour(unordered_set<kmer_t>* N){
    neighbours.erase(N);
    degree--;
  }

};

#endif
