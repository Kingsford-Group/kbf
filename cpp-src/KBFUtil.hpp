/*************************************************************
* Copyright (c) David Pellow, Darya Filippova, Carl Kingsford
*************************************************************/

#ifndef KBF_UTIL
#define KBF_UTIL

#include <unordered_set>
#include <unordered_map>
#include "BaseBloomFilter.hpp"
#include "HSNode.hpp"
#include "FastaReader.h"
#include "JellyfishUtil.h"

bool contains(unordered_set<kmer_t> & set, kmer_t t) {
  return set.find(t) != set.end();
}

////////////////////////////////////////////////
// search for kmer extensions to the left of the input
// kmer in bf up to extensions of length len
// returns true if it found at least one kmer within LEN steps to the left
//  -----
//  -----
//   -----
////////////////////////////////////////////////
bool searchLeft(const kmer_t & kmer, unsigned len, const int k, bf::basic_bloom_filter * bf){
    if (len == 0) return false;
    bool contains = false;
    kmer_t suffix = kmer >> 2;
    int shift = 2 * (k - 1);
    for (kmer_t i = 0; i < 4; i++){
      kmer_t extension = suffix | (i << shift);
      // if this extension is not in BF -- keep looking
      if (bf->lookup(extension) )
        return true;
      else
          contains |= searchLeft(extension, len - 1, k, bf);
      if (contains)
        return contains;
    }
    return contains;
}


////////////////////////////////////////////////
// search for kmers that are extensions of the input
// of length dist
////////////////////////////////////////////////
bool skipAndSearchLeft(const kmer_t& kmer, unsigned dist, const int k, bf::basic_bloom_filter * bf){
  kmer_t suffix = kmer >> (2*dist); //TODO: put in check that this distance is not too far
  for(kmer_t i = 0; i < pow(4,dist); i++){
    kmer_t extension = suffix | (i<<((k-dist)*2));
    if(bf->lookup(extension))
      return true;
  }
  return false;
}

//same for the right
bool searchRight(const kmer_t & kmer, unsigned len, const int k, bf::basic_bloom_filter * bf){
  if (len == 0) return false;
  bool contains = false;
  kmer_t mask = ( ((kmer_t)1) << (2*(k-1) ) ) -1;
  kmer_t prefix = (kmer & mask) << 2;
  for (kmer_t i = 0; i < 4; i++) {
    kmer_t extension = prefix | i;
    if (bf->lookup(extension))
      return true;
    else
      contains |= searchRight(extension, len-1, k, bf);
    if(contains)
      return contains;
  }
  return contains;
}


////////////////////////////////////////////////
bool skipAndSearchRight(const kmer_t& kmer, unsigned dist, const int k, bf::basic_bloom_filter * bf){
//  kmer_t prefix = kmer<<(2*dist);
//  prefix &= (1<<(k*2))-1;
  kmer_t mask = ( ((kmer_t)1) << (2*(k-dist) ) ) -1;
  kmer_t prefix = (kmer & mask) << (2*dist);
  for(kmer_t i = 0; i < pow(4,dist); i++){
    kmer_t extension = prefix | i;
    if(bf->lookup(extension))
      return true;
  }
  return false;
}

////////////////////////////////////////////////
// write out kmers to binary file
// returns false if fopen or fwrite errors
////////////////////////////////////////////////
bool writeKmers(const unordered_set<kmer_t> & kmer_set, const string & fname){
  FILE* outfile;
  size_t buf_size = 1024*1024;
  outfile = fopen((fname+".bin").c_str(),"wb");
  if(!outfile) {
    cerr << "[ERROR] Could not open file for writing" << endl;
    return false; //better error handling?
  }
  size_t i = 0;
  size_t written;
  kmer_t * data_buf = new kmer_t[buf_size];
  for(auto km: kmer_set){
    data_buf[i] = km;
    //filled buffer
    if(i==buf_size-1){
      written = fwrite(data_buf,sizeof(kmer_t),buf_size,outfile);
      if(written != buf_size){ //better error handling
        fclose(outfile);
        delete data_buf;
        return false;
      }
      i = 0;
      continue;
    }
    i++;
  }
  written = fwrite(data_buf,sizeof(kmer_t),i,outfile);
  fclose(outfile);
  delete data_buf;
  return (written == i);
}

////////////////////////////////////////////////
//read in a kmer binary file and convert to kmer set
////////////////////////////////////////////////
bool readKmerFile(unordered_set<kmer_t> & kmer_set, const string & fname){
  FILE* infile;
  size_t buf_size = 1024*1024;
  infile = fopen(fname.c_str(),"rb");
  if(!infile) return false;
  kmer_t * data_buf = new kmer_t[buf_size];
  //get the number of elements
  fseek(infile,0,SEEK_END);
  long long fSize = ftell(infile);
  rewind(infile);
  long long num_kmers = fSize/sizeof(kmer_t);

  size_t numRead;
  while(num_kmers > 0){
    numRead = fread(data_buf,sizeof(kmer_t),buf_size,infile);
    if(numRead < buf_size && numRead < num_kmers){
      delete data_buf;
      fclose(infile);
      return false;
    }
    for(size_t i = 0; i < numRead; i++){
      kmer_set.insert(data_buf[i]);
    }
    num_kmers -= numRead;
  }
  delete data_buf;
  return true;
}

////////////////////////////////////////////////////////////////
// parse a fasta file and return vector of reads
////////////////////////////////////////////////////////////////
vector<string> parseFasta(string const & path) {
    vector<string> reads;
    FastaReader fr(path.c_str());
    kseq_t * seq;
    size_t cnt = 0;
    while ( (seq = fr.nextSequence() ) ) {
      // cerr << seq->seq.s << endl;
      reads.push_back(seq->seq.s);
      cnt++;
    }
    cerr << "(" << cnt << " reads) ";
    return reads;
}

////////////////////////////////////////////////////////////////
// take in reads and output the kmer set
////////////////////////////////////////////////////////////////
unordered_set<kmer_t> getKmers(vector<string> & reads, const int K) {
    unordered_set<kmer_t> kmers;
    for (auto r : reads) {

        if (r.size() < K) continue;
        for (int i = 0; i < r.size() - K + 1; i++) {
            // if (i+K > r.size()) {
            //   cerr << i + K << endl;
            // }
            kmer_t kmer_bin = mer_string_to_binary(&r[i], K);
            // kmers.insert( r.substr(i, K) );
            kmers.insert(kmer_bin);
        }
    }
    // cerr << endl;
    return kmers;
}

// take in reads and output both the kmer set and the edge kmer set
// edge kmers are the kmers that are within extend_len of either end of a read
void getKmersAndEdgeKmers(vector<string> & reads, const int K, const unsigned extend_len, unordered_set<kmer_t> & kmers, unordered_set<kmer_t> & edgeKmers){
  for (auto r : reads) {
    if (r.size() < K) continue;
    for (size_t i = 0; i < r.size() - K + 1; i++) {
      kmer_t kmer_bin = mer_string_to_binary(&r[i], K);
      kmers.insert(kmer_bin);
      if ( i < extend_len || i > r.size()-K-extend_len){
        edgeKmers.insert(kmer_bin);
      }
    }
  }
}

// same, but only taking every (N+1)th kmer
// ** this is to be used for sparse SBFs where the inputs
// are full length sequences, not different reads that could overlap **
void getNthKmersAndEdgeKmers(vector<string> & sequences, const int K, const unsigned skip_len, unordered_set<kmer_t> & kmers, unordered_set<kmer_t> & edgeKmers){
  for (auto seq : sequences) {
    if (seq.size() < K) continue;
    size_t last_inserted = 0;
    for (size_t i = 0; i < seq.size() - K + 1; i++) {
      kmer_t kmer_bin = mer_string_to_binary(&seq[i], K);
      if(!(i%(skip_len+1))){ // N+1th kmer
        kmers.insert(kmer_bin);
        last_inserted=i;
      }
      if (i == 0)
        edgeKmers.insert(kmer_bin);
    }
    for(size_t i=last_inserted; i<seq.size()-K+1;i++)
      edgeKmers.insert(mer_string_to_binary(&seq[i], K));
  }
}

// take in reads and output the sparsified kmer set and the edge kmers from the reads
// for each read check which index to start sparsifying from based on how many
// kmers would match those that have laready been taken
void getBestFitKmersAndEdgeKmers(vector<string> & sequences, const int K, const unsigned skip_len, unordered_set<kmer_t> & kmers, unordered_set<kmer_t> & edgeKmers){
  for(auto seq: sequences) {
    if (seq.size() < K) continue;
    // get all the kmers and count matches to current set
    vector<kmer_t> seqKmers(seq.size()-K+1);
    vector<int> start_inds(skip_len+1,0);
    for (size_t i = 0; i < seq.size() - K + 1; i++) {
      seqKmers[i] = mer_string_to_binary(&seq[i],K);
      if(kmers.find(seqKmers[i]) != kmers.end()){ // this kmer has been seen before
        start_inds[i%(skip_len+1)]++;
      }
    }
    int max_ind = 0; // find the best matching start index
    int max_val = start_inds[0];
    for(size_t i=1; i<start_inds.size();i++){
      if(start_inds[i] > max_val){
        max_ind = i;
        max_val = start_inds[i];
      }
    }
    // get the sparse kmer set starting from max_ind and skipping by skip_len
    size_t last_inserted = 0;
    for (size_t i = 0; i < seq.size() - K + 1; i++){
      if(i <= max_ind) // at the beginning, before we start taking kmers
        edgeKmers.insert(seqKmers[i]);

      if(!((i-max_ind)%(skip_len+1))){
        kmers.insert(seqKmers[i]);
        last_inserted = i;
      }
    }
    for(size_t i=last_inserted; i<seq.size()-K+1;i++)
      edgeKmers.insert(seqKmers[i]);
  }
}

//for node comparison in heap
struct node_compare {
  bool operator()(HSNode * node1, HSNode * node2) const{
    return node1->get_degree() < node2->get_degree();
  }
};

// greedy hitting set sparsification
// create a node for each kmer with the kmer value, and pointers to the neighbour sets that it covers
// greedily add the highest degree kmer to the sparse kmer set
// then remove those nodes and their neighbours from the heap and decrease the degree of their neighbours
void hittingSetKmersAndEdgeKmers(vector<string> & sequences, const int K, unordered_set<kmer_t> & kmers, unordered_set<kmer_t> & edgeKmers){
  // sets of all neighbours for each kmer
  unordered_map<kmer_t, unordered_set<kmer_t> > Left;
  unordered_map<kmer_t, unordered_set<kmer_t> > Right;
  // get the kmers for the sequences
  for (auto seq : sequences) {
    if (seq.size() < K) continue;
    vector<kmer_t> seqKmers(seq.size()-K+1);
    for (size_t i = 0; i < seq.size() - K + 1; i++) {
      kmer_t kmer_bin = mer_string_to_binary(&seq[i], K);
      if(i == 0 || i == seq.size()-K){
        edgeKmers.insert(kmer_bin);
      }
      seqKmers[i] = kmer_bin;
    }
    // insert kmers into sets
    for(size_t i = 0; i < seq.size() - K + 1; i++){
      if(i != 0)
        Left[seqKmers[i]].insert(seqKmers[i-1]);
      if(i != seq.size()-K)
        Right[seqKmers[i]].insert(seqKmers[i+1]);
      Left[seqKmers[i]].insert(seqKmers[i]);
      Right[seqKmers[i]].insert(seqKmers[i]);
    }
  }
  // create the heap of nodes, map each kmer to its handle in the heap, get all the neighbours
  unordered_map<kmer_t, boost::heap::fibonacci_heap<HSNode*, boost::heap::compare<node_compare> >::handle_type> kmer_handles;
  boost::heap::fibonacci_heap<HSNode*, boost::heap::compare<node_compare> > HSNode_heap;
  for(auto& km_set: Left){
    for(auto& km: km_set.second){
      if(kmer_handles.find(km)==kmer_handles.end()){ //new node
        HSNode* kmer_node = new HSNode(km);
        kmer_handles[km] = HSNode_heap.push(kmer_node);
      }
      (*kmer_handles[km])->add_neighbour(&(km_set.second));
      HSNode_heap.update(kmer_handles[km]);
    }
  }
  for(auto& km_set: Right){
    for(auto& km: km_set.second){
      if(kmer_handles.find(km)==kmer_handles.end()){ //new node
        HSNode* kmer_node = new HSNode(km);
        kmer_handles[km] = HSNode_heap.push(kmer_node);
      }
      (*kmer_handles[km])->add_neighbour(&(km_set.second));
      HSNode_heap.update(kmer_handles[km]);
    }
  }

  //TODO: Erase handles from kmer_handles when they won't be needed anymore
  //      Erase kmer sets from Left and Right once they are covered

  while(!HSNode_heap.empty()){
    // while the heap is not empty get the highest degree node
    // put it in the kmer set, remove it and decrement the degrees of all the nodes that hit the same neighbour sets
    HSNode * max_node = HSNode_heap.top();
    kmers.insert(max_node->get_kmer());
    HSNode_heap.pop();
    kmer_handles.erase(max_node->get_kmer());
    vector<HSNode *> nodes_to_remove;
    vector<kmer_t> left_sets_to_remove;
    vector<kmer_t> right_sets_to_remove;
    for(auto& N : max_node->get_neighbours()){ //for every neighbour set
      for(auto& km: *N){ // for every kmer in this set
        //update the nodes pointing to it
        if(km == max_node->get_kmer()) continue;
        (*kmer_handles[km])->remove_neighbour(N);
        HSNode_heap.update(kmer_handles[km]);
        if((*kmer_handles[km])->get_degree() <= 0){
          HSNode_heap.erase(kmer_handles[km]);
          nodes_to_remove.push_back((*kmer_handles[km]));
          kmer_handles.erase(km);
        }
      }
    }
    while (nodes_to_remove.size()) {
      auto n = nodes_to_remove.back();
      nodes_to_remove.pop_back();
      delete n;
    }
    delete max_node;
  }
}

#endif
