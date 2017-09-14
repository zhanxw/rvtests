#ifndef _RINGMEMORYPOOL_H_
#define _RINGMEMORYPOOL_H_

#include <cstddef>  // size_t definnition
#include <vector>

/**
 * This class uses a ring to store indice,
 * and each index points to a chunk of memory of a fixed length.
 * Limitation:
 * when allocate memory, you get a sereis of index: 0, 1, 2, 3...
 * wehn deallocate memory, you will need to deallocate 0, 1, 2, 3, ... in the
 * same order
 */
class RingMemoryPool {
 public:
  RingMemoryPool();
  RingMemoryPool(int numElementsInChunk);
  /**
   * Allocate @param numChunk of chunks, each chunk has @param
   * numElementsInChunk elements.
   */
  RingMemoryPool(int numElementsInChunk, int numChunk);
  void init(int numElementsInChunk, int numChunk);

  int allocate();
  void deallocate(int idx);
  float* chunk(int idx);
  float* firstChunk();
  float* lastChunk();
  size_t capacity() const;
  size_t size() const;
  // internal data will be reset
  void setChunkSize(const size_t s);

 private:
  static const int DefaultCapacity = 64;
  std::vector<float> memory_;  // hold chunks memory
  size_t capacity_;            // how many indice can be hold
  size_t size_;                // number of existing elements
  size_t head_;                // pointer to the ring head
  size_t tail_;                // pointer to the ring tail
  int headIndex_;              // the first unallocated index for the memory
  int tailIndex_;              // the index to be allocated
  size_t numElementsInChunk_;
};

#endif /* _RINGMEMORYPOOL_H_ */
