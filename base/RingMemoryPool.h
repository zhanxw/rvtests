#ifndef _RINGMEMORYPOOL_H_
#define _RINGMEMORYPOOL_H_

#include <assert.h>
#include <string.h>
#include <vector>

class RingMemoryPool {
 public:
  RingMemoryPool() {
    capacity_ = 0;
    base_ = 0;
    head_ = 0;
    tail_ = 0;
    numElementsInChunk_ = 0;
  }
  RingMemoryPool(int numElementsInChunk, int numChunk = 64)
      : numElementsInChunk_(numElementsInChunk) {
    memory_.resize(numElementsInChunk * numChunk);
    capacity_ = numChunk;
    base_ = 0;
    head_ = 0;
    tail_ = 0;
  }
  int allocate() {
    if (tail_ == capacity_) {
      const size_t oldSize = capacity_;
      capacity_ = 2 * oldSize;
      memory_.resize(capacity_ * numElementsInChunk_);

      if (head_ != 0) {  // need to move data from head_ to capacity
        // move [head ... oldSize ] to [ (head+oldSize) ... capacity_];
        memcpy(memory_.data() + (head_ + oldSize) * numElementsInChunk_,
               memory_.data() + head_ * numElementsInChunk_,
               sizeof(float) * (oldSize - head_) * numElementsInChunk_);
        base_ += oldSize;
      }
    }
    const int idx = tail_;
    ++tail_;
    return idx;
  }
  void deallocate(int idx) {
    assert(idx == head_);  // only dealloc the head elment
    if (head_ < tail_) {
      ++head_;
    }
    assert(head_ <= tail_);
  }
  float* chunk(int idx) {
    if (idx >= tail_ || idx < head_) {
      return NULL;
    }
    const int pos = base_ + idx;
    float* ret = NULL;
    if (pos < capacity_) {
      ret = memory_.data() + pos * numElementsInChunk_;
    } else {
      ret = memory_.data() + (pos - capacity_) * numElementsInChunk_;
    }
    return ret;
  }

  // internal data may be not kept
  void setChunkSize(const int s) {
    assert(s >= 0);
    if (s == numElementsInChunk_) {
      return;
    }

    numElementsInChunk_ = s;
    if (memory_.size() == 0) {
      capacity_ = 64;
      memory_.resize(s * capacity_);
      base_ = head_ = tail_ = 0;
    } else if (memory_.size() % s == 0) {
      capacity_ = memory_.size() / s;
      base_ = head_ = tail_ = 0;
    } else {
      memory_.resize(s * (memory_.size() / s));
      capacity_ = memory_.size() / s;
      base_ = head_ = tail_ = 0;
    }
  }

 private:
  std::vector<float> memory_;  // hold some memory
  int capacity_;               // how many indice can be hold
  int base_;                   // the index for the first element
  int head_;                   // index of chunk
  int tail_;                   // index of chunk
  int numElementsInChunk_;
};

#endif /* _RINGMEMORYPOOL_H_ */
