#include "RingMemoryPool.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  // memcpy
#include <cassert>

RingMemoryPool::RingMemoryPool() { init(0, DefaultCapacity); }
RingMemoryPool::RingMemoryPool(int numElementsInChunk) {
  init(numElementsInChunk, DefaultCapacity);
}
/**
 * Allocate @param numChunk of chunks, each chunk has @param
 * numElementsInChunk elements.
 */
RingMemoryPool::RingMemoryPool(int numElementsInChunk, int numChunk) {
  init(numElementsInChunk, numChunk);
}
void RingMemoryPool::init(int numElementsInChunk, int numChunk) {
  assert(numElementsInChunk >= 0);
  assert(numChunk > 0);
  numElementsInChunk_ = numElementsInChunk;
  memory_.resize((size_t)numElementsInChunk * numChunk);
  capacity_ = numChunk;
  size_ = 0;
  head_ = 0;
  tail_ = 0;
  headIndex_ = 0;
  tailIndex_ = 0;
}
int RingMemoryPool::allocate() {
  // Case 1: . . Head x x x Tail . .  (head < tail)
  // Case 2: x x Tail . . . Head x x  (head > tail)
  // Case 3: . . Head/Tail . . . (head == tail)
  // Case 4: x x Head/Tail x x x (head == tail)
  // . => empty slot
  // x => used slot
  assert(numElementsInChunk_ > 0);
  // allocate memory if needed
  assert(capacity_ > 0);
  if (size_ == capacity_) {
    capacity_ = 2 * size_;
    memory_.resize(capacity_ * numElementsInChunk_);

    if (tail_ <= head_ && size_ > 0) {  // case 2 or 4
      // need to move data from head_ to oldCapacity(which is size_)
      // move [head ... oldCapacity ] to [ (head+oldCapacity) ... capacity_];
      memcpy(memory_.data() + (head_ + size_) * numElementsInChunk_,
             memory_.data() + head_ * numElementsInChunk_,
             sizeof(float) * (size_ - head_) * numElementsInChunk_);
      head_ += size_;
    }
  }

  ++tail_;
  if (tail_ == capacity_) {
    tail_ = 0;
  }
  const int ret = tailIndex_;
  tailIndex_++;
  size_++;
  return ret;
}
void RingMemoryPool::deallocate(int idx) {
#ifndef NDEBUG
  if (idx != headIndex_) {
    fprintf(stderr,
            "Cannot deallocate memory %d, headIndex = %d, tailIndex = %d\n",
            idx, headIndex_, tailIndex_);
    return;
  }
#endif
  assert(idx == headIndex_);  // only dealloc the head elment
  ++head_;
  if (head_ == capacity_) {
    head_ = 0;
  }
  headIndex_++;
  size_--;
}
float* RingMemoryPool::chunk(int idx) {
  if (idx >= tailIndex_ || idx < headIndex_) {
    return NULL;
  }
  const size_t pos = head_ + (idx - headIndex_);
  float* ret = NULL;
  if (pos < capacity_) {
    ret = memory_.data() + pos * numElementsInChunk_;
  } else {
    ret = memory_.data() + (pos - capacity_) * numElementsInChunk_;
  }
  return ret;
}
float* RingMemoryPool::firstChunk() { return memory_.data(); }
float* RingMemoryPool::lastChunk() {
  return memory_.data() + (capacity_ - 1) * numElementsInChunk_;
}
size_t RingMemoryPool::capacity() const { return capacity_; }
size_t RingMemoryPool::size() const { return size_; }
// internal data will be reset
void RingMemoryPool::setChunkSize(const size_t s) {
  if (s == numElementsInChunk_) {
    return;
  }

  numElementsInChunk_ = s;
  if (memory_.size() == 0) {
    capacity_ = DefaultCapacity;
    memory_.resize(s * capacity_);
  } else {
    memory_.resize(s * (memory_.size() / s));
    capacity_ = memory_.size() / s;
  }
  head_ = tail_ = headIndex_ = tailIndex_ = 0;
}
