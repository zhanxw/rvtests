#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "base/RingMemoryPool.h"

int main() {
  float* p;

  {
    RingMemoryPool mp(100, 64);
    for (int i = 0; i < 32; ++i) {
      int idx = mp.allocate();
      p = mp.chunk(idx);
      for (int j = 0; j < 100; ++j) {
        p[j] = j;
      }
    }
    for (int i = 0; i < 32; ++i) {
      mp.deallocate(i);
    }
  }

  {
    RingMemoryPool mp(10, 64);
    for (int i = 0; i < 128; ++i) {
      int idx = mp.allocate();
      p = mp.chunk(idx);
      for (int j = 0; j < 10; ++j) {
        p[j] = j + 1;
      }
    }
    // expect sum = 128 * (1+2+ ...+10)
    float sum = 0;
    for (int i = 0; i < 128; ++i) {
      p = mp.chunk(i);
      float subSum = 0;
      for (int j = 0; j < 10; ++j) {
        subSum += p[j];
      }
      assert((int)subSum == 55);
      sum += subSum;
    }
    assert(int(sum) == 128 * 55);

    for (int i = 0; i < 128; ++i) {
      mp.deallocate(i);
    }
  }

  {
    RingMemoryPool mp(4, 2);
    for (int i = 0; i < 8; ++i) {
      int idx = mp.allocate();
      p = mp.chunk(idx);
      for (int j = 0; j < 4; ++j) {
        p[j] = rand() % 1024;
      }
    }
    for (int i = 0; i < 8; ++i) {
      mp.deallocate(i);
    }

    float expectedSum = 0.0;
    for (int i = 0; i < 8; ++i) {
      int idx = mp.allocate();
      assert(idx == 8 + i);
      p = mp.chunk(idx);
      for (int j = 0; j < 4; ++j) {
        p[j] = rand() % 1024;
        expectedSum += p[j];
        // printf("set chunk = %d, offset = %d, value = %g, expectedSum = %g
        // \n",
        //        idx, j, p[j], expectedSum);
      }
    }

    // expect sum = 128 * (1+2+ ...+10)
    float sum = 0;
    for (int i = 0; i < 8; ++i) {
      p = mp.chunk(i);
      assert(p == NULL);
    }
    for (int i = 8; i < 8 + 8; ++i) {  // since we have allocated 8 times before
      p = mp.chunk(i);
      for (int j = 0; j < 4; ++j) {
        sum += p[j];
        // printf("get chunk = %d, offset = %d, value = %g, sum = %g\n", i, j,
        //        p[j], sum);
      }
    }
    assert(sum == expectedSum);

    for (int i = 0; i < 8; ++i) {
      mp.deallocate(8 + i);
    }
  }

  {
    RingMemoryPool mp(1, 4);
    for (int i = 0; i < 3; ++i) {
      int idx = mp.allocate();
      p = mp.chunk(idx);
      *p = i;
    }
    assert(mp.size() == 3);  // mp: 0 1 2 .

    mp.deallocate(0);
    mp.deallocate(1);
    assert(mp.size() == 1);  // mp: . . 2 .

    assert(mp.allocate() == 3);  // mp: . . 2 3
    *(mp.chunk(3)) = 3;
    assert(*(mp.lastChunk()) == 3);

    assert(mp.allocate() == 4);  // mp: 4 . 2 3
    *(mp.chunk(4)) = 4;
    assert(*(mp.firstChunk()) == 4);

    assert(mp.allocate() == 5);  // mp: 4 5 2 3
    *(mp.chunk(5)) = 5;
    assert(mp.size() == 4);
    assert(mp.capacity() == 4);

    mp.deallocate(2);
    mp.deallocate(3);
    assert(mp.allocate() == 6);
    assert(mp.allocate() == 7);  // mp: 4 5 6 7
    *(mp.chunk(6)) = 6;
    *(mp.chunk(7)) = 7;
    assert(*(mp.firstChunk()) == 4);
    assert(*(mp.lastChunk()) == 7);
    assert(mp.size() == 4);
    assert(mp.capacity() == 4);
  }
}
