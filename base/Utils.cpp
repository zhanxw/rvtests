#include "Utils.h"

#ifdef __SSE2__
#pragma message "Enable SSE2 optimized ssechr"
// copy from:
// https://mischasan.wordpress.com/2011/06/22/what-the-is-sse2-good-for-char-search-in-long-strings/

#include <emmintrin.h>
char const* ssechr(char const* s, char ch) {
  __m128i zero = _mm_setzero_si128();  // set zero 16 times
  __m128i cx16 = _mm_set1_epi8(ch);    // (ch) replicated 16 times.
  while (1) {
    // load 128 bit, @param s does not need to be aligned
    // on little endian system, s[0] is the least significant part
    // on memory, it looks like s[15], s[14], ..., s[0]
    __m128i x = _mm_loadu_si128((__m128i const*)s);
    // _mm_cmpeq_epi8: compare 16 times on the 8 bit number, set ff for equal or
    // 00 for not equal
    // _mm_movemask_epi8: extract each of the highest significnat bit of the 16
    // 8bit number
    // this command identify the location of '\0'
    unsigned u = _mm_movemask_epi8(_mm_cmpeq_epi8(zero, x));
    //  ~u, change all 1 to 0, all 0 to 1, e.g. 0110 0000 -> 1001 1111
    // (u-1), least significant 0s and 1 will be flipped, 0110 0000 -> 0101 1111
    // ~u & (u-1), the large significant part -> 0, the least significat 1 and
    // all traiting zeros -> 1
    // e.g. 0110 0000 -> 0001 1111
    // _mm_movemask_epi8(_mm_cmpeq_epi8(cx16, x)): will set the bit where s[i]
    // == ch to 1
    unsigned v = _mm_movemask_epi8(_mm_cmpeq_epi8(cx16, x)) & ~u & (u - 1);
    // ffs find the first bit in a word, e.g. 0110 0000 -> 6
    if (v) return s + ffs(v) - 1;  //
    if (u) return NULL;            // does not find char
    s += 16;
  }
}

#else
#pragma message "Disabled SSE2 => no optimized ssechr"
#define ssechr strchr
#endif
