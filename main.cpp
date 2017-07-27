#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <inttypes.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include <vector>
#include <climits>
#include <algorithm>

#include "wrap_traits.hpp"
#include "tables.hpp"
#include "tables256.hpp"

using std::vector;

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) < (b) ? (b) : (a))

#define VERBOSE 0
#define MHZ 2600
#define NUM_PRIME 30
#define CALLBACK &dummy_callback
#define BLOCK_COUNT 32 // how many blocks of the underlying READ_SIZE to process before running callbacks, you can't change this w/o fixing up the asm

typedef uint64_t ULONGLONG;
typedef unsigned char BYTE;

typedef void (callback_fn)(uint64_t, void*);
typedef struct {
    callback_fn* f;
    void* param;
} callback;
typedef uint64_t (process_fn)(uint64_t, uint64_t , const BYTE *, callback *);

extern "C" process_fn ProcessA_asm;
extern "C" process_fn ProcessA2_asm;

typedef void (kfunc)(uint8_t *offset_indexes, uint8_t *block_result);
extern "C" kfunc kernel64_asm;
extern "C" kfunc kernel256_asm;
extern "C" kfunc kernel512_asm;


// first 30 primes, minus one
const BYTE primes[] = {
        2,   4,   6,  10,  12,  16,  18,  22,
        28,  30,  36,  40,  42,  46,  52,  58,
        60,  66,  70,  72,  78,  82,  88,  96,
        100, 102, 106, 108, 112, 126,   0,   0 };

// a callback that just prints the prime
void print_callbackf(uint64_t prime, void *param) {
    printf("Found prime: %lu\n", prime);
}
callback print_callback = { print_callbackf, NULL };

void dummy_callbackf(uint64_t prime, void *param) {}
callback dummy_callback = { dummy_callbackf, NULL };

typedef struct {
    callback c;
    uint64_t *array;
    size_t size, capacity;
    size_t idx;
} recording_callback;

void record(uint64_t prime, void *param) {
    recording_callback *c = (recording_callback *)param;
    size_t sz = c->size;
    if (sz + 1 > c->capacity) {
        size_t newcap = c->capacity * 2 + 16;
        c->array = (uint64_t*)realloc(c->array, newcap * sizeof(uint64_t));
        c->capacity = newcap;
        //        printf("expand nc=%lu\n", newcap);
    }
    c->array[c->size++] = prime;
    //    printf("i=%lu p=%lu\n", c->idx++, prime);
}

void record_init(recording_callback *c) {
    *c = (recording_callback){ { &record, c }, NULL, 0, 0 };
}

/*
 * Subtract one from each byte.
 */
void sub1(BYTE *state) {
    for (int i = 0; i < NUM_PRIME; i++) {
        state[i]--;
    }
}

bool any_wrap(const BYTE *state) {
    for (int i = 0; i < NUM_PRIME; i++) {
        if (state[i] & 0x80) {
            return true;
        }
    }
    return false;
}

BYTE single_state(uint64_t start, BYTE prime) {
    // really dumb, keep subtracting 2 until you hit a multiple of the prime
    int distance = 0;
    start -= 2;
    while (start % prime != 0) {
        distance++;
        start -= 2;
    }
    assert(distance <= 255);
    return (BYTE)distance;
}

BYTE distance(uint64_t a, uint8_t b) {
    assert((a & 1) && (b && 1));  // a and b both odd
    assert(a > 2 && b != 0);
    // really dumb, keep subtracting 2 until you hit a multiple of the b
    uint64_t distance = 0;
    a -= 2;
    assert(a >= b);
    while (a % b != 0) {
        distance++;
        a -= 2;
    }
    assert(distance <= 255);
    return (uint8_t)distance;
}

void make_state_table(uint64_t start, BYTE *state) {
    assert(start & 1); // must be odd
    for (int i = 0; i < NUM_PRIME; i++) {
        uint64_t prime = primes[i] + 1;
        state[i] = single_state(start, prime);
    }
}

void print_state(uint64_t start) {
    BYTE state[NUM_PRIME];
    make_state_table(start, state);
    printf("BYTE state[32] = {\n\t");
    for (int i = 0; i < NUM_PRIME; i++) {
        printf("0x%02x, ", state[i]);
        if (i % 8 == 7)
            printf("\n\t");
    }
    printf("0x00, 0x00\n}\n");
}

uint64_t ProcessA(uint64_t start, uint64_t end, const BYTE *initial, callback* fn) {
    assert(NUM_PRIME == 30); // doens't work with other NUM_PRIME values
    return ProcessA_asm(start, end, initial, fn);
}

uint64_t ProcessA2(uint64_t start, uint64_t end, const BYTE *initial, callback* fn) {
    assert(NUM_PRIME == 30); // doens't work with other NUM_PRIME values
    return ProcessA2_asm(start, end, initial, fn);
}

uint64_t ProcessC(uint64_t start, uint64_t end, const BYTE *initial, callback* fn) {
    uint64_t count = 0;
    uint64_t cur = start;
    BYTE state[NUM_PRIME];
    memcpy(state, initial, NUM_PRIME);

    for(;;) {

        uint32_t loop_count = 0;

        for (;;) {

            sub1(state);

            if (!any_wrap(state))
                break; // found a prime

            loop_count++;

            // any counters that have wrapped to 0xFF, set them back to their prime-1
            for (int i = 0; i < NUM_PRIME; i++) {
                BYTE val   = state[i];
                BYTE prime = primes[i];
                state[i] = (val > prime ? prime : val);
            }

        }

        uint32_t delta = loop_count * 2 + 2;

        cur -= delta;

        if (cur < end || (int64_t)cur < 0) {
            break;
        }

        if (fn) fn->f(cur, fn->param);

        // found a prime
        count++;

    }

    return count;
}

static const uint64_t BITMAPS[] = {
        /*   3 */  0x9249249249249249 ,
        /*   5 */  0x1084210842108421 ,
        /*   7 */  0x8102040810204081 ,
        /*  11 */  0x0080100200400801 ,
        /*  13 */  0x0010008004002001 ,
        /*  17 */  0x0008000400020001 ,
        /*  19 */  0x0200004000080001 ,
        /*  23 */  0x0000400000800001 ,
        /*  29 */  0x0400000020000001 ,
        /*  31 */  0x4000000080000001 ,
        /*  37 */  0x0000002000000001 ,
        /*  41 */  0x0000020000000001 ,
        /*  43 */  0x0000080000000001 ,
        /*  47 */  0x0000800000000001 ,
        /*  53 */  0x0020000000000001 ,
        /*  59 */  0x0800000000000001 ,
        /*  61 */  0x2000000000000001 ,
        /*  67 */  0x0000000000000001 ,
        /*  71 */  0x0000000000000001 ,
        /*  73 */  0x0000000000000001 ,
        /*  79 */  0x0000000000000001 ,
        /*  83 */  0x0000000000000001 ,
        /*  89 */  0x0000000000000001 ,
        /*  97 */  0x0000000000000001 ,
        /* 101 */  0x0000000000000001 ,
        /* 103 */  0x0000000000000001 ,
        /* 107 */  0x0000000000000001 ,
        /* 109 */  0x0000000000000001 ,
        /* 113 */  0x0000000000000001 ,
        /* 127 */  0x0000000000000001 ,
};

/* interpret a bitmap starting at cur and issue callbacks for the likely prime */
void do_bitmap_callback(uint64_t cur, uint64_t bitmap, callback* fn) {
    assert(fn != NULL);
    for (int bit = 0; bit < 64; bit++, cur -= 2) {
        if ((bitmap & (1UL << bit)) == 0) {
            fn->f(cur, fn->param);
        }
    }
}

BYTE ss1(uint64_t a, uint8_t b) {
    a -= 2;
    if ((a % b) % 2 == 0) {
        return (a % b) / 2;
    } else {
        return (a % b + b) / 2;
    }
}


void test_state() {
    for (uint8_t b = 1; b <= 255; b += 2) {
        for (uint64_t a = (uint64_t)b + 2; a < 100001; a += 2) {
            uint8_t r1 = single_state(a, b);
            uint8_t r2 = ss1(a, b);
            if (r1 != r2) {
                printf("Failed for a=%d, b=%d, r1=%d, r2=%d", (int)a, (int)b, r1, r2);
                exit(1);
            }
        }
    }
    printf("OK!");
}

inline uint64_t shlz(uint64_t x, uint32_t shift) {
    return shift >= 64 ? 0 : (x << shift);
}

inline uint64_t read64(const uint8_t *array, size_t i) {
    uint64_t v;
    memcpy(&v, array + i, 8);
    return v;
}

inline void write64(uint64_t value, uint8_t *array, size_t i) {
    memcpy(array + i, &value, 8);
}

int32_t find_offset_for_shift(uint8_t primeidx, uint8_t shift, const uint8_t *bytemap, size_t bytemap_size) {

    uint8_t shiftedmap[8];
    write64(shlz(BITMAPS[primeidx], shift), shiftedmap, 0);

    // essentially we do a search to find the offset with the bitmap to find the byte at
    // which the next 8 bytes match the shifted bitmap. No doubt this could be calculated
    // directly with some modular arithmetic, but bleh
    for (int32_t offset = 0; offset + 7U < bytemap_size; offset++) {
        if (memcmp(shiftedmap, bytemap + offset, 8) == 0) {
            //printf("p=%d found match for shift %d at offset %d", prime, shift, offset);
            return offset;
        }
    }
    return -1;
}

vector<uint8_t> to_vector(const vector<bool>& input) {
    assert(input.size() % 8 == 0);

    vector<uint8_t> output;
    for (size_t i = 0; i < input.size(); i += 8) {
        output.push_back(
                (input[i + 0] << 0) |
                (input[i + 1] << 1) |
                (input[i + 2] << 2) |
                (input[i + 3] << 3) |

                (input[i + 4] << 4) |
                (input[i + 5] << 5) |
                (input[i + 6] << 6) |
                (input[i + 7] << 7)
        );
    }
    return output;
}

/* covert a byte vector to a bit-vector of bool starting with the LSB of the first byte */
vector<bool> to_vector(const vector<uint8_t>& input) {
    vector<bool> output;
    output.reserve(input.size() * 8);
    for (uint8_t b : input) {
        output.push_back(b & 1);
        output.push_back(b & 2);
        output.push_back(b & 4);
        output.push_back(b & 8);

        output.push_back(b & 16);
        output.push_back(b & 32);
        output.push_back(b & 64);
        output.push_back(b & 128);
    }

    assert(to_vector(output) == input);
    return output;
}


/*
 * Finds the byte offset within bitmap such that the byte_count bytes starting at offset are equal to
 * the byte_count at the start of bitmap shifted left by shift bits.
 */
int32_t find_offset_for_shift2(const vector<uint8_t> &bitmap, unsigned shift, unsigned byte_count) {
    vector<uint8_t> shifted(bitmap.begin(), bitmap.begin() + byte_count);
    vector<bool> bitvector = to_vector(shifted);
    bitvector.insert(bitvector.begin(), shift, false);
    assert(bitvector.size() == (byte_count * 8) + shift);
    bitvector.resize(byte_count * 8);
    shifted = to_vector(bitvector);

    auto res = std::search(bitmap.begin(), bitmap.end(), shifted.begin(), shifted.end());
    return res == bitmap.end() ? -1 : res - bitmap.begin();
}

/* use bitmap manipulation rather than math to sieve out the multiples */
uint64_t ProcessBitmap(uint64_t start, uint64_t end, const BYTE initial[NUM_PRIME], callback* fn) {

    // decrement by 2 initially because the other functions find the first possible prime not from "start"
    // but from start - 2
    start -= 2;

    uint64_t count = 0;

    // hardcode for 3 primes
    unsigned shift3 = initial[0];
    unsigned shift5 = initial[1];
    unsigned shift7 = initial[2];

    for (int64_t cur = start, i = 0; cur >= (int64_t)end; cur -= 128, i++) {
#if VERBOSE
        printf("PB  iter=%lu cur=%3lu end=%lu s3=%d s5=%d s7=%d\n",
                i, cur, end, shift3, shift5, shift7);
#endif

        uint64_t bitmap3 = BITMAPS[0] << shift3;
        uint64_t bitmap5 = BITMAPS[1] << shift5;
        uint64_t bitmap7 = BITMAPS[2] << shift7;

        uint64_t accum = bitmap3 | bitmap5 | bitmap7;

        shift3 = (shift3 + (3 - 64 % 3)) % 3;
        shift5 = (shift5 + (5 - 64 % 5)) % 5;
        shift7 = (shift7 + (7 - 64 % 7)) % 7;

        // the last iteration might have "incomplete" bitmap, we need to mask off values that are lower
        // than end
        uint64_t rem = (cur - end) / 2 + 1;
        if (rem < 64) {
            accum |= ~0UL << rem;
        }

        if (fn != NULL) {
            do_bitmap_callback(cur, accum, fn);
        }
        count += 64 - __builtin_popcountl(accum);
    }

    return count;
}

/* generic bitmap function */
uint64_t AlgoBitmap1(uint64_t start, uint64_t end, const BYTE initial[NUM_PRIME], callback* fn) {

    // decrement by 2 initially because the other functions find the first possible prime not from "start"
    // but from start - 2
    start -= 2;

    uint64_t count = 0;

    unsigned shift[NUM_PRIME];
    for (int i = 0; i < NUM_PRIME; i++)
        shift[i] = initial[i];

    for (int64_t cur = start, iter = 0; cur >= (int64_t)end; cur -= 128, iter++) {
#if VERBOSE
        printf("PB2 iter=%lu cur=%3lu end=%lu s3=%d s5=%d s7=%d\n",
                iter, cur, end, shift[0], shift[1], shift[2]);
#endif

        uint64_t accum = 0;

        for (int i = 0; i < NUM_PRIME; i++) {
            unsigned s = shift[i];
            uint64_t bitmap = s <= 63 ? BITMAPS[i] << s : 0;
            accum |= bitmap;
            uint8_t prime = primes[i] + 1;
            shift[i] = (s + (prime - 64 % prime)) % prime;
        }

        // the last iteration might have "incomplete" bitmap, we need to mask off values that are lower
        // than end
        uint64_t rem = (cur - end) / 2 + 1;
        if (rem < 64) {
            accum |= ~0UL << rem;
        }

        if (fn != NULL) {
            do_bitmap_callback(cur, accum, fn);
        }
        count += 64 - __builtin_popcountl(accum);
    }

    return count;
}

/* generic bitmap function */
uint64_t Bitmap2(uint64_t start, uint64_t end, const BYTE initial[NUM_PRIME], callback* fn) {

    // decrement by 2 initially because the other functions find the first possible prime not from "start"
    // but from start - 2
    start -= 2;

    uint64_t count = 0;

    uint8_t offset_indexes[NUM_PRIME];
    for (int i = 0; i < NUM_PRIME; i++) {
        offset_indexes[i] = SHIFT_TO_OFFSET_IDX[i][initial[i]];
    }

    for (int64_t cur = start, iter = 0; cur >= (int64_t)end; cur -= 128, iter++) {
#if VERBOSE
        printf("PB3 iter=%lu cur=%3lu end=%lu\n",
                iter, cur, end);
#endif

        uint64_t accum = 0;

        for (unsigned i = 0; i < NUM_PRIME; i++) {
            uint8_t oidx = offset_indexes[i];
            uint8_t offset = BYTE_OFFSETS[i][oidx];
            uint64_t bitmap = read64(BYTE_BITMAPS[i], offset);
            accum |= bitmap;
            oidx++;
            offset_indexes[i] = (oidx == OFFSET_PERIODS[i] ? 0 : oidx);
        }

#if VERBOSE
        printf("ac=%lX", accum);
#endif

        // the last iteration might have "incomplete" bitmap, we need to mask off values that are lower
        // than end
        uint64_t rem = (cur - end) / 2 + 1;
        if (rem < 64) {
            accum |= ~0UL << rem;
        }

#if VERBOSE
        printf(" ac-masked=%lX\n", accum);
#endif

        if (fn != NULL) {
            do_bitmap_callback(cur, accum, fn);
        }
        count += 64 - __builtin_popcountl(accum);
    }

    return count;
}

uint64_t handle_primes(uint64_t accum, int64_t cur, uint64_t end, callback *fn) {
    // the last iteration might have "incomplete" bitmap, we need to mask off values that are lower
    // than end
    int64_t rem = (cur - (int64_t)end) / 2 + 1;

    if (rem <= 0) {
        return 0;
    }

    if (rem < 64) {
        accum |= ~0UL << rem;
    }

    if (fn != NULL) {
        do_bitmap_callback(cur, accum, fn);
    }
    return 64 - __builtin_popcountl(accum);
}



/* does 4 consecutive 64-bit reads to emulate a 256-bit read */
uint64_t ProcessBitmap256(uint64_t start, uint64_t end, const BYTE initial[NUM_PRIME], callback* fn) {

    // decrement by 2 initially because the other functions find the first possible prime not from "start"
    // but from start - 2
    start -= 2;

    uint64_t count = 0;

    uint8_t offset_indexes[NUM_PRIME];
    for (int i = 0; i < NUM_PRIME; i++) {
        offset_indexes[i] = WrapTraits<32>::SHIFT_TO_OFFSET_IDX[i][initial[i]];
    }

    for (int64_t cur = start, iter = 0; cur >= (int64_t)end; cur -= 128 * 4, iter++) {

        uint64_t a0 = 0, a1 = 0, a2 = 0, a3 = 0;

        for (unsigned i = 0; i < (unsigned)NUM_PRIME; i++) {
            uint8_t oidx = offset_indexes[i];
            uint8_t offset = BYTE_OFFSETS256[i][oidx];
            //            assert((int)offset + 24 < 256);
            a0 |= read64(BYTE_BITMAPS256[i], offset +  0U);
            a1 |= read64(BYTE_BITMAPS256[i], offset +  8U);
            a2 |= read64(BYTE_BITMAPS256[i], offset + 16U);
            a3 |= read64(BYTE_BITMAPS256[i], offset + 24U);
            //            assert(BYTE_BITMAPS256[i][offset +  0] != 0xFF);
            //            assert(BYTE_BITMAPS256[i][offset + 31] != 0xFF);
            oidx++;
            offset_indexes[i] = (oidx == OFFSET_PERIODS256[i] ? 0 : oidx);
        }

#if VERBOSE
        printf("a0=%lX\na1=%lX\na2=%lX\na3=%lX\n", a0, a1, a2, a3);
#endif

        count += handle_primes(a0, cur,       end, fn);
        count += handle_primes(a1, cur - 128, end, fn);
        count += handle_primes(a2, cur - 256, end, fn);
        count += handle_primes(a3, cur - 384, end, fn);
    }

    return count;
}

uint64_t do_primes(uint64_t start, uint64_t end, process_fn *method, callback* c, bool message = true) {

    BYTE state[NUM_PRIME];
    make_state_table(start, state);

    int64_t count = method(start, end, state, c);

    if (message) {
        printf("Count for thread %lu, density %0.3f\n", (uint64_t)count,
                (double)count / (start - end));
    }

    return count;
}

/* handle BLOCK_COUNT worth of primes */
template <size_t SIZE>
uint64_t handle_primes_array(const uint8_t *bitmap, int64_t cur, callback *fn) {
    if (fn) {
        for (uint64_t offset = 0; offset < SIZE; offset += 8) {
            do_bitmap_callback(cur, read64(bitmap, offset), fn);
            cur -= 128;
        }
    }

    uint64_t popcnt = 0;
    for (uint64_t offset = 0; offset < SIZE; offset += 8) {
        popcnt += 64 - __builtin_popcountl(read64(bitmap, offset));
    }
    return popcnt;
}

__attribute__ ((noinline))
void kernel_c(uint8_t *offset_indexes, uint8_t *block_result) {
    for (uint64_t block = 0; block < BLOCK_COUNT; block++) {
        uint64_t a0 = 0, a1 = 0, a2 = 0, a3 = 0;
        for (unsigned i = 0; i < (unsigned) (NUM_PRIME); i++) {
            uint8_t oidx = offset_indexes[i];
            uint8_t offset = BYTE_OFFSETS256[i][oidx];
            a0 |= read64(BYTE_BITMAPS256[i], offset + 0U);
            a1 |= read64(BYTE_BITMAPS256[i], offset + 8U);
            a2 |= read64(BYTE_BITMAPS256[i], offset + 16U);
            a3 |= read64(BYTE_BITMAPS256[i], offset + 24U);
            //            assert(BYTE_BITMAPS256[i][offset +  0] != 0xFF);
            //            assert(BYTE_BITMAPS256[i][offset + 31] != 0xFF);
            oidx++;
            offset_indexes[i] = (oidx == OFFSET_PERIODS256[i] ? 0 : oidx);
        }
        size_t block_offset = block * 32;
        write64(a0, block_result, block_offset + 0U);
        write64(a1, block_result, block_offset + 8U);
        write64(a2, block_result, block_offset + 16U);
        write64(a3, block_result, block_offset + 24U);
    }
}

/* a version of ProcessBitmap256 which only works in chunks of 512 candidates and delegates the tail
 * to a simpler, slower function */
template<kfunc K, size_t READ_BYTES>
uint64_t WrapBitmap(uint64_t start_, uint64_t end_, const BYTE initial[NUM_PRIME], callback* fn) {

    constexpr size_t VALUES_PER_BLOCK = READ_BYTES * 8 * 2; // 1 bit per value, and *2 since we skip even numbers

    int64_t start = start_, end = end_;

    // decrement by 2 initially because the other functions find the first possible prime not from "start"
    // but from start - 2
    start -= 2;

    uint64_t count = 0;

    uint8_t offset_indexes[NUM_PRIME];
    for (int i = 0; i < NUM_PRIME; i++) {
        //        offset_indexes[i] = SHIFT_TO_OFFSET_IDX256[i][initial[i]];
        offset_indexes[i] = WrapTraits<READ_BYTES>::SHIFT_TO_OFFSET_IDX[i][initial[i]];
    }

    int64_t cur = start;
    for (uint64_t iter = 0; cur - BLOCK_COUNT * (ssize_t)VALUES_PER_BLOCK >= end; cur -= BLOCK_COUNT * VALUES_PER_BLOCK, iter++) {
#if VERBOSE
        printf("PB256 iter=%lu cur=%3lu end=%lu\n",
                iter, cur, end);
#endif

        uint8_t block_result[BLOCK_COUNT * READ_BYTES];  // 8 bits * 2 skipping evens

        K(offset_indexes, block_result);

#if VERBOSE
        printf("a0=%lX\na1=%lX\na2=%lX\na3=%lX\n", a0, a1, a2, a3);
#endif

        count += handle_primes_array<BLOCK_COUNT * READ_BYTES>(block_result, cur, fn);
    }

    // we he handle the tail be delegating to ProcessBitmap2
    if (cur >= end) {
        cur += 2;
        count += do_primes(cur, end, AlgoBitmap1, fn, false);
    }

    return count;
}

void make_bitmaps() {
    printf("static const uint64_t BITMAPS[] = {\n");
    for (int i = 0; i < NUM_PRIME; i++) {
        unsigned prime = primes[i] + 1;
        uint64_t bitmap = 0;
        for (unsigned bit = 0; bit < 64; bit += prime) {
            bitmap |= (1UL << bit);
        }
        //        printf("\t/* %3d */ { 0x%016lx },\n", prime, bitmap);
        printf("\t/* %3d */  0x%016lx ,\n", prime, bitmap);
    }
    printf("};\n");
}

/* how many contiguous bytes a read will do */
/* this is the unit of granularity for offsets as well, i.e., each offset corresponds to a read of this size */
#define PB3_READ_BYTES 64

/* how many PB3_READ_BYTES quantities will be available in each bitmap */
#define PB3_BITMAP_UNROLL 1

/* how many offsets you can access without checking for wrap */
#define PB3_OFFSET_UNROLL  1

vector< vector<uint8_t> > byte_bitmaps() {
    vector< vector<uint8_t> > bitmaps;
    for (int i = 0; i < NUM_PRIME; i++) {
        bitmaps.emplace_back();
        auto& bitmap = bitmaps.back();
        unsigned prime = primes[i] + 1;
        unsigned byte_count = PB3_BITMAP_UNROLL * PB3_READ_BYTES + prime - 1;
        bitmap.resize(byte_count);
        for (unsigned bit = 0; bit < byte_count * 8; bit += prime) {
            unsigned byte_index = bit / 8;
            unsigned bit_index = bit % 8;
            assert((bitmap.at(byte_index) & (1 << bit_index)) == 0);
            bitmap.at(byte_index) |= (1 << bit_index);
        }
    }
    return bitmaps;
}

/* a shift that saturates to zero (rather than UB) for shift amounts >= 64 */
/* read a 64-bit value from a char array, (native endianness) */
/* write a 64-bit value to a char array at i */
void make_byte_bitmaps() {
    std::string suffix_str = std::to_string(PB3_READ_BYTES * 8);
    const char *suffix = suffix_str.c_str();

    vector< vector<uint8_t> > bitmaps = byte_bitmaps();
    assert(bitmaps.size() == NUM_PRIME);
    unsigned longest_bitmap = bitmaps.back().size();
    printf("extern \"C\" const uint8_t BYTE_BITMAPS%s[][%d] = {\n", suffix, longest_bitmap);
    for (int i = 0; i < NUM_PRIME; i++) {
        unsigned prime = primes[i] + 1;
        vector<uint8_t> bitmap = bitmaps.at(i);
        printf("\t/* i=%2d p=%3d */  {", i, prime);
        for (uint8_t byte : bitmap) {
            printf("0x%02x, ", byte);
        }
        for (unsigned j = 0; j < longest_bitmap - bitmap.size(); j++) {
            printf("0x%02x, ", 0xFF);
        }
        printf("},\n");

    }
    printf("};\n");

    vector<int32_t> shift_to_offset_idx[NUM_PRIME] = {};
    uint8_t periods[NUM_PRIME] = {};

    unsigned largest_prime = primes[NUM_PRIME-1] + 1;
    unsigned longest_period =  largest_prime + PB3_OFFSET_UNROLL;
    // now make the rotating offset amounts
    printf("\n\nextern \"C\" const uint8_t BYTE_OFFSETS%s[][%d] = {\n", suffix, longest_period);
    for (int i = 0; i < NUM_PRIME; i++) {
        unsigned prime = primes[i] + 1;
        unsigned period = prime;
        while (period < PB3_OFFSET_UNROLL) {
            period += prime;
        }
        assert(period <= 255);
        periods[i] = period;
        unsigned offset_count = period + PB3_OFFSET_UNROLL;

        shift_to_offset_idx[i].resize(prime, -1);

        vector<uint8_t> bitmap = bitmaps.at(i);

        printf("\t/* %3d */  {", prime);

        for (unsigned j = 0, shift = 0; j < offset_count; j++) {
            assert(shift < prime);
            ////////////////////
            int32_t offset = find_offset_for_shift2(bitmap, shift, PB3_READ_BYTES);

            assert(offset + (PB3_READ_BYTES - 1U) < bitmap.size());
            assert(offset >= 0 && offset <= 255);
            printf("%3d,", offset);

            // we record the offset index for this shift (if it hasn't already been recorded)
            int32_t odix = shift_to_offset_idx[i].at(shift);
            if (odix == -1) {
                assert(j < prime);
                shift_to_offset_idx[i].at(shift) = j;
            }

            shift = (shift + (prime - (PB3_READ_BYTES * 8) % prime)) % prime;
        }

        printf("},\n");

        // check that all shift -> offset entires were filled
        for (auto oidx : shift_to_offset_idx[i]) {
            assert(oidx >= 0);
        }
    }

    printf("};\n");

    /* finally, record the periods */
    printf("\n\nextern \"C\" const uint8_t OFFSET_PERIODS%s[] = {\n", suffix);
    for (int i = 0; i < NUM_PRIME; i++) {
        printf("\t/* %3d */ %d,\n", primes[i] + 1, periods[i]);
    }
    printf("};\n");


    /* and really finally, record the shift -> offset table for converting the initial state */
    printf("\n\nextern \"C\" const uint8_t SHIFT_TO_OFFSET_IDX%s[][%d] = {\n", suffix, largest_prime);
    for (int i = 0; i < NUM_PRIME; i++) {
        printf("\t/* %3d */ { ", primes[i] + 1);
        for (auto oidx : shift_to_offset_idx[i]) {
            printf("%d, ", oidx);
        }
        printf("},\n");
    }
    printf("};\n");


}

unsigned int RunCpu(int threadnum)
{
    ULONGLONG b;

    switch (threadnum)
    {
    case 0:
    {
        BYTE state[32] = {
                0x00, 0x03, 0x00, 0x00, 0x07, 0x08, 0x0c, 0x0f,
                0x0d, 0x0d, 0x1c, 0x00, 0x10, 0x06, 0x20, 0x24,
                0x29, 0x18, 0x3d, 0x1c, 0x10, 0x04, 0x37, 0x53,
                0x3d, 0x40, 0x45, 0x39, 0x55, 0x1e, 0x62, 0x40 };

        b = ProcessA(69780348563, 131, state, CALLBACK);
        break;
    }
    case 1:
    {
        BYTE state[32] = {
                0x00, 0x02, 0x05, 0x0a, 0x0b, 0x0e, 0x0a, 0x03,
                0x08, 0x0b, 0x0d, 0x00, 0x23, 0x15, 0x1d, 0x28,
                0x33, 0x14, 0x19, 0x1f, 0x0a, 0x47, 0x04, 0x38,
                0x0a, 0x0f, 0x17, 0x6b, 0x34, 0x3e, 0x45, 0x87 };

        b = ProcessA(139560697001, 69780348563, state, CALLBACK);
        break;
    }
    case 2:
    {
        BYTE state[32] = {
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40,
                0x20, 0x10, 0x08, 0x04, 0x02, 0x01, 0x00, 0x00 };

        // Note: this isn't the correct start/end for this state, so
        // the 'likely' primes produced in ProcessA won't actually
        // be primes.  However this is still useful for testing since
        // it produces the largest loop count I know how to create
        // using 30 primes.
        b = ProcessA(501, 305, state, CALLBACK);
        break;
    }
    default:
    {
        b = 0;
        printf("Invalid thread id: %d\n", threadnum);
        break;
    }
    }

    printf("Count for thread %d: %lu\n", threadnum, (uint64_t)b);

    return 0;
}



/* given two prime finding algorithms, check that they produce the same results */
bool check_primes(uint64_t start, uint64_t end, process_fn *method1, process_fn *method2) {
    recording_callback cb1, cb2;
    record_init(&cb1);
    record_init(&cb2);

    do_primes(start, end, method1, (callback *)&cb1, false);
    do_primes(start, end, method2, (callback *)&cb2, false);

    // now validate that each algorithm got the same set of primes
    ssize_t sz1 = cb1.size, sz2 = cb2.size;
    ssize_t min = MIN(sz1, sz2);
    bool same = cb1.size == cb2.size;
    for (ssize_t i = 0; i < min; i++) {
        uint64_t r1 = cb1.array[i];
        uint64_t r2 = cb2.array[i];
        if (r1 != r2) {
            printf("Value mismatch at prime %lu: r1=%lu r2=%lu\n", i, r1, r2);
            same = false;
            // print some values around the mismatch
            for (ssize_t j = std::max(i - 5, 0L); j < std::min(i + 5, min); j++) {
                printf("%s", i ==j ? ">>>" : "   ");
                printf("i=%lu p1=%lu p2=%lu\n", j, cb1.array[j], cb2.array[j]);
            }
            break;
        }
    }

    if (same) {
#if VERBOSE
        printf("check_primes OK!\n"); // sz1=%lu sz2=%lu sz1=%lu sz2=%lu", sz1, sz2, cb1.size, cb2.size);
#endif
    } else {
        printf("mismatched lists start=%lu end=%lu sz1=%lu sz2=%lu", start, end, sz1, sz2);
    }

    return same;
}

void do_time(uint64_t start, uint64_t end, process_fn *method) {
    if ((start & 1) == 0)
        start++;
    clock_t before, after;
    BYTE state[NUM_PRIME];
    make_state_table(start, state);
    before = clock();
    uint64_t count = method(start, end, state, NULL);
    after = clock();
    double delta = ((double)after - before) / CLOCKS_PER_SEC;
    printf("Count %ld, density=%5.3f elapsed %5.3fs, cycles/candidate %4.2f, cycles/prime %4.2f\n",
            count,
            (double)count / (start - end),
            delta,
            delta * MHZ * 1000000 / (start - end),
            delta * MHZ * 1000000 / count
    );
}

int main(int argc, char **argv)
{
    assert(CHAR_BIT == 8);

    if (argc < 2) {
        fprintf(stderr, "Usage: %s [check | time | tables] arguments...\n"
                "\t time ALGORITHM: time the given algorithm\n"
                "\t tables: generate some tables needed by various algorithms\n",
                argv[0]);
        exit(1);
    }

    std::string command = argv[1];

    process_fn *method = nullptr;
    if (argc >= 3) {
        std::string algo = argv[2];
        if (algo == "ProcessA") {
            method = ProcessA;
        } else if (algo == "ProcessA2") {
            method = ProcessA2;
        } else if (algo == "ProcessA2") {
            method = ProcessA2;
        } else if (algo == "ProcessC") {
            method = ProcessC;
        } else if (algo == "Bitmap1") {
            method = AlgoBitmap1;
        } else if (algo == "Bitmap2") {
            method = Bitmap2;
        } else if (algo == "asm256") {
            method = WrapBitmap<kernel256_asm, 32>;
        } else if (algo == "asm512") {
            method = WrapBitmap<kernel512_asm, 64>;
        } else if (algo == "Bitmap256-kernel") {
            method = WrapBitmap<kernel_c, 32>;
        } else if (algo == "Bitmap256") {
            method = ProcessBitmap256;
        } else {
            fprintf(stderr, "Unknown algorithm: %s\n", algo.c_str());
        }

    }

    if (command == "tables") {
        make_byte_bitmaps();
        exit(0);
    }

    if (command == "time") {
        if (!method) {
            fprintf(stderr, "No algorithm provided!\n");
        } else {
            size_t start = (argc == 4 ? std::stol(argv[3]) : 15111111111L);
            size_t end = 131;
            printf("Finding all likely primes between %lu and %lu...\n", start, end);
            do_time(start, end, method);
        }
        exit(0);
    }

    //    test_state();
    //    return 0;
    //    make_bitmaps();
    //    make_byte_bitmaps();
    //    return 0;

    //    69780348563, end = 69780340563

    //    check_primes(101010101, 131, ProcessC, AlgoBitmap1);

    //    print_state(15);
    int64_t start =     331;
    int64_t   end =     131;
    //
    callback cb = print_callback;
    //
    ////    do_primes(start, end, ProcessC, (callback *)&cb);
    //    do_primes(start, end, Bitmap2, (callback *)&cb);


    if (argc == 1) {
        for (start = 133; start < 100000; start += 100) {
            if (!check_primes(start, std::max(131L, start - 90000), ProcessBitmap256, WrapBitmap<kernel512_asm, 64>)) {
                exit(1);
            }
        }
        printf("Verified OK!\n");
    }

    //    return 0;

    uint64_t iters = 1511111111;

    if (argc == 2 || true) {
        do_time(iters, 131, WrapBitmap<kernel512_asm, 64>);
    } else {
        do_time(iters, 131, ProcessBitmap256);
        //        do_time(iters, 131, ProcessA);
        //        do_time(iters, 131, ProcessA2);
    }

    //        do_time(iters, 131, ProcessA2);
    //        do_time(iters, 131, Bitmap2);
    return 0;
}
