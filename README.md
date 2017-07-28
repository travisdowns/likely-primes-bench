# likely-primes-bench
Code from https://codereview.stackexchange.com/questions/169819/sse-loop-to-walk-likely-primes

# Make

Run `make`.

# Running 
Run it like `likely-primes` to see some basic usage info.

## Time Command
To time an algorithm you can use the `time` command:

    ./likely-primes time Bitmap2 69780348563
    
`Bitmap2` is the algorithm, and `69780348563` is the prime to start at (it will find all the primes between this value and `131`).

## Check Command

You can check that the specified algorithm returns exactly the same primes as a reference algorithm with `check`:

```
$./likely-primes check asm512
Verified OK!
```

It is not exhaustive, but it checks a large number of primes with various start and end values.

## Tables Command

You can generate the set of tables used by a specific algorithm using the `tables` command. You adjust the `RB_*` constants to pick what tables are generated and to control various factors such as the maximum supported unroll factor. 

# Algorithms

The following algorithms are supported:

 - ProcessA
 
   The original algorithm posted by the OP, with small edits to convert it to NASM syntax and Linux calling convention.
 - ProcessA2
 
   A branch-free version of the orignial algorithm which just does the counting in the inner loop and avoids the unpredictable exits to the outer loop
  - ProcessC
  
   A C version of the original algorithm.
  - Bitmap1
   The original 64-bit bitmap + shifts version, slow with two divisions and branch mispredicts.
  - Bitmap2
   A table-driven version of Bitmap1.
  - Bitmap256
   Kind of a C version of the later 32-byte reading SIMD versions that use lookup tables only (no shifts). 
  - Bitmap256-kernel
   The same as `Bitmap256` except using the generic wrapper around the core code, as used also by the SIMD versions.
  - asm256
   A fast lookup-based AVX2 bitmap algorithm that reads 32B in its inner loop.
  - asm512
   Same as `asm256` but reads 64 bytes and so needs different lookup tables.
    
