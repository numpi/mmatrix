GTH algorithm implemented in Fortran + Matlab mex-file interface.

To build the Mex files (with a working Fortran setup and a modern Matlab):
```matlab
>> build_mex
```

To run automated tests and check everything is working:
```matlab
>> runtests('tests/');
```

To run a more detailed interactive demo:
```matlab
>> main_matlab
```

Descriptions of the variants:

* `dmmtf2.f`: basic unblocked (BLAS level-2) version.
* `dmmtrf2.f`: recursive cache-oblivious version. Needs a temporary buffer `WORK` of size `2n`.
* `dmmtrf.f`: blocked version (BLAS level-3) which calls `dmmtf2.f` on the single blocks. Needs a temporary buffer `WORK` of size `nb`, where `nb` is the user-specified block size.
* `dmmtrff.f`: blocked version (BLAS level-3) which calls `dmmtrf2.f` on the single blocks (this one should be the fastest). Needs a temporary buffer `WORK` of size `2*nb`, where `nb` is the user-specified block size.
