## GTH algorithm implemented in Fortran + Matlab mex-file interface.

If you have a working Fortran compiler setup, you can compile the Mex files yourself with
```matlab
>> build_mex
```
Otherwise, you can download pre-compiled Linux and Windows versions from [here](https://github.com/numpi/mmatrix/releases/). Unzip the files in the `+mmex` folder to use them.

To run automated tests and check everything is working:
```matlab
>> runtests('tests/');
```

To run a more detailed interactive demo:
```matlab
>> main_matlab
```

To download and install the toolbox directly from GitHub, run:

```matlab
unzip('https://github.com/numpi/mmatrix/archive/main.zip');
movefile('mmatrix-main', 'mmatrix-toolbox');
addpath(fullfile('mmatrix-toolbox'));
savepath;
```

### Descriptions of the variants:

#### Right triplet algorithms:

* `dmmtf2.f`: basic unblocked (BLAS level-2) version.
* `dmmtrf2.f`: recursive cache-oblivious version. Needs a temporary buffer `WORK` of size `2n`.
* `dmmtrf.f`: blocked version (BLAS level-3) which calls `dmmtf2.f` on the single blocks. Needs a temporary buffer `WORK` of size `nb`, where `nb` is the user-specified block size.
* `dmmtrff.f`: blocked version (BLAS level-3) which calls `dmmtrf2.f` on the single blocks (this one should be the fastest). Needs a temporary buffer `WORK` of size `2*nb`, where `nb` is the user-specified block size.


#### Left triplet algorithms:

* `dmltf2.f`: basic unblocked (BLAS level-2) version.
* `dmltrf2.f`: recursive cache-oblivious version.
* `dmltrf.f`: blocked version (BLAS level-3) which calls `dmltf2.f` on the single blocks. Needs a temporary buffer `WORK` of size `nb`, where `nb` is the user-specified block size.
* `dmltrff.f`: blocked version (BLAS level-3) which calls `dmltrf2.f` on the single blocks. Needs a temporary buffer `WORK` of size `nb`, where `nb` is the user-specified block size.
