# MMatrix Matlab toolbox

A toolbox for fast and accurate computations using M-matrices in Matlab, by Bruno Iannazzo, Mehdi Najafi-Kalyani, and Federico Poloni.

## Installation

To install the toolbox, download the source files from the [Releases page](https://github.com/numpi/mmatrix/releases/) into a folder, and add that folder to the Matlab search path.

If you have a working Fortran compiler set up on your computer to use with Matlab, you can compile the Mex files yourself with
```matlab
>> build_mex
```
Otherwise, you can download pre-compiled versions for Windows and Linux from the [Releases page](https://github.com/numpi/mmatrix/releases/). (Producing a pre-compiled version for Mac is more complicated, unfortunately, since Matlab does not support building Mex-files using GFortran on Mac machines.) Download the `mex-???.zip` file that matches your operating system, and unzip its contents into the `+mmex` folder.

After installation, you can run the automated tests to check everything is working:
```matlab
>> runtests('tests/');
```

### Automated installation inside Matlab

To download and install the most recent development version of the toolbox automatically, you can also run:

```matlab
unzip('https://github.com/numpi/mmatrix/archive/main.zip');
movefile('mmatrix-main', 'mmatrix-toolbox');
addpath(fullfile('mmatrix-toolbox'));
savepath;
```

### Descriptions of the Fortran files

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
