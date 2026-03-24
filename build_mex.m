function build_mex()
%Compile Fortran/C files 

    fprintf('=== Compiling Fortran to objects ===\n');
    fortran_files = { ...
        "dmmtf2.f", "dmmtrf.f", "dmmtrff.f", "dmmtrf2.f", ...
        "dmltf2.f", "dmltrf.f", "dmltrff.f", "dmltrf2.f", ...
        "dgetf2_reference.f", "dgetrf_reference.f", "dgetrf2_reference.f"};

    for i = 1:length(fortran_files)
        fprintf('Compiling %s...\n', fortran_files{i});
        mex("-c", "-largeArrayDims", "-outdir", "+mmex", "fortran/"+fortran_files{i});
    end

    if isfile('+mmex/dmmtf2.o')
        objext = 'o';
    elseif isfile('+mmex/dmmtf2.obj')
        objext = 'obj';
    else
        error('unknown compiled object file extension');
    end

    fprintf('\n=== Linking MEX files ===\n');

    link_gateway("dmmtf2", {"dmmtf2"}, {"HAS_TRIPLET"});
    link_gateway("dmmtrf2", {"dmmtrf2"}, {"HAS_TRIPLET", "HAS_WORK"});
    link_gateway("dmmtrf", {"dmmtf2", "dmmtrf"}, {"HAS_TRIPLET", "HAS_WORK", "HAS_NB"});
    link_gateway("dmmtrff", {"dmmtrf2", "dmmtrff"}, {"HAS_TRIPLET", "HAS_WORK", "HAS_NB", "WORK_SIZE_2NB"});

    link_gateway("dmltf2", {"dmltf2"}, {"HAS_TRIPLET", "LEFT_TRIPLET"});
    link_gateway("dmltrf2", {"dmltrf2"}, {"HAS_TRIPLET"});
    link_gateway("dmltrf", {"dmltf2", "dmltrf"}, {"HAS_TRIPLET", "LEFT_TRIPLET", "HAS_WORK", "HAS_NB"});
    link_gateway("dmltrff", {"dmltrf2", "dmltrff"}, {"HAS_TRIPLET", "LEFT_TRIPLET", "HAS_WORK",  "HAS_NB"});

    % Native Lapack versions
    link_gateway("dgetf2", {}, {"HAS_IPIV"});
    link_gateway("dgetrf2", {}, {"HAS_IPIV"});
    link_gateway("dgetrf", {}, {"HAS_IPIV"});

    % Built from sources
    link_gateway("dgetf2_reference", {"dgetf2_reference"}, {"HAS_IPIV"});
    link_gateway("dgetrf2_reference", {"dgetrf2_reference"}, {"HAS_IPIV"});
    link_gateway("dgetrf_reference", {"dgetrf2_reference", "dgetrf_reference"}, {"HAS_IPIV"});

    fprintf('\n=== All compilations complete! ===\n');

    function link_gateway(funname, dependencies, switches)
        fprintf('Linking %s...\n', funname);
        dep_objects = cellfun(@(x) "+mmex/"+x+"."+objext, dependencies, UniformOutput=false);
        additional_parameters = cellfun(@(x) "-D"+x, switches, UniformOutput=false);
        mex("-largeArrayDims", "-outdir", "+mmex", "mex_gateway.c", dep_objects{:}, "-lmwlapack", "-lmwblas", "-output", funname, "-DFUNNAME="+funname, additional_parameters{:});
        % build also the _solve version
        solvername = funname + "_solve";
        fprintf('Linking %s...\n', solvername);
        mex("-largeArrayDims", "-outdir", "+mmex", "mex_gateway.c", dep_objects{:}, "-lmwlapack", "-lmwblas", "-output", solvername, "-DFUNNAME="+funname, additional_parameters{:}, "-DSOLVE_RHS");
    end
end
