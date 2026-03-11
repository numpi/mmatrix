function build_mex()
%Compile Fortran/C files 

    fprintf('=== Compiling Fortran to objects ===\n');
    fortran_files = { ...
        "dmmtf2.f", "dmmtrf.f", "dmmtrff.f", "dmmtrf2.f", ...
        "dmltf2.f", "dmltrf.f", "dmltrff.f", "dmltrf2.f", ...
        "dgetf2.f", "dgetrf.f", "dgetrf2.f"};

    for i = 1:length(fortran_files)
        fprintf('Compiling %s...\n', fortran_files{i});
        mex("-c", "-largeArrayDims", fortran_files{i});
    end

    if isfile('dmmtf2.o')
        objext = 'o';
    elseif isfile('dmmtf2.obj')
        objext = 'obj';
    else
        error('unknown compiled object file extension');
    end

    fprintf('\n=== Linking MEX files ===\n');

    link_gateway("dmmtf2", {}, {"HAS_TRIPLET"});
    link_gateway("dmmtrf2", {}, {"HAS_TRIPLET", "HAS_WORK"});
    link_gateway("dmmtrf", {"dmmtf2"}, {"HAS_TRIPLET", "HAS_WORK", "HAS_NB"});
    link_gateway("dmmtrff", {"dmmtrf2"}, {"HAS_TRIPLET", "HAS_WORK", "HAS_NB", "WORK_SIZE_2NB"});

    link_gateway("dmltf2", {}, {"HAS_TRIPLET", "LEFT_TRIPLET"});
    link_gateway("dmltrf2", {}, {"HAS_TRIPLET"});
    link_gateway("dmltrf", {"dmltf2"}, {"HAS_TRIPLET", "LEFT_TRIPLET", "HAS_WORK", "HAS_NB"});
    link_gateway("dmltrff", {"dmltrf2"}, {"HAS_TRIPLET", "LEFT_TRIPLET", "HAS_WORK",  "HAS_NB"});

    link_gateway("dgetf2", {}, {"HAS_IPIV"});
    link_gateway("dgetrf2", {}, {"HAS_IPIV"});
    link_gateway("dgetrf", {"dgetrf2"}, {"HAS_IPIV"});
    
    fprintf('\n=== All compilations complete! ===\n');

    function link_gateway(funname, dependencies, switches)
        fprintf('Linking %s_mex...\n', funname);
        dep_objects = cellfun(@(x) x+"."+objext, dependencies, UniformOutput=false);
        additional_parameters = cellfun(@(x) "-D"+x, switches, UniformOutput=false);
        mex("-largeArrayDims", "mex_gateway.c", funname+"."+objext, dep_objects{:}, "-lmwlapack", "-lmwblas", "-output", funname+"_mex", "-DFUNNAME="+funname, additional_parameters{:});
        % build also the _solve version
        solvername = funname + "_mexsolve";
        fprintf('Linking %s...\n', solvername);
        mex("-largeArrayDims", "mex_gateway.c", funname+"."+objext, dep_objects{:}, "-lmwlapack", "-lmwblas", "-output", solvername, "-DFUNNAME="+funname, additional_parameters{:}, "-DSOLVE_RHS");
    end
end
