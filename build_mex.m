function build_mex()
%Compile Fortran/C files 

    % List of required MEX files
    mexFiles = {'dmmtf2_mex.mexw64', 'dmmtrf_mex.mexw64', 'dmmtrff_mex.mexw64'};
    
    % Check which files are missing
    missing = cellfun(@(f) ~isfile(f), mexFiles);
    
    if ~any(missing)
        fprintf('All MEX files are already available. Skipping compilation.\n');
        return
    end
    
    fprintf('Some MEX files are missing. Starting compilation...\n');
    
    % Delete old object/MEX files (optional)
    % delete('*.obj');
    % delete('*.mexw64');

    fprintf('=== Compiling Fortran to objects ===\n');
    fortran_files = {'dmmtf2.f', 'dmmtrf.f', 'dmmtrff.f'};

    % Add dgetrf files if they exist
    if exist('dgetrf.f', 'file'),   fortran_files{end+1} = 'dgetrf.f'; end
    if exist('dgetrf2.f', 'file'),  fortran_files{end+1} = 'dgetrf2.f'; end

    % Compile .f files to .obj
    for i = 1:length(fortran_files)
        fprintf('Compiling %s...\n', fortran_files{i});
        eval(['mex -c -largeArrayDims ' fortran_files{i}]);
    end

    fprintf('\n=== Linking MEX files ===\n');

    % dmmtf2_mex
    if isfile('dmmtf2_mex.c')
        fprintf('Linking dmmtf2_mex...\n');
        mex -largeArrayDims dmmtf2_mex.c dmmtf2.obj -lmwlapack -lmwblas
    end

    % dmmtrf_mex (needs dmmtf2.obj)
    if isfile('dmmtrf_mex.c')
        fprintf('Linking dmmtrf_mex...\n');
        mex -largeArrayDims dmmtrf_mex.c dmmtrf.obj dmmtf2.obj -lmwlapack -lmwblas
    end

    % dmmtrff_mex
    if isfile('dmmtrff_mex.c')
        fprintf('Linking dmmtrff_mex...\n');
        mex -largeArrayDims dmmtrff_mex.c dmmtrff.obj -lmwlapack -lmwblas
    end

    % dgetf2_mex
    if exist('dgetf2_mex.c', 'file')
        fprintf('Linking dgetf2_mex...\n');
        mex -largeArrayDims dgetf2_mex.c -lmwlapack -lmwblas
    end

    % dgetrf2_mex
    if exist('dgetrf2_mex.c', 'file')
        fprintf('Linking dgetrf2_mex...\n');
        mex -largeArrayDims dgetrf2_mex.c -lmwlapack -lmwblas
    end

    % dgetrf_mex
    if exist('dgetrf_mex.c', 'file') && exist('dgetrf.obj', 'file')
        fprintf('Linking dgetrf_mex...\n');
        mex -largeArrayDims dgetrf_mex.c dgetrf.obj dgetrf2.obj -lmwlapack -lmwblas
    end

    fprintf('\n=== All compilations complete! ===\n');
end
