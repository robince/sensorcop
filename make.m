%
% make.m - build script for MatlabAPI_lite
%

% set to false if you don't have PyF95++ and are just using checked out
% sources
PROCESS_TEMPLATES = false;
% mex options
DEBUG = false;
VERBOSE = true;
LARGEARRAY = true;

% path to matlab api
% MATLABAPI_DIR = '/Users/robince/git/MatlabAPI_lite';
% MATLABAPI_DIR = '/home/robini/code/MatlabAPI_lite_grid';
% MATLABAPI_DIR = '/home/robini/code/MatlabAPI_lite';
% MATLABAPI_DIR = '/home/robini/code/MatlabAPI_lite_gfortran44';
% MATLABAPI_DIR = '/home/robini/code/MatlabAPI_lite_gfortran';
MATLABAPI_DIR = '/home/robini/code/MatlabAPI_lite_ifort15';
% MATLABAPI_DIR = '/Users/robince/git/MatlabAPI_lite_gfortran';

% needed for some compilers (eg intel)
EXTRA_LIBS = '';
% EXTRA_LIBS = '/opt/intel/composerxe/lib/intel64/libifcoremt_pic.a';

% Matlab doesn't use the normal system path so need the full path to
% executables
PYF95 = '/Users/robince/slash/bin/PyF95++';
% PyF95++ installs to a different python directory so need to set
% PYTHONPATH appropriately
PYF95_PP = '/Users/robince/slash/lib/python2.7/site-packages';
setenv('PYTHONPATH', PYF95_PP);
% I had to change the !# line in the PyF95++ to explicitly point to the
% current python since env inside Matlab gives the wrong one because of the
% path

if ispc
    OBJEXT = 'obj';
else
    OBJEXT = 'o';
end

%% Clean

if PROCESS_TEMPLATES
    delete *.F90
    delete *.f90
end
delete *.mod
delete *.i90
delete(['*.' OBJEXT])
delete(['*.' mexext])
% cd tests
% if PROCESS_TEMPLATES
%     delete *.F90
%     delete *.f90
% end
% delete *.mod
% delete(['*.' OBJEXT])
% delete(['*.' mexext])
% cd ..
%% Process templates

if PROCESS_TEMPLATES
    st = system([ PYF95 ' --sources="fcinfo.F90T instantiate.F90T copnorm.F90T copnorm_slice_omp.F90T" --std="f03" -x'],'-echo');
    if st ~= 0
        error('Problem processing MatlabAPImx.F90T template')
    end
end
%% Build mex

ARGS = {};
if LARGEARRAY,  ARGS{end+1} = '-largeArrayDims';
else            ARGS{end+1} = '-compatibleArrayDims'; end
if DEBUG,       ARGS{end+1} = '-g'; end
if VERBOSE,     ARGS{end+1} = '-v'; end
ARGS{end+1} = ['-I' MATLABAPI_DIR];
% ARGS{end+1} = ['-I/opt/intel/composer_xe_2011_sp1.9.289/mkl/include/'];

typekind = {'c_double' 'c_float'};
%%
% lib_array
MEXARGS = ARGS;
MEXARGS{end+1} = '-c';
MEXARGS{end+1} = 'lib_array.f';
mex(MEXARGS{:})

%%
MEXARGS = ARGS;
MEXARGS{end+1} = 'rocarea.f';
MEXARGS{end+1} = fullfile(MATLABAPI_DIR,['MatlabAPImx.' OBJEXT]);
MEXARGS{end+1} = fullfile(MATLABAPI_DIR,['MatlabAPImex.' OBJEXT]);
MEXARGS{end+1} = ['lib_array.' OBJEXT];
MEXARGS{end+1} = '-lmwblas';
MEXARGS{end+1} = '-lmwlapack';
if ~isempty(EXTRA_LIBS)
    MEXARGS{end+1} = EXTRA_LIBS;
end
mex(MEXARGS{:})

%%
MEXARGS = ARGS;
MEXARGS{end+1} = 'kstest_slice_omp.f';
MEXARGS{end+1} = fullfile(MATLABAPI_DIR,['MatlabAPImx.' OBJEXT]);
MEXARGS{end+1} = fullfile(MATLABAPI_DIR,['MatlabAPImex.' OBJEXT]);
MEXARGS{end+1} = ['lib_array.' OBJEXT];
MEXARGS{end+1} = '-lmwblas';
MEXARGS{end+1} = '-lmwlapack';
if ~isempty(EXTRA_LIBS)
    MEXARGS{end+1} = EXTRA_LIBS;
end
mex(MEXARGS{:})

%%
MEXARGS = ARGS;
MEXARGS{end+1} = 'kstest2d_slice_omp.f';
MEXARGS{end+1} = fullfile(MATLABAPI_DIR,['MatlabAPImx.' OBJEXT]);
MEXARGS{end+1} = fullfile(MATLABAPI_DIR,['MatlabAPImex.' OBJEXT]);
MEXARGS{end+1} = ['lib_array.' OBJEXT];
MEXARGS{end+1} = '-lmwblas';
MEXARGS{end+1} = '-lmwlapack';
if ~isempty(EXTRA_LIBS)
    MEXARGS{end+1} = EXTRA_LIBS;
end
mex(MEXARGS{:})

%%
MEXARGS = ARGS;
MEXARGS{end+1} = 'kstest2d_slice_omp.f';
MEXARGS{end+1} = fullfile(MATLABAPI_DIR,['MatlabAPImx.' OBJEXT]);
MEXARGS{end+1} = fullfile(MATLABAPI_DIR,['MatlabAPImex.' OBJEXT]);
MEXARGS{end+1} = ['lib_array.' OBJEXT];
MEXARGS{end+1} = '-lmwblas';
MEXARGS{end+1} = '-lmwlapack';
mex(MEXARGS{:})
