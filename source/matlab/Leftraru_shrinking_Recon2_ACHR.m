addpath(genpath('/home/aacevedo'));
initCobraToolbox;
changeCobraSolver('gurobi');
myPath = '/home/aacevedo/NLHPC/models';

load /home/aacevedo/NLHPC/models/small_metabolism.mat

%recon = readCbModel([myPath filesep 'small_metabolism.mat']);
%model = recon;
model = toy_metabolism
sol = optimizeCbModel(model);sol.f

sampleFile                    =   [];
samplerName                   =  'ACHR';
options.nStepsPerPoint        =  20;%- Number of sampler steps per point saved (200)
options.nPointsReturned       =  500;%- Number of points loaded for analysis (2000)
options.nWarmupPoints         =  10;%. - Number of warmup points (5000). ACHR only.
options.nFiles                =  50;%. - Number of output files (10). ACHR only.
options.nPointsPerFile        =  500;%. - Number of points per file (1000). ACHR only.
%options.nFilesSkipped        =  %. - Number of output files skipped when loading points to avoid potentially biased initial samples (2) loops (true). ACHR only.
options.maxTime               =  36000; %. - Maximum time limit (Default = 36000 s). ACHR only.
%options.toRound              =  %. - Option to round the model before sampling (true). CHRR only.
%options.lambda               =  %. - the bias vector for exponential sampling. CHRR_EXP only.
modelSampling                 =  [];
[modelSampled,samples_leftraru ,volume_leftraru ] = sampleCbModel(model, sampleFile, samplerName, options, modelSampling);

%ACHR_shrunk_leftraru          = modelSampled;

%clear  myPath modelSampled model sampleFile samplerName options modelSampling
scriptPath = '/home/aacevedo/NLHPC/output';
csvwrite([scriptPath filesep 'samples_leftraru.csv'], samples_leftraru)
