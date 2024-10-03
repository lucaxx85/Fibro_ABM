classdef abm_run < handle
    %comment two
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
properties (SetAccess = public)
    Models = {};
    Rule;
    rule;
    generations;  %number of iterations
    GenerationSize = 20; %duration in minutes of each iteration, must also change in InflammatoryDataFitting.m
    Runs = 1;  %number of simulations
    hours = 64; 
    gridSize = 120;% default: 9  (120 corresponds to a grid of 40x40 patches)
    % SA grid size: 9 (3x3), 18 (6x6), 36 (12x12),  72 (24x24) (%added)
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
methods
    
function sim = abm_run()
    %size, runs, generations
    sim.generations = sim.hours*60/sim.GenerationSize;  %tempo in min/durata di ciascuna simulazione = total number of iterations
    for i = 1: sim.Runs
        sim.Models{i} = MacModel_singlemac(sim.gridSize);
        sim.Models{i}.setGeneration(sim.generations);
        sim.Rule = MyRules_singlemac(sim.Models{i});
    end
    sim.generations = sim.Models{1}.MaxGenerations;

    sim.setParameters();
    sim.run();
end

function setParameters(this)
    % young
    for i = 1:length(this.Models)
        % no parameters needed yet
    end 
end

function run(this)
    display('running ...');
    for i = 1:length(this.Models)
        this.Models{i}.run();
        model = this.Models{i};

        % Append rows of data to matrices for each outcome.
        this.Rule.Results{model.Outcome}.Runs = this.Rule.Results{model.Outcome}.Runs + 1;
        this.Rule.Results{model.Outcome}.ImmuneCellCounts = this.Rule.Results{model.Outcome}.ImmuneCellCounts + model.Rules{1}.ImmuneCellCounts;
        this.Rule.Results{model.Outcome}.ImmuneCellCountsSquared = this.Rule.Results{model.Outcome}.ImmuneCellCountsSquared + (model.Rules{1}.ImmuneCellCounts .^ 2);
        this.Rule.Results{model.Outcome}.IntermediateCounts = this.Rule.Results{model.Outcome}.IntermediateCounts + model.Rules{1}.IntermediateCounts;
        this.Rule.Results{model.Outcome}.IntermediateCountsSquared = this.Rule.Results{model.Outcome}.IntermediateCountsSquared + (model.Rules{1}.IntermediateCounts .^ 2);
        this.Rule.Results{model.Outcome}.M1Counts = this.Rule.Results{model.Outcome}.M1Counts + model.Rules{1}.M1Counts;
        this.Rule.Results{model.Outcome}.M1CountsSquared = this.Rule.Results{model.Outcome}.M1CountsSquared + (model.Rules{1}.M1Counts .^ 2);
        this.Rule.Results{model.Outcome}.M2Counts = this.Rule.Results{model.Outcome}.M2Counts + model.Rules{1}.M2Counts;
        this.Rule.Results{model.Outcome}.M2CountsSquared = this.Rule.Results{model.Outcome}.M2CountsSquared + (model.Rules{1}.M2Counts .^ 2);
        this.Rule.Results{model.Outcome}.TotalMacs = this.Rule.Results{model.Outcome}.TotalMacs + model.Rules{1}.TotalMacs;
        this.Rule.Results{model.Outcome}.TotalMacsSquared = this.Rule.Results{model.Outcome}.TotalMacsSquared + (model.Rules{1}.TotalMacs .^ 2);
        this.Rule.Results{model.Outcome}.ProInflammatoryCounts = this.Rule.Results{model.Outcome}.ProInflammatoryCounts + model.Rules{1}.ProInflammatoryCounts;
        this.Rule.Results{model.Outcome}.ProInflammatoryCountsSquared = this.Rule.Results{model.Outcome}.ProInflammatoryCountsSquared + (model.Rules{1}.ProInflammatoryCounts .^ 2);
        this.Rule.Results{model.Outcome}.AntiInflammatoryCounts = this.Rule.Results{model.Outcome}.AntiInflammatoryCounts + model.Rules{1}.AntiInflammatoryCounts;
        this.Rule.Results{model.Outcome}.AntiInflammatoryCountsSquared = this.Rule.Results{model.Outcome}.AntiInflammatoryCountsSquared + (model.Rules{1}.AntiInflammatoryCounts .^ 2);
        this.Rule.Results{model.Outcome}.SOCSCounts = this.Rule.Results{model.Outcome}.SOCSCounts + model.Rules{1}.SOCSCounts;
        this.Rule.Results{model.Outcome}.SOCSCountsSquared = this.Rule.Results{model.Outcome}.SOCSCountsSquared + (model.Rules{1}.SOCSCounts .^ 2);
        this.Rule.Results{model.Outcome}.AverageM1Activation = this.Rule.Results{model.Outcome}.AverageM1Activation + model.Rules{1}.AverageM1Activation;
        this.Rule.Results{model.Outcome}.AverageM1ActivationSquared = this.Rule.Results{model.Outcome}.AverageM1ActivationSquared + (model.Rules{1}.AverageM1Activation .^ 2);
        this.Rule.Results{model.Outcome}.AverageM2Activation = this.Rule.Results{model.Outcome}.AverageM2Activation + model.Rules{1}.AverageM2Activation;
        this.Rule.Results{model.Outcome}.AverageM2ActivationSquared = this.Rule.Results{model.Outcome}.AverageM2ActivationSquared + (model.Rules{1}.AverageM2Activation .^ 2);
        this.Rule.Results{model.Outcome}.RecruitedCells = this.Rule.Results{model.Outcome}.RecruitedCells + model.Rules{1}.RecruitedCells;
        this.Rule.Results{model.Outcome}.RecruitedCellsSquared = this.Rule.Results{model.Outcome}.RecruitedCellsSquared + (model.Rules{1}.RecruitedCells .^ 2);
        this.Rule.Results{model.Outcome}.ProbRecruited = this.Rule.Results{model.Outcome}.ProbRecruited + model.Rules{1}.ProbRecruited;
        this.Rule.Results{model.Outcome}.ProbRecruitedSquared = this.Rule.Results{model.Outcome}.ProbRecruitedSquared + (model.Rules{1}.ProbRecruited .^ 2);
        this.Rule.Results{model.Outcome}.F1Counts = this.Rule.Results{model.Outcome}.F1Counts + model.Rules{1}.F1Counts;
        this.Rule.Results{model.Outcome}.F1CountsSquared = this.Rule.Results{model.Outcome}.F1CountsSquared + (model.Rules{1}.F1Counts .^ 2);
        this.Rule.Results{model.Outcome}.AverageF1Activation = this.Rule.Results{model.Outcome}.AverageF1Activation + model.Rules{1}.AverageF1Activation;
        this.Rule.Results{model.Outcome}.AverageF1ActivationSquared = this.Rule.Results{model.Outcome}.AverageF1ActivationSquared + (model.Rules{1}.AverageF1Activation .^ 2);

        this.Rule.Results{model.Outcome}.RecruitedFibroblasts = this.Rule.Results{model.Outcome}.RecruitedFibroblasts + model.Rules{1}.RecruitedFibroblasts;
        this.Rule.Results{model.Outcome}.RecruitedFibroblastsSquared = this.Rule.Results{model.Outcome}.RecruitedFibroblastsSquared + (model.Rules{1}.RecruitedFibroblasts .^ 2);

        this.Rule.Results{model.Outcome}.ProbRecruited_f = this.Rule.Results{model.Outcome}.ProbRecruited_f + model.Rules{1}.ProbRecruited_f;
        this.Rule.Results{model.Outcome}.ProbRecruitedSquared_f = this.Rule.Results{model.Outcome}.ProbRecruitedSquared_f + (model.Rules{1}.ProbRecruited_f .^ 2);


    end
    
    % calculate averages (among simulations)
        this.Rule.Results{model.Outcome}.ImmuneCellCounts = this.Rule.Results{model.Outcome}.ImmuneCellCounts / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.ImmuneCellCountsSquared = this.Rule.Results{model.Outcome}.ImmuneCellCountsSquared / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.IntermediateCounts = this.Rule.Results{model.Outcome}.IntermediateCounts / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.IntermediateCountsSquared = this.Rule.Results{model.Outcome}.IntermediateCountsSquared / this.Rule.Results{model.Outcome}.Runs;        
        this.Rule.Results{model.Outcome}.M1Counts = this.Rule.Results{model.Outcome}.M1Counts / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.M1CountsSquared = this.Rule.Results{model.Outcome}.M1CountsSquared / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.M2Counts = this.Rule.Results{model.Outcome}.M2Counts / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.M2CountsSquared = this.Rule.Results{model.Outcome}.M2CountsSquared / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.TotalMacs = this.Rule.Results{model.Outcome}.TotalMacs / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.TotalMacsSquared = this.Rule.Results{model.Outcome}.TotalMacsSquared / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.ProInflammatoryCounts = this.Rule.Results{model.Outcome}.ProInflammatoryCounts / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.ProInflammatoryCountsSquared = this.Rule.Results{model.Outcome}.ProInflammatoryCountsSquared / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.AntiInflammatoryCounts = this.Rule.Results{model.Outcome}.AntiInflammatoryCounts / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.AntiInflammatoryCountsSquared = this.Rule.Results{model.Outcome}.AntiInflammatoryCountsSquared / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.SOCSCounts = this.Rule.Results{model.Outcome}.SOCSCounts / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.SOCSCountsSquared = this.Rule.Results{model.Outcome}.SOCSCountsSquared / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.AverageM1Activation = this.Rule.Results{model.Outcome}.AverageM1Activation / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.AverageM1ActivationSquared = this.Rule.Results{model.Outcome}.AverageM1ActivationSquared / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.AverageM2Activation = this.Rule.Results{model.Outcome}.AverageM2Activation / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.AverageM2ActivationSquared = this.Rule.Results{model.Outcome}.AverageM2ActivationSquared / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.RecruitedCells = this.Rule.Results{model.Outcome}.RecruitedCells / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.RecruitedCellsSquared = this.Rule.Results{model.Outcome}.RecruitedCellsSquared / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.ProbRecruited = this.Rule.Results{model.Outcome}.ProbRecruited / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.ProbRecruitedSquared = this.Rule.Results{model.Outcome}.ProbRecruitedSquared / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.F1Counts = this.Rule.Results{model.Outcome}.F1Counts / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.F1CountsSquared = this.Rule.Results{model.Outcome}.F1CountsSquared / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.AverageF1Activation = this.Rule.Results{model.Outcome}.AverageF1Activation / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.AverageFctivationSquared = this.Rule.Results{model.Outcome}.AverageF1ActivationSquared / this.Rule.Results{model.Outcome}.Runs;

        this.Rule.Results{model.Outcome}.RecruitedFibroblasts = this.Rule.Results{model.Outcome}.RecruitedFibroblasts / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.RecruitedFibroblastsSquared = this.Rule.Results{model.Outcome}.RecruitedFibroblastsSquared / this.Rule.Results{model.Outcome}.Runs;

        this.Rule.Results{model.Outcome}.ProbRecruited_f = this.Rule.Results{model.Outcome}.ProbRecruited_f / this.Rule.Results{model.Outcome}.Runs;
        this.Rule.Results{model.Outcome}.ProbRecruitedSquared_f = this.Rule.Results{model.Outcome}.ProbRecruitedSquared_f / this.Rule.Results{model.Outcome}.Runs;

        % convert to vectors
        %%% healthy outcome
        avgm0count_h = [this.Models{i}.InitialImmuneCount (this.Rule.Results{1}.ImmuneCellCounts)];
        avgintcount_h = [0 (this.Rule.Results{1}.IntermediateCounts)];
        avgm1count_h = [this.Models{i}.InitialM1Count (this.Rule.Results{1}.M1Counts)];
        avgm2count_h = [this.Models{i}.InitialM2Count (this.Rule.Results{1}.M2Counts)];
        avgfcount_h = [this.Models{i}.InitialF1Count (this.Rule.Results{1}.F1Counts)];

        avgtotalmacs_h = [avgm0count_h(1)+avgm1count_h(1)+avgm2count_h(1) (this.Rule.Results{1}.TotalMacs)];
        avgpimcount_h = [this.Models{i}.InitTotalPIM (this.Rule.Results{1}.ProInflammatoryCounts)];
        avgaimcount_h = [this.Models{i}.InitTotalAIM (this.Rule.Results{1}.AntiInflammatoryCounts)];
        avgsocscount_h = [this.Models{i}.InitialSOCS (this.Rule.Results{1}.SOCSCounts)];
        avgm1act_h = [mean(mean(this.Models{i}.InitialM1ActivationLattice)) (this.Rule.Results{1}.AverageM1Activation)];
        avgm2act_h = [mean(mean(this.Models{i}.InitialM2ActivationLattice)) (this.Rule.Results{1}.AverageM2Activation)];
        avgfact_h = [mean(mean(this.Models{i}.InitialF1ActivationLattice)) (this.Rule.Results{1}.AverageF1Activation)];

        avgrecruit_h = [0 (this.Rule.Results{1}.RecruitedCells)];
        avgprobrecruit_h = [0 (this.Rule.Results{1}.ProbRecruited)];

        avgrecruit_f_h = [0 (this.Rule.Results{1}.RecruitedFibroblasts)];
        avgprobrecruit_f_h = [0 (this.Rule.Results{1}.ProbRecruited_f)];

        avgm0count_sq_h = [avgm0count_h(1)^2 (this.Rule.Results{1}.ImmuneCellCountsSquared)];
        avgintcount_sq_h = [avgintcount_h(1)^2 (this.Rule.Results{1}.IntermediateCountsSquared)];
        avgm1count_sq_h = [avgm1count_h(1)^2 (this.Rule.Results{1}.M1CountsSquared)];
        avgm2count_sq_h = [avgm2count_h(1)^2 (this.Rule.Results{1}.M2CountsSquared)];
        avgfcount_sq_h = [avgfcount_h(1)^2 (this.Rule.Results{1}.F1CountsSquared)];

        avgtotalmacs_sq_h = [avgtotalmacs_h(1)^2 (this.Rule.Results{1}.TotalMacsSquared)];
        avgpimcount_sq_h = [avgpimcount_h(1)^2 (this.Rule.Results{1}.ProInflammatoryCountsSquared)];
        avgaimcount_sq_h = [avgaimcount_h(1)^2 (this.Rule.Results{1}.AntiInflammatoryCountsSquared)];
        avgsocscount_sq_h = [avgsocscount_h(1)^2 (this.Rule.Results{1}.SOCSCountsSquared)];
        avgm1act_sq_h = [avgm1act_h(1)^2 (this.Rule.Results{1}.AverageM1ActivationSquared)];
        avgm2act_sq_h = [avgm2act_h(1)^2 (this.Rule.Results{1}.AverageM2ActivationSquared)];
       avgfact_sq_h = [avgfact_h(1)^2 (this.Rule.Results{1}.AverageF1ActivationSquared)];

        avgrecruit_sq_h = [avgrecruit_h(1)^2 (this.Rule.Results{1}.RecruitedCellsSquared)];
        avgprobrecruit_sq_h = [avgprobrecruit_h(1)^2 (this.Rule.Results{1}.ProbRecruitedSquared)];
       
        avgrecruit_f_sq_h = [avgrecruit_h(1)^2 (this.Rule.Results{1}.RecruitedFibroblastsSquared)];
        avgprobrecruit_f_sq_h = [avgprobrecruit_f_h(1)^2 (this.Rule.Results{1}.ProbRecruitedSquared_f)];


        %%% inflamed outcome
%         avgm0count_i = [this.Models{i}.InitialImmuneCount (this.Rule.Results{2}.ImmuneCellCounts)];
%         avgintcount_i = [0 (this.Rule.Results{2}.IntermediateCounts)];
%         avgm1count_i = [this.Models{i}.InitialM1Count (this.Rule.Results{2}.M1Counts)];
%         avgm2count_i = [this.Models{i}.InitialM2Count (this.Rule.Results{2}.M2Counts)];
%         avgtotalmacs_i = [avgm0count_i(1)+avgm1count_i(1)+avgm2count_i(1) (this.Rule.Results{2}.TotalMacs)];
%         avgpimcount_i = [this.Models{i}.InitTotalPIM (this.Rule.Results{2}.ProInflammatoryCounts)];
%         avgaimcount_i = [this.Models{i}.InitTotalAIM (this.Rule.Results{2}.AntiInflammatoryCounts)];
%         avgsocscount_i = [this.Models{i}.InitialSOCS (this.Rule.Results{2}.SOCSCounts)];
%         avgm1act_i = [mean(mean(this.Models{i}.InitialM1ActivationLattice)) (this.Rule.Results{2}.AverageM1Activation)];
%         avgm2act_i = [mean(mean(this.Models{i}.InitialM2ActivationLattice)) (this.Rule.Results{2}.AverageM2Activation)];
%         avgrecruit_i = [0 (this.Rule.Results{2}.RecruitedCells)];
%         avgprobrecruit_i = [0 (this.Rule.Results{2}.ProbRecruited)];
%         
%         avgm0count_sq_i = [avgm0count_i(1)^2 (this.Rule.Results{2}.ImmuneCellCountsSquared)];
%         avgintcount_sq_i = [avgintcount_i(1)^2 (this.Rule.Results{2}.IntermediateCountsSquared)];
%         avgm1count_sq_i = [avgm1count_i(1)^2 (this.Rule.Results{2}.M1CountsSquared)];
%         avgm2count_sq_i = [avgm2count_i(1)^2 (this.Rule.Results{2}.M2CountsSquared)];
%         avgtotalmacs_sq_i = [avgtotalmacs_i(1)^2 (this.Rule.Results{2}.TotalMacsSquared)];
%         avgpimcount_sq_i = [avgpimcount_i(1)^2 (this.Rule.Results{2}.ProInflammatoryCountsSquared)];
%         avgaimcount_sq_i = [avgaimcount_i(1)^2 (this.Rule.Results{2}.AntiInflammatoryCountsSquared)];
%         avgsocscount_sq_i = [avgsocscount_i(1)^2 (this.Rule.Results{2}.SOCSCountsSquared)];
%         avgm1act_sq_i = [avgm1act_i(1)^2 (this.Rule.Results{2}.AverageM1ActivationSquared)];
%         avgm2act_sq_i = [avgm2act_i(1)^2 (this.Rule.Results{2}.AverageM2ActivationSquared)];
%         avgrecruit_sq_i = [avgrecruit_i(1)^2 (this.Rule.Results{2}.RecruitedCellsSquared)];
%         avgprobrecruit_sq_i = [avgprobrecruit_i(1)^2 (this.Rule.Results{2}.ProbRecruitedSquared)];
        
% compute standard deviations
%%% healthy outcome
sdm0count_h = sqrt(avgm0count_sq_h - avgm0count_h.^2);
sdintcount_h = sqrt(avgintcount_sq_h - avgintcount_h.^2);
sdm1count_h = sqrt(avgm1count_sq_h - avgm1count_h.^2);
sdm2count_h = sqrt(avgm2count_sq_h - avgm2count_h.^2);
sdfcount_h = sqrt(avgfcount_sq_h - avgfcount_h.^2);

sdtotalmacs_h = sqrt(avgtotalmacs_sq_h - avgtotalmacs_h.^2);
sdpimcount_h = sqrt(avgpimcount_sq_h - avgpimcount_h.^2);
sdaimcount_h = sqrt(avgaimcount_sq_h - avgaimcount_h.^2);
sdsocscount_h = sqrt(avgsocscount_sq_h - avgsocscount_h.^2);
sdm1act_h = sqrt(avgm1act_sq_h - avgm1act_h.^2);
sdm2act_h = sqrt(avgm2act_sq_h - avgm2act_h.^2);
sdfact_h = sqrt(avgfact_sq_h - avgfact_h.^2);

sdrecruit_h = sqrt(avgrecruit_sq_h - avgrecruit_h.^2);
sdprobrecruit_h = sqrt(avgprobrecruit_sq_h - avgprobrecruit_h.^2);
sdrecruit_f_h = sqrt(avgrecruit_f_sq_h - avgrecruit_f_h.^2);
sdprobrecruit_h = sqrt(avgprobrecruit_f_sq_h - avgprobrecruit_f_h.^2);
%%% inflamed outcome
%         sdm0count_i = sqrt(avgm0count_sq_i - avgm0count_i.^2);
%         sdintcount_i = sqrt(avgintcount_sq_i - avgintcount_i.^2);
%         sdm1count_i = sqrt(avgm1count_sq_i - avgm1count_i.^2);
%         sdm2count_i = sqrt(avgm2count_sq_i - avgm2count_i.^2);
%         sdtotalmacs_i = sqrt(avgtotalmacs_sq_i - avgtotalmacs_i.^2);
%         sdpimcount_i = sqrt(avgpimcount_sq_i - avgpimcount_i.^2);
%         sdaimcount_i = sqrt(avgaimcount_sq_i - avgaimcount_i.^2);
%         sdsocscount_i = sqrt(avgsocscount_sq_i - avgsocscount_i.^2);
%         sdm1act_i = sqrt(avgm1act_sq_i - avgm1act_i.^2);
%         sdm2act_i = sqrt(avgm2act_sq_i - avgm2act_i.^2);
%         sdrecruit_i = sqrt(avgrecruit_sq_i - avgrecruit_i.^2);
%         sdprobrecruit_i = sqrt(avgprobrecruit_sq_i - avgprobrecruit_i.^2);
%         
        % save some results for healthy - can do the same for inflamed
        this.Rule.Results{1}.avgm1act = avgm1act_h;
        this.Rule.Results{1}.avgm2act = avgm2act_h;
        this.Rule.Results{1}.avgfact = avgfact_h;
        this.Rule.Results{1}.avgpimcount = avgpimcount_h;
        this.Rule.Results{1}.avgaimcount = avgaimcount_h;
        this.Rule.Results{1}.sdm1act = sdm1act_h;
        this.Rule.Results{1}.sdm2act = sdm2act_h;
                this.Rule.Results{1}.sdfact = sdfact_h;

        this.Rule.Results{1}.sdpimcount = sdpimcount_h;
        this.Rule.Results{1}.sdaimcount = sdaimcount_h;
        %added
        this.Rule.Results{1}.avgm1act_tot = avgm1count_h;
        this.Rule.Results{1}.sdm1count = sdm1count_h;
        this.Rule.Results{1}.avgm2act_tot = avgm2count_h;
        this.Rule.Results{1}.sdm2count = sdm2count_h;
        this.Rule.Results{1}.avgfact_tot = avgfcount_h;
        this.Rule.Results{1}.sdfcount = sdfcount_h;

        this.Rule.Results{1}.avgm0count_tot = avgm0count_h;
        this.Rule.Results{1}.sdm0count = sdm0count_h;
        this.Rule.Results{1}.avgintact_tot = avgintcount_h;
        this.Rule.Results{1}.sdmintcount = sdintcount_h;

        %inflammed (added)
%         this.Rule.Results{2}.avgm1act = avgm1act_i;
%         this.Rule.Results{2}.avgm2act = avgm2act_i;
%         this.Rule.Results{2}.avgpimcount = avgpimcount_i;
%         this.Rule.Results{2}.avgaimcount = avgaimcount_i;
%         this.Rule.Results{2}.sdm1act = sdm1act_i;
%         this.Rule.Results{2}.sdm2act = sdm2act_i;
%         this.Rule.Results{2}.sdpimcount = sdpimcount_i;
%         this.Rule.Results{2}.sdaimcount = sdaimcount_i;


        t = 1:this.generations;       
        t = t .* (this.GenerationSize/60);  %20 min = 1/3 hour
        t=[0 t]; % include initial conditions in plot  %vettore dei tempi in h con time step in h

        % plot results
%         if(model.togglePlot) % full macrophage model  (plots representing results related to all macros)
%             %%% healthy outcome
%             figure('name','Results: Healthy Outcome');          
%             subplot(2,3,1)
%             boundedline(t,avgpimcount_h, sdpimcount_h);ylabel('Pro-inflammatory count');xlabel('hours');
%             title('Healthy outcome')
% 
%             subplot(2,3,2)
%             boundedline(t,avgaimcount_h, sdaimcount_h);ylabel('Anti-inflammatory count');xlabel('hours');
% 
%             subplot(2,3,3)
%             boundedline(t,avgm1act_h, sdm1act_h,'r',t,avgm2act_h, sdm2act_h,'b');ylabel('Average activation');xlabel('hours');
%             legend('M1','M2')
%             
%             subplot(2,3,4)
%             boundedline(t,avgrecruit_h,sdrecruit_h);ylabel('Macrophages recruited');xlabel('hours');
% 
%             subplot(2,3,5)
%             boundedline(t,avgprobrecruit_h,sdprobrecruit_h);ylabel('Average probability of mac recruitment');xlabel('hours');
% 
%             subplot(2,3,6)
%             boundedline(t,avgm1count_h,sdm1count_h,'r',t,avgm2count_h,sdm2count_h,'b',t,avgintcount_h,sdintcount_h,'y');ylabel('Macrophage count');xlabel('hours');
%             legend('M1','M2','Intermediate')
%             
%             
% %             figure;
% %             plot(t,avgtotalmacs_h./1600*100);ylabel('Confluency (%)');xlabel('hours')
%             %%% inflamed outcome
% %             figure('name','Results: Inflamed Outcome');          
% %             subplot(2,3,1)
% %             boundedline(t,avgpimcount_i, sdpimcount_i);ylabel('Pro-inflammatory count');xlabel('hours');
% %             title('Inflamed outcome')
% %             
% %             subplot(2,3,2)
% %             boundedline(t,avgaimcount_i, sdaimcount_i);ylabel('Anti-inflammatory count');xlabel('hours');
% % 
% %             subplot(2,3,3)
% %             boundedline(t,avgm1act_i, sdm1act_i,'r',t,avgm2act_i, sdm2act_i,'b');ylabel('Average activation');xlabel('hours');
% %             legend('M1','M2')
% %             
% %             subplot(2,3,4)
% %             boundedline(t,avgrecruit_i,sdrecruit_i);ylabel('Macrophages recruited');xlabel('hours');
% % 
% %             subplot(2,3,5)
% %             boundedline(t,avgprobrecruit_i,sdprobrecruit_i);ylabel('Average probability of mac recruitment');xlabel('hours');
% % 
% %             subplot(2,3,6)
% %             boundedline(t,avgm1count_i,sdm1count_i,'r',t,avgm2count_i,sdm2count_i,'b',t,avgintcount_i,sdintcount_i,'y');ylabel('Macrophage count');xlabel('hours');
% %             legend('M1','M2','Intermediate')
% %             
% %             figure
% %             boundedline(t,avgm1act_h, sdm1act_h,'r',t,avgm2act_h, sdm2act_h,'b');ylabel('Average activation');xlabel('hours');
% %             legend('M1','M2')
% %             set(gca,'fontsize',16)
%         end
        
        % plot results
%         if(model.togglePlotSingleMac) % single macrophage model
%             %%% healthy outcome
%             figure('name','Results: Healthy Outcome');          
%             subplot(2,3,1)
%             boundedline(t,avgpimcount_h, sdpimcount_h);ylabel('Pro-inflammatory count');xlabel('hours');
%             title('Healthy outcome')
% 
%             subplot(2,3,2)
%             boundedline(t,avgaimcount_h, sdaimcount_h);ylabel('Anti-inflammatory count');xlabel('hours');
% 
%             subplot(2,3,3)
%             boundedline(t,avgm1act_h, sdm1act_h,'r',t,avgm2act_h, sdm2act_h,'b');ylabel('Average activation');xlabel('hours');
%             legend('M1','M2')
%             
%             subplot(2,3,4)
%             boundedline(t,avgsocscount_h,sdsocscount_h);ylabel('Total SOCS');xlabel('hours');
% 
%             subplot(2,3,5) % make sure no macs were recruited
%             boundedline(t,avgrecruit_h,sdrecruit_h);ylabel('Macrophages recruited');xlabel('hours');
%         end
        
        if model.toggleLayeredFigure
            %%% healthy outcome
            figure('name','Results: Healthy Outcome')
            a=area(t,[avgm0count_h; avgm1count_h; avgintcount_h; avgm2count_h; avgfcount_h]');
            a(1).FaceColor = [199 199 199]/255; % M0
            a(2).FaceColor = [255 133 194]/255; % M1
            a(3).FaceColor = [255 253 128]/255; % intermediate
            a(4).FaceColor = [161 176 255]/255; % M2
            a(5).FaceColor = [255 0 0]/255; %F
            xlabel('hours')
            ylabel('Macrophages')
            legend('M0','M1','Intermediate','M2','F')
            set(gca,'fontsize',16)
            %%% inflamed outcome
%             figure('name','Results: Inflamed Outcome')
%             a=area(t,[avgm0count_i; avgm1count_i; avgintcount_i; avgm2count_i]');
%             a(1).FaceColor = [199 199 199]/255; % M0 light grey
%             a(2).FaceColor = [255 133 194]/255; % M1 pink
%             a(3).FaceColor = [255 253 128]/255; % intermediate yellow
%             a(4).FaceColor = [161 176 255]/255; % M2 light blue
%             xlabel('hours')
%             ylabel('Macrophages')
%             legend('M0','M1','Intermediate','M2')
        end
        
        % plot results
%         if(model.SingleMacWriteUpFigures) % single macrophage model
            %%% healthy outcome
%             figure
%             boundedline(t,avgpimcount_h, sdpimcount_h,'r');ylabel('Pro-inflammatory count');xlabel('Time (hours)');
%             set(gca,'fontsize',16)
%             
%             figure
%             boundedline(t,avgaimcount_h, sdaimcount_h,'b');ylabel('Anti-inflammatory count');xlabel('Time (hours)');
%             set(gca,'fontsize',16)
% 
%             figure
%             boundedline(t,avgm1act_h, sdm1act_h,'r');ylabel('Average M1 activation');xlabel('Time (hours)');
%             set(gca,'fontsize',16)
% 
%             figure
%             boundedline(t,avgm2act_h, sdm2act_h,'b');ylabel('Average M2 activation');xlabel('Time (hours)');
%             set(gca,'fontsize',16)
%             
%             figure
%             boundedline(t,avgsocscount_h,sdsocscount_h);ylabel('Fibroblast count');xlabel('Time (hours)');
%             set(gca,'fontsize',16)
%         end

end

end
    
end