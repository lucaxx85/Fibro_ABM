classdef MyRules_singlemac < LungModelRule
    

    
properties (SetAccess = public)
    Model;
    GenerationSize = 20; %duration in minutes, must also change in abm_run.m, MacModel.m
    ImmuneCellCounts;  %M0
    ImmuneCellCountsSquared;
    IntermediateCounts;
    IntermediateCountsSquared;
%     %F0Counts;  %fibrocytes
    %F0CountsSquared;
    F1Counts;  %myofibroblasts (F1)
    F1CountsSquared;
    M1Counts;
    M2Counts;
    M1CountsSquared;
    M2CountsSquared;
    TotalMacs;
    TotalMacsSquared;
    %add also Total Fibroblasts (including both F0 and F1)
    TotalFibro;
    TotalFibroSquared;
    ProInflammatoryCounts;   
    ProInflammatoryCountsSquared;
    AntiInflammatoryCounts;
    AntiInflammatoryCountsSquared;   %dobbiamo poterlo osservare anche senza mettere AIM in ingresso, è TGF-beta
    SOCSCounts;
    SOCSCountsSquared;
    AverageM1Activation;
    AverageM1ActivationSquared;
    AverageM2Activation;
    AverageM2ActivationSquared;
    AverageF1Activation;
    AverageF1ActivationSquared;
    RecruitedCells; %macrofagi reclutati 
    RecruitedCellsSquared;
    RecruitedFibroblasts;  %fibrociti (F0), che sono reclutati in base agli AIM prodotti dagli M2 prevalentemente
    RecruitedFibroblastsSquared;
    ProbRecruited;
    ProbRecruitedSquared;
    ProbRecruited_f;
    ProbRecruitedSquared_f;
    Results;
    ProbImmuneCellArrives; %probabilitá nel singolo pixel
    ProbFCellArrives;%probabilitá nel singolo pixel

    % parameters

    %add age for F0 and F1
    AgeMeanM0 = 24; % change here & in MacModel.m
    AgeStDevM0 = 6; % change here & in MacModel.m
    AgeMeanActivated = 48; % change here & in MacModel.m
    AgeStDevActivated = 12; % change here & in MacModel.m
    AgeMeanF0 = 24; %first approximation
    AgeStDevF0=6;%first approximation
    AgeMeanF1 = 48; %first approximation
    AgeStDevF1=12;%first approximation

    %add rate for production of tgf-b and PDGF by F1 (not for now)
    ImmuneProInflammatoryRate = 0.35; % M1s produce pro-inflammatories
    ImmuneAntiInflammatoryRate = 0.85; % M2s produce anti-inflammatories
    ImmuneM1AntiInflammatoryRate=0.175; % M1s produce anti-inflammatories
    ProInflammatoryDecayRate = 0.03;
    AntiInflammatoryDecayRate = 0.03;
    SOCSDecayRate = 0.03;
    FibroAntiInflammatoryRate = 0.85;%myofibroblasts produce anti-inflammatories (first approximation)
    FibroProInflammatoryRate = 0.35; %myofibroblasts produce pro-inflammatories (first approximation) 

    PIMNegativeFeedbackRate=0.002;  %decadimento dell'ativazione 
    AIMNegativeFeedbackRate=0.003;   %decadimento dell'ativazione 
    FNegativeFeedbackRate = 0.003;      %decadimento dell'ativazione (first approximation)

    RecruitmentMMTerm = 30;
    RecruitmentFFTerm = 30;
    AIMRecruitScale=0.1;
    PIMActivationScale=0.75;
    AIMActivationScale=0.75;
    AIMInfinity=1; % Recruitment: M1 activation inhibited by AIM
    
    M1ActivationRate=0.05; % M1 activation increased by PIM (0,1)
    M1ActHillParameter=1; % increase M1 expression via PIM
    M2ActScalar=0.065; % increase M2 activation by AIM
    F1ActScalar=0.065; % increase F activation by AIM
    M2ActHillParameter=0.85; % Hill: increase M2 activation via AIM
    F1ActHillParameter=0.85; % Hill: increase F activation via AIM

    M1AIMInfinity=0.05; % Regulates effectiveness of AIM in inhibiting M1 activation of local cells by PIM
    F1AIMInfinity  = 0.5;  %inhibition of myofibroblasts production of pro-inflammatories by AIM (first approximation)
    M1DecreaseViaAIM=0.01; % AIM decreases M1 activation
    M1DecreaseViaAIMHill=0.4; % Hill parameter - AIM decreases M1 activation

    SOCSProductionRate=4; % rate at which AIM produce SOCS
    AIMSOCSHill=5; % AIM production of SOCS, larger -> more accumulation to be effective
    M1SOCSInfinity=4; % SOCS inhibition of M1 activation
    M2SOCSInfinity=7; % SOCS inhibition of M2 activation
    F1SOCSInfinity = 7; %SOCS inhibition of F1 activation (first approximation)
    AIMSOCSInfinity=0.01; % SOCS inhibition of AIM production
end

methods
    
    function rule = MyRules_singlemac(model)
        rule.Model = model;
        rule.reset();

        rule.Results = cell(1, length(enumeration('Outcomes')));

        % We have to initialize all of the matrices where our counts
        % get appended to. These don't get wiped out between runs.
        % They stick around even after a reset().
        for i = 1:length(enumeration('Outcomes'))
            rule.Results{i}.Runs = 0;
            rule.Results{i}.ImmuneCellCounts = zeros(1, model.MaxGenerations);
            rule.Results{i}.ImmuneCellCountsSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.IntermediateCounts = zeros(1, model.MaxGenerations);
            rule.Results{i}.IntermediateCountsSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.M1Counts = zeros(1, model.MaxGenerations);
            rule.Results{i}.M1CountsSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.M2Counts = zeros(1, model.MaxGenerations);
            rule.Results{i}.M2CountsSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.TotalMacs = zeros(1, model.MaxGenerations);
            rule.Results{i}.TotalMacsSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.ProInflammatoryCounts = zeros(1, model.MaxGenerations);
            rule.Results{i}.ProInflammatoryCountsSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.AntiInflammatoryCounts = zeros(1, model.MaxGenerations);
            rule.Results{i}.AntiInflammatoryCountsSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.SOCSCounts = zeros(1, model.MaxGenerations);
            rule.Results{i}.SOCSCountsSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.AverageM1Activation = zeros(1, model.MaxGenerations);
            rule.Results{i}.AverageM1ActivationSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.AverageM2Activation = zeros(1, model.MaxGenerations);
            rule.Results{i}.AverageM2ActivationSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.RecruitedCells = zeros(1, model.MaxGenerations);
            rule.Results{i}.RecruitedCellsSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.ProbRecruited = zeros(1, model.MaxGenerations);
            rule.Results{i}.ProbRecruitedSquared = zeros(1, model.MaxGenerations);
            %Fibroblasts
            rule.Results{i}.F1Counts = zeros(1, model.MaxGenerations);
            rule.Results{i}.F1CountsSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.AverageF1Activation = zeros(1, model.MaxGenerations);
            rule.Results{i}.AverageF1ActivationSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.RecruitedFibroblasts = zeros(1, model.MaxGenerations);
            rule.Results{i}.RecruitedFibroblastsSquared = zeros(1, model.MaxGenerations);
            rule.Results{i}.ProbRecruited_f = zeros(1, model.MaxGenerations);
            rule.Results{i}.ProbRecruitedSquared_f = zeros(1, model.MaxGenerations);
            rule.Results{i}.ProbFCellArrives = zeros(1, model.MaxGenerations);


        end
    end    

    function apply_rule(this)   
        % move macrophages, determine M1/M2 activation, update age
        model = this.Model;
        if(model.toggleImmune)
            this.moveImmuneCells();
        end
                       
        % diffuse pro- and anti-inflammatory mediators
        pimlattice = model.ProInflammatoryLattice;
        model.ProInflammatoryLattice = diffuse(model.ProInflammatoryLattice); %  aggiungere la chemotassi qui ;
        model.AntiInflammatoryLattice = diffuse(model.AntiInflammatoryLattice);
        %         model.ChemotaxisLattice = model.ProInflammatoryLattice - model.ProInflammatoryLattice(size(model.ProInflammatoryLattice,1)/2,size(model.ProInflammatoryLattice,1)/2);
        %         model.CitokinesLattice = model.ProInflammatoryLattice +  model.ChemotaxisLattice;
        pimdiffuse = model.ProInflammatoryLattice;
        [m,n] = size(pimdiffuse);
        for i = 1:m-2
            for j = 1:n-2
                local_matrix = pimlattice(i:i+2,j:j+2);  %provisional 3x3 matrix with PIM values
                chemotaxis_3 = local_matrix - local_matrix(2,2);
                diffused_3 = pimdiffuse(i:i+2,j:j+2);
                pimcytokines(i:i+2,j:j+2) = chemotaxis_3 + diffused_3;
            end
        end
       model.ProInflammatoryCytokinesLattice = pimcytokines;
        % The amount of pro- and anti-inflammatories in each cell decays at a set rate
        model.ProInflammatoryLattice = model.ProInflammatoryLattice .* (1 - this.ProInflammatoryDecayRate);
        model.AntiInflammatoryLattice = model.AntiInflammatoryLattice .* (1 - this.AntiInflammatoryDecayRate);
        model.SOCSLattice = model.SOCSLattice .* (1 - this.SOCSDecayRate);
        
        
        % recruit macrophages
        if(model.toggleImmune)
            %%% calculate the probability a macrophage will be recruited
            if(model.ToggleRecruitment)
                this.ProbImmuneCellArrives = (model.ProInflammatoryLattice+this.AIMRecruitScale*model.AntiInflammatoryLattice).^2./((model.ProInflammatoryLattice+this.AIMRecruitScale*model.AntiInflammatoryLattice).^2+this.RecruitmentMMTerm^2);
                this.ProbFCellArrives = (model.ProInflammatoryLattice+this.AIMRecruitScale*model.AntiInflammatoryLattice).^2./((model.ProInflammatoryLattice+this.AIMRecruitScale*model.AntiInflammatoryLattice).^2+this.RecruitmentFFTerm^2);
            else
                this.ProbImmuneCellArrives = zeros(size(model.ImmuneLattice));
                this.ProbFCellArrives = zeros(size(model.ImmuneLattice));

            end
            this.ProbRecruited(model.CurrentGeneration)=mean(mean(this.ProbImmuneCellArrives));
            this.ProbRecruited_f(model.CurrentGeneration)=mean(mean(this.ProbFCellArrives));

            temp = (rand(size(model.ImmuneLattice)) < this.ProbImmuneCellArrives) & (~ismember(model.ImmuneLattice,[1 4 6 8 10 12]));
            this.RecruitedCells(model.CurrentGeneration)=sum(sum(temp)); % number of macrophages recruited
            % M1 & M2 activation
            temp_m1act=model.ProInflammatoryLattice./(model.ProInflammatoryLattice+this.PIMActivationScale).*(1/1+(model.AntiInflammatoryLattice/this.AIMInfinity));
            temp_m2act=model.AntiInflammatoryLattice./(model.AntiInflammatoryLattice+this.AIMActivationScale);
            % if M1act+M2act>1, scaling is needed
            if temp_m1act+temp_m2act>1
                temp_sum=temp_m1act+temp_m2act;
                temp_m1act=temp_m1act./temp_sum;
                temp_m2act=temp_m2act./temp_sum;
            end 
            model.M1ActivationLattice(temp) = temp_m1act(temp);
            model.M2ActivationLattice(temp) = temp_m2act(temp);
            temp_m1=temp_m1act>0.5;
            temp_m2=temp_m2act>0.5;
            temp_int=((temp_m1act+temp_m2act)>=0.25) & (temp_m1act<0.5) & (temp_m2act<0.5);
            temp_naive=(temp_m1act+temp_m2act)<0.25;
            % update immune state
            model.ImmuneLattice(temp_naive & temp) = ImmuneStates.M0Static;
            model.ImmuneLattice(temp_m1 & temp) = ImmuneStates.M1Static;
            model.ImmuneLattice(temp_m2 & temp) = ImmuneStates.M2Static;
            model.ImmuneLattice(temp_int & temp) = ImmuneStates.MIntStatic;
            % define ages (naive/activated)
            temp_age_act=ismember(model.ImmuneLattice,[4 6 8]); % activated macrophages
            ages=round(this.AgeStDevActivated.*randn(size(model.ImmuneLattice)) + this.AgeMeanActivated).*60/this.GenerationSize;
            model.ImmuneAge(temp_age_act & temp)=ages(temp_age_act & temp);
            temp_age_m0=model.ImmuneLattice==ImmuneStates.M0Static; 
            ages=round(this.AgeStDevM0.*randn(size(model.ImmuneLattice)) + this.AgeMeanM0).*60/this.GenerationSize;
            model.ImmuneAge(temp_age_m0 & temp)=ages(temp_age_m0 & temp);
           
            %aggiungi qui il reclutamento dei fibroblasti 
            temp = (rand(size(model.ImmuneLattice)) < this.ProbFCellArrives) & (~ismember(model.ImmuneLattice,[1 4 6 8 10 12]));  %10 fibrocytes; 12 myofibroblasts
            this.RecruitedFibroblasts(model.CurrentGeneration)=sum(sum(temp)); % number of fibrocytes recruited
            % F1 activation
            temp_f1act=model.AntiInflammatoryLattice./(model.AntiInflammatoryLattice+this.AIMActivationScale); %the same as m2act, so far
            % if Fact>1, scaling is needed
            if temp_f1act>1
                temp_f1act = 1; %lo taglio a 1
            end 
            model.F1ActivationLattice(temp) = temp_f1act(temp);
            temp_f1=temp_f1act>=0 & temp_f1act<=1;
            temp_f0=temp_f1act<0.25;  %condition over fibrocytes (I chose it)
            % update f state
            model.ImmuneLattice(temp_f0 & temp) = ImmuneStates.F0Static;
            model.ImmuneLattice(temp_f1 & temp) = ImmuneStates.F1Static;
            %restart from here 
            % define ages (naive/activated)

             temp_age_f0=model.ImmuneLattice==ImmuneStates.F0Static; %fybrocytes
            ages_f=round(this.AgeStDevF0.*randn(size(model.ImmuneLattice)) + this.AgeMeanF0).*60/this.GenerationSize;
            model.ImmuneAge(temp_age_f0 & temp)=ages_f(temp_age_f0 & temp);

            temp_age_f1act=ismember(model.ImmuneLattice,12); % activated F
            ages_f=round(this.AgeStDevF1.*randn(size(model.ImmuneLattice)) + this.AgeMeanF1).*60/this.GenerationSize;
            model.ImmuneAge(temp_age_f1act & temp)=ages_f(temp_age_f1act & temp);

            
            clear temp temp_m1act temp_m2act temp_age_act temp_age_m0 ages;
            clear temp temp_f1act temp_age_f1act ages_f;
            %%%% all immune cells recruit & produce inflammatories
            temp=ismember(model.ImmuneLattice,[1 4 6 8 10 12]); % any kind of macrophage and fibroblasts, proportional to activation (see below)
            % M1s produce pro-inflammatories, inhibited by AIM
            model.ProInflammatoryLattice(temp) = model.ProInflammatoryLattice(temp) + this.ImmuneProInflammatoryRate*normrnd(1,0.25,size(model.ImmuneLattice(temp))).*...
                model.M1ActivationLattice(temp).*(1./(1+(model.AntiInflammatoryLattice(temp)./this.M1AIMInfinity)));
            % M1s produce anti-inflammatories
            model.AntiInflammatoryLattice(temp) = model.AntiInflammatoryLattice(temp) + this.ImmuneM1AntiInflammatoryRate*normrnd(1,0.25,size(model.ImmuneLattice(temp))).*...
                model.M1ActivationLattice(temp);
            % M2s produce anti-inflammatories, inhibited by SOCS
            model.AntiInflammatoryLattice(temp) = model.AntiInflammatoryLattice(temp) + this.ImmuneAntiInflammatoryRate*normrnd(1,0.25,size(model.ImmuneLattice(temp))).*...
                model.M2ActivationLattice(temp).*1./(1+(model.SOCSLattice(temp)./this.AIMSOCSInfinity).^3);
            % SOCS are produced by AIM
            model.SOCSLattice(temp) = model.SOCSLattice(temp) + this.SOCSProductionRate*normrnd(1,0.25,size(model.ImmuneLattice(temp))).*...
                model.M2ActivationLattice(temp).^2./(model.M2ActivationLattice(temp)+this.AIMSOCSHill^2);
            % Macrophages which interact with myofibroblasts via CSF1-CSF1R axis ( both M1 and M2 or just M1?) produce anti-inflammatories (like TGF-b and PDGF),
            % inhibited by SOCS (this is not sure)
            model.AntiInflammatoryLattice(temp) = model.AntiInflammatoryLattice(temp) + this.ImmuneM1AntiInflammatoryRate*normrnd(1,0.25,size(model.ImmuneLattice(temp))).*...
                model.F1ActivationLattice(temp).*1./(1+(model.SOCSLattice(temp)./this.F1SOCSInfinity).^3);
            %add here that F1 produce CSF1 acting as PIM (I removed the
            %multiplication factor like in 244-246 lines since here there shouldn't be inhibition
            %by AIM)
            model.ProInflammatoryLattice(temp) = model.ProInflammatoryLattice(temp) + this.FibroProInflammatoryRate*normrnd(1,0.25,size(model.ImmuneLattice(temp))).*...
                model.F1ActivationLattice(temp);
            clear temp;
        end
        
        % Save the plot data
        this.ImmuneCellCounts(model.CurrentGeneration) = sum(model.ImmuneLattice(:) == 1);
        this.ProInflammatoryCounts(model.CurrentGeneration) = sum(model.ProInflammatoryLattice(:));
        this.AntiInflammatoryCounts(model.CurrentGeneration) = sum(model.AntiInflammatoryLattice(:));  
        this.SOCSCounts(model.CurrentGeneration) = sum(model.SOCSLattice(:));  
        this.IntermediateCounts(model.CurrentGeneration) = sum(model.ImmuneLattice(:) == 8);
        this.M1Counts(model.CurrentGeneration) = sum(model.ImmuneLattice(:) == 4);
        this.M2Counts(model.CurrentGeneration) = sum(model.ImmuneLattice(:) == 6);
        this.F1Counts(model.CurrentGeneration) = sum(model.ImmuneLattice(:) == 12);
        % non ho variabili di tipo Counts per i fibrociti
        this.TotalMacs(model.CurrentGeneration) = sum(model.ImmuneLattice(:) == 1) + sum(model.ImmuneLattice(:) == 4) +...
            sum(model.ImmuneLattice(:) == 6) + sum(model.ImmuneLattice(:) == 8);
         this.TotalFibro(model.CurrentGeneration) = sum(model.ImmuneLattice(:) == 10) + sum(model.ImmuneLattice(:) == 12);
        temp=(model.ImmuneLattice(:) == 1) | (model.ImmuneLattice(:) == 4) | (model.ImmuneLattice(:) == 6) | (model.ImmuneLattice(:) == 8);
        this.AverageM1Activation(model.CurrentGeneration) = mean(mean(model.M1ActivationLattice(temp)));
        this.AverageM2Activation(model.CurrentGeneration) = mean(mean(model.M2ActivationLattice(temp)));
        temp=model.ImmuneLattice(:) == 10 | (model.ImmuneLattice(:) == 12);
        this.AverageF1Activation(model.CurrentGeneration) = mean(mean(model.F1ActivationLattice(temp)));

        % Determine the model's outcome. If ??? then the outcome is Inflamed.
        %             if this.ImmuneCellCounts(model.CurrentGeneration) + this.IntermediateCounts(model.CurrentGeneration) + ...
        %                     this.M1Counts(model.CurrentGeneration) + this.M2Counts(model.CurrentGeneration) >= 400
        %                 model.Outcome = Outcomes.Inflamed;
        %             else
        model.Outcome = Outcomes.Healthy;
        %             end
    end

    function reset(this) 
        model = this.Model;

        % Note that the Results DO NOT get reset.
        % Set up the arrays holding the transient plot data.
        this.ImmuneCellCounts = zeros(1, model.MaxGenerations);
        this.ProInflammatoryCounts = zeros(1, model.MaxGenerations);
        this.AntiInflammatoryCounts = zeros(1, model.MaxGenerations);
        this.SOCSCounts = zeros(1, model.MaxGenerations);
        this.IntermediateCounts = zeros(1, model.MaxGenerations);
        this.M1Counts = zeros(1, model.MaxGenerations);
        this.M2Counts = zeros(1, model.MaxGenerations);
        this.TotalMacs = zeros(1, model.MaxGenerations);
        this.AverageM1Activation = zeros(1, model.MaxGenerations);
        this.AverageM2Activation = zeros(1, model.MaxGenerations);
        this.RecruitedCells = zeros(1, model.MaxGenerations);
        this.ProbRecruited = zeros(1, model.MaxGenerations);
        this.F1Counts = zeros(1, model.MaxGenerations);
%         this.F1CountsSquared = zeros(1, model.MaxGenerations);
        this.AverageF1Activation = zeros(1, model.MaxGenerations);
        this.RecruitedFibroblasts = zeros(1, model.MaxGenerations);
        this.ProbRecruited_f = zeros(1, model.MaxGenerations);
%         this.ProbRecruitedSquared_f = zeros(1, model.MaxGenerations);
        this.ProbFCellArrives = zeros(1, model.MaxGenerations);

    end
    
    function moveImmuneCells(this)
        %change this function, temp is useless, leave the for loops and
        %create a probability matrix containing the weight to move
        %macrophages
        model = this.Model;
        %change the space to where the immune cells moves to 3. If the
        %cell cannot move, change the space its on to 3. at the end,
        %change all 3s to 1s.
        temp = randi([1,9],size(model.ImmuneLattice));
        %1 2 3
        %4 5 6
        %7 8 9
        [m,n] = size(model.ImmuneLattice);
        ii=0;
        jj=0;
        for i = 1:m
            for j = 1:n
                if (model.ImmuneLattice(i,j) == ImmuneStates.M0Static) || ...
                        (model.ImmuneLattice(i,j) == ImmuneStates.M1Static) || ...
                        (model.ImmuneLattice(i,j) == ImmuneStates.M2Static) || ...
                        (model.ImmuneLattice(i,j) == ImmuneStates.MIntStatic) || ...
                        (model.ImmuneLattice(i,j) == ImmuneStates.F0Static) || ...   %specified
                        (model.ImmuneLattice(i,j) == ImmuneStates.F1Static)    %added
                    %move up
                    if temp(i,j) < 3
                        ii = i - 1;
                        %move down
                    elseif temp(i,j) > 6
                        ii = i + 1;
                    else
                        ii = i;
                    end
                    %move right
                    if mod(temp(i,j),3) == 0
                        jj = j + 1;
                        %move left
                    elseif mod(temp(i,j),3) == 1
                        jj = j - 1;
                    else
                        jj = j;
                    end
                    %adjust values at edges to spill over to other side
                    if ii <= 0
                        ii = m;
                    end
                    if ii > m
                        ii = 1;
                    end
                    if jj <= 0
                        jj = n;
                    end
                    if jj > n
                        jj = 1;
                    end
                    % move cells
                    if model.ImmuneLattice(ii,jj) ~= ImmuneStates.Empty  %se la destinazione è vuota (cambio con se non è vuoto)
                        ii = i;
                        jj = j;
                    end
                    % move age cell
                    model.ImmuneAge(ii,jj) = model.ImmuneAge(i,j)-1;
                    old_age=model.ImmuneAge(i,j);
                    if i~=ii || j ~= jj
                        model.ImmuneAge(i,j) = 0;
                    end
                    fibroflag = model.ImmuneLattice(i,j) == ImmuneStates.F0Static;
                    if fibroflag ~= 1
                        % move SOCS
                        model.SOCSLattice(ii,jj) = model.SOCSLattice(i,j);
                        if i~=ii || j ~= jj
                            model.SOCSLattice(i,j) = 0;
                        end

                        % update M1/M2 activation
                        pim_fun=@(x) x^2/(x^2+this.M1ActHillParameter^2);
                        aim_fun=@(x,hill) x^2/(x^2+hill^2);

                        % get old M1 activation
                        oldm1act=model.M1ActivationLattice(i,j);
                        model.M1ActivationLattice(i,j)=0; % macrophage no longer there
                        % get old M2 activation
                        oldm2act=model.M2ActivationLattice(i,j);
                        model.M2ActivationLattice(i,j)=0; % macrophage no longer there
                        %get olf F activation
                        %                             oldfact=model.F1ActivationLattice(i,j);
                        %                             model.F1ActivationLattice(i,j)=0; % fibroblast no longer there

                        % increase M1 expression via PIM, inhibited by SOCS
                        model.M1ActivationLattice(ii,jj)=oldm1act+min([this.M1ActivationRate*pim_fun(model.ProInflammatoryLattice(ii,jj))*normrnd(1,0.25)...
                            *1/(1+(model.SOCSLattice(ii,jj)/this.M1SOCSInfinity)^2), 1-oldm1act-oldm2act]);

                        % decrease M1 expression via AIM
                        model.M1ActivationLattice(ii,jj)=max([model.M1ActivationLattice(ii,jj)-this.M1DecreaseViaAIM*aim_fun(model.AntiInflammatoryLattice(ii,jj),this.M1DecreaseViaAIMHill), 0]);

                        % increase M2 expression via AIM
                        model.M2ActivationLattice(ii,jj)=oldm2act+min([this.M2ActScalar*model.AntiInflammatoryLattice(ii,jj)^4/(model.AntiInflammatoryLattice(ii,jj)^4+this.M2ActHillParameter^4)*normrnd(1,0.25)...
                            *1/(1+(model.SOCSLattice(ii,jj)/this.M2SOCSInfinity)^2), 1-oldm1act-oldm2act]);
                        %                         model.M2ActivationLattice(ii,jj)=oldm2act+min([this.M2ActScalar*aim_fun(model.AntiInflammatoryLattice(ii,jj),this.M2ActHillParameter)*normrnd(1,0.25)...

                        % increase F expression via AIM
                        %                             model.F1ActivationLattice(ii,jj)=oldfact+min([this.F1ActScalar*model.AntiInflammatoryLattice(ii,jj)^4/(model.AntiInflammatoryLattice(ii,jj)^4+this.F1ActHillParameter^4)*normrnd(1,0.25)...
                        %                                 *1/(1+(model.SOCSLattice(ii,jj)/this.F1SOCSInfinity)^2), 1-oldfact]);

                        % decrease M1 & M2 expression (natural decay)
                        model.M1ActivationLattice(ii,jj)=model.M1ActivationLattice(ii,jj).*(1-this.PIMNegativeFeedbackRate);
                        model.M2ActivationLattice(ii,jj)=model.M2ActivationLattice(ii,jj).*(1-this.AIMNegativeFeedbackRate);
                        %                             model.F1ActivationLattice(ii,jj)=model.F1ActivationLattice(ii,jj).*(1-this.FNegativeFeedbackRate);

                        % make old space empty
                        %                         fibroflag = model.ImmuneLattice(i,j) == ImmuneStates.F1Static;
                        model.ImmuneLattice(i,j) = ImmuneStates.Empty;

                        % change new state

                        % was original state M0?
                        oldstate=model.ImmuneLattice(ii,jj) == ImmuneStates.M0Moving;
                        %                             if fibroflag~=1
                        if model.M1ActivationLattice(ii,jj)>0.5
                            model.ImmuneLattice(ii,jj) = ImmuneStates.M1Moving;
                            if oldstate==1 % if M0 -> M1, change age to 12 hours
                                model.ImmuneAge(ii,jj)=min(old_age,(round(this.AgeStDevActivated.*randn(1,1) + this.AgeMeanActivated).*60/this.GenerationSize));
                            end
                        elseif model.M2ActivationLattice(ii,jj)>0.5
                            model.ImmuneLattice(ii,jj) = ImmuneStates.M2Moving;
                            if oldstate==1 % if M0 -> M2, change age to 12 hours
                                model.ImmuneAge(ii,jj)=min(old_age,(round(this.AgeStDevActivated.*randn(1,1) + this.AgeMeanActivated).*60/this.GenerationSize));
                            end
                        elseif model.M1ActivationLattice(ii,jj)+model.M2ActivationLattice(ii,jj)>0.25
                            model.ImmuneLattice(ii,jj) = ImmuneStates.MIntMoving;
                            if oldstate==1 % if M0 -> intermediate, change age to 12 hours
                                model.ImmuneAge(ii,jj)=min(old_age,(round(this.AgeStDevActivated.*randn(1,1) + this.AgeMeanActivated).*60/this.GenerationSize));
                            end
                        else
                            model.ImmuneLattice(ii,jj) = ImmuneStates.M0Moving;
                        end
                       
                    elseif fibroflag == 1
                        oldf1act=model.F1ActivationLattice(i,j);
                        model.F1ActivationLattice(i,j)=0; % fibroblast no longer there
                        model.ImmuneLattice(ii,jj) = ImmuneStates.F1Moving;
                        %increase of F1activation when myofibroblasts are
                        %obtained through FMT
                        model.F1ActivationLattice(ii,jj)=oldf1act+min([this.F1ActScalar*model.AntiInflammatoryLattice(ii,jj)^4/(model.AntiInflammatoryLattice(ii,jj)^4+this.F1ActHillParameter^4)*normrnd(1,0.25)...
                            *1/(1+(model.SOCSLattice(ii,jj)/this.F1SOCSInfinity)^2), 1-oldf1act]);

                        %natural decay (not sure to include it)
                        %model.F1ActivationLattice(ii,jj)=model.F1ActivationLattice(ii,jj).*(1-this.FNegativeFeedbackRate);

                        model.ImmuneLattice(i,j) = ImmuneStates.Empty;
                        

                    end
                    %                     end
                end
            end
        end
        %finalize cells from Moving to Static
        % M0
        temp = model.ImmuneLattice == ImmuneStates.M0Moving;
        model.ImmuneLattice(temp) = ImmuneStates.M0Static;
        % M1
        temp = model.ImmuneLattice == ImmuneStates.M1Moving;
        model.ImmuneLattice(temp) = ImmuneStates.M1Static;
        % M2
        temp = model.ImmuneLattice == ImmuneStates.M2Moving;
        model.ImmuneLattice(temp) = ImmuneStates.M2Static;
        % intermediate
        temp = model.ImmuneLattice == ImmuneStates.MIntMoving;
        model.ImmuneLattice(temp) = ImmuneStates.MIntStatic;
        % F0 (fybrocytes)
        temp = model.ImmuneLattice == ImmuneStates.F0Moving;
        model.ImmuneLattice(temp) = ImmuneStates.F0Static;
        %fibroblasts
        temp = model.ImmuneLattice == ImmuneStates.F1Moving;
        model.ImmuneLattice(temp) = ImmuneStates.F1Static;

        %kill cells
        tempDeath = (model.ImmuneAge <= 0);
        model.ImmuneLattice(tempDeath) = ImmuneStates.Empty;
        model.M1ActivationLattice(tempDeath) = 0;
        model.M2ActivationLattice(tempDeath) = 0;
        model.SOCSLattice(tempDeath) = 0;
        model.ImmuneAge(tempDeath) = 0;
        model.F1ActivationLattice(tempDeath)=0;
    end
    
    function finalize(this)
        model = this.Model;
        % Append rows of data to matrices for each outcome.
        this.Results{model.Outcome}.Runs = this.Results{model.Outcome}.Runs + 1;

        this.Results{model.Outcome}.ImmuneCellCounts = this.Results{model.Outcome}.ImmuneCellCounts + model.Rules{1}.Results{model.Outcome}.ImmuneCellCounts;
        this.Results{model.Outcome}.ImmuneCellCountsSquared = this.Results{model.Outcome}.ImmuneCellCountsSquared + (model.Rules{1}.Results{model.Outcome}.ImmuneCellCounts .^ 2);
        this.Results{model.Outcome}.IntermediateCounts = this.Results{model.Outcome}.IntermediateCounts + model.Rules{1}.Results{model.Outcome}.IntermediateCounts;
        this.Results{model.Outcome}.IntermediateCountsSquared = this.Results{model.Outcome}.IntermediateCountsSquared + (model.Rules{1}.Results{model.Outcome}.IntermediateCounts .^ 2);        
        this.Results{model.Outcome}.M1Counts = this.Results{model.Outcome}.M1Counts + model.Rules{1}.Results{model.Outcome}.M1Counts;
        this.Results{model.Outcome}.M1CountsSquared = this.Results{model.Outcome}.M1CountsSquared + (model.Rules{1}.Results{model.Outcome}.M1Counts .^ 2);
        this.Results{model.Outcome}.M2Counts = this.Results{model.Outcome}.M2Counts + model.Rules{1}.Results{model.Outcome}.M2Counts;
        this.Results{model.Outcome}.M2CountsSquared = this.Results{model.Outcome}.M2CountsSquared + (model.Rules{1}.Results{model.Outcome}.M2Counts .^ 2);
        this.Results{model.Outcome}.TotalMacs = this.Results{model.Outcome}.TotalMacs + model.Rules{1}.Results{model.Outcome}.TotalMacs;
        this.Results{model.Outcome}.TotalMacsSquared = this.Results{model.Outcome}.TotalMacsSquared + (model.Rules{1}.Results{model.Outcome}.TotalMacs .^ 2);        
        this.Results{model.Outcome}.ProInflammatoryCounts = this.Results{model.Outcome}.ProInflammatoryCounts + model.Rules{1}.Results{model.Outcome}.ProInflammatoryCounts;
        this.Results{model.Outcome}.ProInflammatoryCountsSquared = this.Results{model.Outcome}.ProInflammatoryCountsSquared + (model.Rules{1}.Results{model.Outcome}.ProInflammatoryCounts .^ 2);
        this.Results{model.Outcome}.AntiInflammatoryCounts = this.Results{model.Outcome}.AntiInflammatoryCounts + model.Rules{1}.Results{model.Outcome}.AntiInflammatoryCounts;
        this.Results{model.Outcome}.AntiInflammatoryCountsSquared = this.Results{model.Outcome}.AntiInflammatoryCountsSquared + (model.Rules{1}.Results{model.Outcome}.AntiInflammatoryCounts .^ 2);
        this.Results{model.Outcome}.SOCSCounts = this.Results{model.Outcome}.SOCSCounts + model.Rules{1}.Results{model.Outcome}.SOCSCounts;
        this.Results{model.Outcome}.SOCSCountsSquared = this.Results{model.Outcome}.SOCSCountsSquared + (model.Rules{1}.Results{model.Outcome}.SOCSCounts .^ 2);
        this.Results{model.Outcome}.AverageM1Activation = this.Results{model.Outcome}.AverageM1Activation + model.Rules{1}.Results{model.Outcome}.AverageM1Activation;
        this.Results{model.Outcome}.AverageM1ActivationSquared = this.Results{model.Outcome}.AverageM1ActivationSquared + (model.Rules{1}.Results{model.Outcome}.AverageM1Activation.^ 2);
        this.Results{model.Outcome}.AverageM2Activation = this.Results{model.Outcome}.AverageM2Activation + model.Rules{1}.Results{model.Outcome}.AverageM2Activation;
        this.Results{model.Outcome}.AverageM2ActivationSquared = this.Results{model.Outcome}.AverageM2ActivationSquared + (model.Rules{1}.Results{model.Outcome}.AverageM2Activation.^ 2);
        this.Results{model.Outcome}.RecruitedCells = this.Results{model.Outcome}.RecruitedCells + model.Rules{1}.Results{model.Outcome}.RecruitedCells;
        this.Results{model.Outcome}.RecruitedCellsSquared = this.Results{model.Outcome}.RecruitedCellsSquared + (model.Rules{1}.Results{model.Outcome}.RecruitedCells.^ 2);
        this.Results{model.Outcome}.ProbRecruited = this.Results{model.Outcome}.ProbRecruited + model.Rules{1}.Results{model.Outcome}.ProbRecruited;
        this.Results{model.Outcome}.ProbRecruitedSquared = this.Results{model.Outcome}.ProbRecruitedSquared + (model.Rules{1}.Results{model.Outcome}.ProbRecruited.^ 2);

        this.Results{model.Outcome}.F1Counts = this.Results{model.Outcome}.F1Counts + model.Rules{1}.Results{model.Outcome}.F1Counts;
        this.Results{model.Outcome}.F1CountsSquared = this.Results{model.Outcome}.F1CountsSquared + (model.Rules{1}.Results{model.Outcome}.F1Counts.^ 2);
        this.Results{model.Outcome}.AverageF1Activation = this.Results{model.Outcome}.AverageF1Activation + model.Rules{1}.Results{model.Outcome}.AverageF1Activation;
        this.Results{model.Outcome}.AverageF1ActivationSquared = this.Results{model.Outcome}.AverageF1ActivationSquared + (model.Rules{1}.Results{model.Outcome}.AverageF1Activation.^ 2);
        this.Results{model.Outcome}.RecruitedFibroblasts = this.Results{model.Outcome}.RecruitedFibroblasts + model.Rules{1}.Results{model.Outcome}.RecruitedFibroblasts;
        this.Results{model.Outcome}.RecruitedFibroblasts = this.Results{model.Outcome}.RecruitedFibroblastsSquared + (model.Rules{1}.Results{model.Outcome}.RecruitedFibroblasts.^ 2);
        this.Results{model.Outcome}.ProbRecruited_f = this.Results{model.Outcome}.ProbRecruited_f + model.Rules{1}.Results{model.Outcome}.ProbRecruited_f;
        this.Results{model.Outcome}.ProbRecruitedSquared_f = this.Results{model.Outcome}.ProbRecruitedSquared_f + (model.Rules{1}.Results{model.Outcome}.ProbRecruited_f.^ 2);

    end
    
    
end
    

    
end