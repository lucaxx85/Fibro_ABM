
classdef MacModel_singlemac < handle

    properties (SetAccess = public)
        CurrentGeneration = 0;
        PreviousGeneration = 0;
        CurrentCounter = 0;
        MaxGenerations;
        GenerationSize = 20; %duration in minutes, must also change in abm_run, MyRules
        ImmuneLattice; %qui dentro aggiungiamo i fibroblasti
        ProInflammatoryLattice;
%         ProInflammatoryLattice_fixed;
        AntiInflammatoryLattice;
        SOCSLattice;
        M1ActivationLattice;
        M2ActivationLattice;
        F1ActivationLattice;
        InitialM1ActivationLattice;
        InitialM2ActivationLattice;
        InitialF1ActivationLattice;
%         ImmuneCellCount;
        InitialImmuneAge;
        ImmuneAge;  %è l'etá di tutto il lattice (M + F)
        InitialImmuneMatrix;
        InitialImmuneCount=0;  %numero iniziale di macrofagi (default 1) sono in M0
        InitialM1Count=0;  %numero di macrofagi che si trovano dall'inizio in M1
        InitialM2Count=0;   %numero di macrofagi che si trovano dall'inizio in M2
        InitialIntCount=0;    %numero di macrofagi che si trovano dall'inizio in M intermedio
        InitialF0Count=0;
        InitialF1Count=0;
        InitialPIM;   %Pro inflammatory stimulus (a.u.)  %10 usato da Luca
        InitTotalPIM = 0;  %somma sui pixel di initialPIM
        InitTotalPIM_fixed = 0;  %somma sui pixel di initialPIM
        InitialAIM = 0;
        InitTotalAIM = 0;
        InitialSOCS = 0;
        AgeMeanM0 = 24;   
        AgeStDevM0 = 6;
        AgeMeanActivated = 48; 
        AgeStDevActivated = 12;
        AgeMeanF0=24*7;  %1368
        AgeStDevF0=6*7;
        AgeMeanF1=336;  %1368
        AgeStDevF1=84;   %72
        Outcome = Outcomes.Healthy;
        ShowLattices = true;  %false
        toggleImmune = true;
        ToggleRecruitment = true;  
        ToggleRecruitment_f = true; %in vivo sempre true
        togglePlot = true;  %false
        togglePlotSingleMac = true;  %false
        toggleLayeredFigure = true;   %false
        SingleMacWriteUpFigures = true; %false
        InitialMatrix;
        Rules;
        contatore = 0;
        RuleSet = {'MyRules_singlemac'};
        Debug = false;
        %                     M0          blank   problem  M1           problem  M2           problem  intermediate   problem
        ImmuneColorMap = [255 255 255; 0 0 0; 0 255 0; 255 0 255; 0 255 0; 0 0 255; 0 255 0; 255 255 0; 0 255 0; 255 0 0; 0 255 0; 128 0 0] / 255;
%            ImmuneColorMap = [255 255 255; 0 0 0; 255 133 194;  161 176 255; 255 253 128; 0 255 0; 128 0 0] / 255;

    end


    methods

        function model = MacModel_singlemac(size)
            % Create macrophage and fibroblast matrix
            
            model.RandomizeMatrix(size);
            model.InitialMatrix = ones(size,size);
            % Create all the rule objects.
            for i = 1:length(model.RuleSet)
                model.Rules{i} = feval(str2func(model.RuleSet{i}), model);
            end
            model.reset();
        end

        function RandomizeMatrix(this,size)
            this.InitialImmuneMatrix = ones(size/3,size/3);
            this.InitialImmuneAge = ones(size/3,size/3);
%             this.InitialFibroAge = ones(size/3,size/3);

            this.InitialImmuneMatrix = this.InitialImmuneMatrix.*2;
            this.InitialM1ActivationLattice = zeros(size/3,size/3);
            this.InitialM2ActivationLattice = zeros(size/3,size/3);
            this.InitialF1ActivationLattice = zeros(size/3,size/3);

            
            %Add random initial immune cells
            n_m0=0;
            n_m1=0;
            n_m2=0;
            n_mint=0;
            n_f0 = 0;
            n_f1 = 0;

            % set up M0
            while n_m0 < this.InitialImmuneCount
                i = randi([1 size/3]);
                j = randi([1 size/3]);
                this.InitialImmuneMatrix(i,j) = 1;
                n_m0 = n_m0 + 1;
                this.InitialImmuneAge(i,j) = (round(this.AgeStDevM0*randn(1)) + this.AgeMeanM0)*60/this.GenerationSize;
                this.InitialM1ActivationLattice(i,j)=random_in_range(0,0.25);
                this.InitialM2ActivationLattice(i,j)=random_in_range(0,0.25-this.InitialM1ActivationLattice(i,j));
            end

            % set up M1
            while n_m1 < this.InitialM1Count
                i = randi([1 size/3]);
                j = randi([1 size/3]);
                if this.InitialImmuneMatrix(i,j)~=1
                    this.InitialImmuneMatrix(i,j) = 4;
                    n_m1 = n_m1 + 1;
                    this.InitialImmuneAge(i,j) = (round(this.AgeStDevActivated*randn(1)) + this.AgeMeanActivated)*60/this.GenerationSize;
                    this.InitialM1ActivationLattice(i,j)=random_in_range(0.5,1);
                    this.InitialM2ActivationLattice(i,j)=random_in_range(0,0.49);
                end
            end

            % set up M2
            while n_m2 < this.InitialM2Count
                i = randi([1 size/3]);
                j = randi([1 size/3]);
                if this.InitialImmuneMatrix(i,j)~=1 || this.InitialImmuneMatrix(i,j)~=4
                    this.InitialImmuneMatrix(i,j) = 6;
                    n_m2 = n_m2 + 1;
                    this.InitialImmuneAge(i,j) = (round(this.AgeStDevActivated*randn(1)) + this.AgeMeanActivated)*60/this.GenerationSize;
                    this.InitialM2ActivationLattice(i,j)=random_in_range(0.5,1);
                    this.InitialM1ActivationLattice(i,j)=random_in_range(0,0.49);
                end
            end
            % set up intermediate  (??)
            while n_mint < this.InitialIntCount
                i = randi([1 size/3]);
                j = randi([1 size/3]);
                if this.InitialImmuneMatrix(i,j)~=1 || this.InitialImmuneMatrix(i,j)~=4 || this.InitialImmuneMatrix(i,j)~=6
                    this.InitialImmuneMatrix(i,j) = 8;
                    n_mint = n_mint + 1;
                    this.InitialImmuneAge(i,j) = (round(this.AgeStDevActivated*randn(1)) + this.AgeMeanActivated)*60/this.GenerationSize;
                    %                     this.ImmuneAge(i,j) = (round(this.AgeStDevActivated*randn(1)) + this.AgeMeanActivated)*60/this.GenerationSize;
                    this.InitialM1ActivationLattice(i,j)=random_in_range(0,0.49);
                    this.InitialM2ActivationLattice(i,j)=random_in_range(max([0, 0.25-this.InitialM1ActivationLattice(i,j)]),0.49);
                end
            end

            while n_f0 < this.InitialF0Count
                i = randi([1 size/3]);
                j = randi([1 size/3]);
                if this.InitialImmuneMatrix(i,j)~=1 || this.InitialImmuneMatrix(i,j)~=4 || this.InitialImmuneMatrix(i,j)~=6 || this.InitialImmuneMatrix(i,j)~=8
                    this.InitialImmuneMatrix(i,j) = 10;
                    n_f0 = n_f0 + 1;
                    this.InitialImmuneAge(i,j) = (round(this.AgeStDevF0*randn(1)) + this.AgeMeanF0)*60/this.GenerationSize;
                    this.InitialF1ActivationLattice(i,j)=random_in_range(0,0.25); %attivazione del F
                    %                     this.InitialM2ActivationLattice(i,j)=random_in_range(max([0, 0.25-this.InitialM1ActivationLattice(i,j)]),0.49);
                end
            end

            while n_f1 < this.InitialF1Count
                i = randi([1 size/3]);
                j = randi([1 size/3]);
                if this.InitialImmuneMatrix(i,j)~=1 || this.InitialImmuneMatrix(i,j)~=4 || this.InitialImmuneMatrix(i,j)~=6 || this.InitialImmuneMatrix(i,j)~=8 || this.InitialImmuneMatrix(i,j)~=10
                    this.InitialImmuneMatrix(i,j) = 12;
                    n_f1 = n_f1 + 1;
                    this.InitialImmuneAge(i,j) = (round(this.AgeStDevF1*randn(1)) + this.AgeMeanF1)*60/this.GenerationSize;
                    this.InitialF1ActivationLattice(i,j)=random_in_range(0.25,1); %attivazione del F
                    %                     this.InitialM2ActivationLattice(i,j)=random_in_range(max([0, 0.25-this.InitialM1ActivationLattice(i,j)]),0.49);
                end
            end

            this.ImmuneLattice = this.InitialImmuneMatrix;
        end

        function setGeneration(this, n)
            this.MaxGenerations = n;
        end

        function reset(this)
           load chemical_pim_gaussian_distribution.mat
            this.ImmuneLattice = this.InitialImmuneMatrix;
            this.ImmuneAge = this.InitialImmuneAge;
            this.ProInflammatoryLattice = zeros(size(this.InitialImmuneMatrix));
            [n,~]=size(this.InitialImmuneMatrix);
            this.ProInflammatoryLattice(((n/3)+1):(2*n/3),((n/3)+1):(2*n/3))= F_norm;  %location of PIM stimulus
            this.InitTotalPIM=sum(sum(this.ProInflammatoryLattice));
            this.AntiInflammatoryLattice = zeros(size(this.InitialImmuneMatrix));
            this.AntiInflammatoryLattice(((n/3)+1):(2*n/3),((n/3)+1):(2*n/3)) = this.InitialAIM;
            this.InitTotalAIM=sum(sum(this.AntiInflammatoryLattice));
            this.SOCSLattice = zeros(size(this.InitialImmuneMatrix));
            this.SOCSLattice(((n/3)+1):(2*n/3),((n/3)+1):(2*n/3)) = this.InitialSOCS;
            this.M1ActivationLattice = this.InitialM1ActivationLattice;
            this.M2ActivationLattice = this.InitialM2ActivationLattice;
            this.F1ActivationLattice = this.InitialF1ActivationLattice;

            this.CurrentGeneration = 0;
            this.PreviousGeneration = 0;
            this.CurrentCounter = 0;
            
            for i = 1:length(this.Rules)
                rule = this.Rules{i};
                rule.reset();
            end
        end

        function run(this)
            % Run the simulation using the rules in the model's RuleSet
            % parameter.

            % Splitting up this figure into separate figures since subimage
            % requires the Image Processing Toolkit.

            if this.ShowLattices
                immune_h = figure;
                set(immune_h, 'Name', 'Macrophages anf fibroblasts');
                title('Immune Cells')


                proinflammatory_h = figure;
                set(proinflammatory_h, 'Name', 'Pro-inflammatory Mediators');
                title('Pro-inflammatory Cells')

%                 proinflammatory_h_fix = figure;
%                 set(proinflammatory_h_fix, 'Name', 'Pro-inflammatory Fixed Stimulus');
%                 title('Biomaterial')


                antiinflammatory_h = figure;
                set(antiinflammatory_h, 'Name', 'Anti-inflammatory Mediators');
                title('Anti-inflammatory Cells')

                m1activation_h = figure;
                set(m1activation_h, 'Name', 'M1 Activation');
                title('M1 Activation')

                m2activation_h = figure;
                set(m2activation_h, 'Name', 'M2 Activation');
                title('M2 Activation')


                factivation_h = figure;
                set(factivation_h, 'Name', 'F Activation');
                title('F Activation')

                %
                %         socs_h = figure;
                %         set(socs_h, 'Name', 'M2 Activation');
                %         title('SOCS')
            end

            if this.Debug
                pause on;
            else
                pause off;
            end

            for generation = 1:this.MaxGenerations
                this.CurrentGeneration = generation;
                this.PreviousGeneration = generation-1;
                this.CurrentCounter = this.contatore;

                if this.ShowLattices
                    % Update the immune cell lattice window
                    set(0, 'CurrentFigure', immune_h);
                    image(this.ImmuneLattice);
                    colormap(this.ImmuneColorMap);
                    axis square;
                    set(gca, 'XTick', [], ...
                        'YTick', [], ...
                        'XTickLabel', '', ...
                        'YTickLabel', '','FontSize',16);
                    title(sprintf('Total time (Hours) %f', generation/(60/this.GenerationSize))); %MUST CHANGE IF YOU CHANGE GENERATION SIZE
                    colorbar('TickLabels',{'M0','M1','M2','Interm.','F0','F1'},'FontSize',14)
%                     savefig('Frame_M0_400_30h_patho_k_2.fig')
%                     saveas(gca,'Frame_M0_400_30h_patho_k_2','epsc');


                    set(0, 'CurrentFigure', proinflammatory_h);
                    imagesc(this.ProInflammatoryLattice);
                    colormap(autumn);
                    axis square;
                    set(gca, 'XTick', [], ...
                        'YTick', [], ...
                        'XTickLabel', '', ...
                        'YTickLabel', '');
                    title('PIM')

%                     set(0, 'CurrentFigure', proinflammatory_h_fix);
%                     imagesc(this.ProInflammatoryLattice_fixed,[0 4]);
%                     colormap(autumn);
%                     axis square;
%                     set(gca, 'XTick', [], ...
%                         'YTick', [], ...
%                         'XTickLabel', '', ...
%                         'YTickLabel', '');
%                     title('fixed PIM')

                    set(0, 'CurrentFigure', antiinflammatory_h);
                    imagesc(this.AntiInflammatoryLattice,[0 10]);
                    colormap(cool);
                    axis square;
                    set(gca, 'XTick', [], ...
                        'YTick', [], ...
                        'XTickLabel', '', ...
                        'YTickLabel', '');
                    title('AIM')

                    %             set(0, 'CurrentFigure', socs_h);
                    %             imagesc(this.SOCSLattice,[0 4]);
                    %             colormap(parula);
                    %             axis square;
                    %             set(gca, 'XTick', [], ...
                    %                      'YTick', [], ...
                    %                      'XTickLabel', '', ...
                    %                      'YTickLabel', '');

                    set(0, 'CurrentFigure', m1activation_h);
                    imagesc(this.M1ActivationLattice,[0 1]);
                    colormap(autumn);
                    axis square;
                    set(gca, 'XTick', [], ...
                        'YTick', [], ...
                        'XTickLabel', '', ...
                        'YTickLabel', '');
                    title('M1 activation')

                    set(0, 'CurrentFigure', m2activation_h);
                    imagesc(this.M2ActivationLattice,[0 1]);
                    colormap(cool);
                    axis square;
                    set(gca, 'XTick', [], ...
                        'YTick', [], ...
                        'XTickLabel', '', ...
                        'YTickLabel', '');
                    title('M2 activation')



                    set(0, 'CurrentFigure', factivation_h);
                    imagesc(this.F1ActivationLattice,[0 1]);
                    colormap(jet); %change color
                    axis square;
                    set(gca, 'XTick', [], ...
                        'YTick', [], ...
                        'XTickLabel', '', ...
                        'YTickLabel', '');
                    title('F activation')

                    drawnow;

                    if this.Debug
                        pause;
                    end
                end

                for i = 1:length(this.Rules)
                    rule = this.Rules{i};
                    rule.apply_rule();
                end

            end

            % After the run we call finalize() on each rule.
            for i = 1:length(this.Rules)
                rule = this.Rules{i};
                rule.finalize();
            end
        end

    end

end