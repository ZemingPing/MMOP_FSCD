%% Add path

addpath(genpath('MM_testfunctions/'));
addpath(genpath('Indicator_calculation/'));

clear all
clc

global fname
                                                                                                                                                                                                                                                                                                                                                                                                   
N_function=22;% number of test function
runtimes=1;  % odd number
local_function = [10,11,12,13,15,19]; % the indicator values of these functions are calculated using only global PS and PF
%% Initialize the parameters in MMO test functions
for i_func=1:N_function
    
    switch i_func
        case 1
            fname='MMF1';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[1 -1];     % the low bounds of the decision variables
            xu=[3 1];      % the up bounds of the decision variables
            repoint=[1.1,1.1]; % reference point used to calculate the Hypervolume, it is set to 1.1*(max value of f_i)
        case 2
            fname='MMF2';
            n_obj=2;
            n_var=2;
            xl=[0 0];
            xu=[1 2];
            repoint=[1.1,1.1];
        case 3
            fname='MMF3';
            n_obj=2;
            n_var=2;
            xl=[0 0];
            xu=[1 1.5];
            repoint=[1.1,1.1];
        case 4
            fname='MMF4';
            n_obj=2;
            n_var=2;
            xl=[-1 0];
            xu=[1 2];
            repoint=[1.1,1.1];
        case 5
            fname='MMF5';
            n_obj=2;
            n_var=2;
            xl=[1 -1];
            xu=[3 3];
            repoint=[1.1,1.1];
        case 6
            fname='MMF6';
            n_obj=2;
            n_var=2;
            xl=[1 -1];
            xu=[3 2];
            repoint=[1.1,1.1];
        case 7
            fname='MMF7';
            n_obj=2;
            n_var=2;
            xl=[1 -1];
            xu=[3 1];
            repoint=[1.1,1.1];
        case 8
            fname='MMF8';
            n_obj=2;
            n_var=2;
            xl=[-pi 0];
            xu=[pi 9];
            repoint=[1.1,1.1];
        case 9
            fname='MMF9';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[0.1 0.1];     % the low bounds of the decision variables
            xu=[1.1 1.1];      % the up bounds of the decision variables
            repoint=[1.21,11]; % reference point used to calculate the Hypervolume
        case 10
            fname='MMF10';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[0.1 0.1];     % the low bounds of the decision variables
            xu=[1.1 1.1];      % the up bounds of the decision variables
            repoint=[1.21,13.2]; % reference point used to calculate the Hypervolume
        case 11
            fname='MMF11';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[0.1 0.1];     % the low bounds of the decision variables
            xu=[1.1 1.1];      % the up bounds of the decision variables
            repoint=[1.21,15.4];
        case 12
            fname='MMF12';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[0 0];     % the low bounds of the decision variables
            xu=[1 1];      % the up bounds of the decision variables
            repoint=[1.54,1.1];
        case 13
            %*need to be modified
            fname='MMF13';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=3;       % the dimensions of the objective space
            xl=[0.1 0.1 0.1];     % the low bounds of the decision variables
            xu=[1.1 1.1 1.1];      % the up bounds of the decision variables
            repoint=[1.54,15.4];
        case 14
            fname='MMF14';  % function name
            n_obj=3;       % the dimensions of the decision space
            n_var=3;       % the dimensions of the objective space
            xl=[0 0 0];     % the low bounds of the decision variables
            xu=[1 1 1];      % the up bounds of the decision variables
            repoint=[2.2,2.2,2.2];
        case 15
            fname='MMF15';  % function name
            n_obj=3;       % the dimensions of the decision space
            n_var=3;       % the dimensions of the objective space
            xl=[0 0 0];     % the low bounds of the decision variables
            xu=[1 1 1];      % the up bounds of the decision variables
            repoint=[2.5,2.5,2.5];
        case 16
            fname='MMF1_z';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[1 -1];     % the low bounds of the decision variables
            xu=[3 1];      % the up bounds of the decision variables
            repoint=[1.1,1.1];
        case 17
            fname='MMF1_e';  % function name
            n_obj=2;       % the dimensions of the decision space
            n_var=2;       % the dimensions of the objective space
            xl=[1 -20];     % the low bounds of the decision variables
            xu=[3 20];      % the up bounds of the decision variables
            repoint=[1.1,1.1];
        case 18
            fname='MMF14_a';  % function name
            n_obj=3;
            n_var=3;
            xl=[0 0 0];
            xu=[1 1 1];
            repoint=[2.2,2.2,2.2];
        case 19
            fname='MMF15_a';  % function name
            n_obj=3;
            n_var=3;
            xl=[0 0 0];
            xu=[1 1 1];
            repoint=[2.5,2.5,2.5];
        case 20
            fname='SYM_PART_simple';
            n_obj=2;
            n_var=2;
            xl=[-20 -20];
            xu=[20 20];
            repoint=[4.4,4.4];
        case 21
            fname='SYM_PART_rotated';
            n_obj=2;
            n_var=2;
            xl=[-20 -20];
            xu=[20 20];
            repoint=[4.4,4.4];
        case 22
            fname='Omni_test';
            n_obj=2;
            n_var=3;
            xl=[0 0 0];
            xu=[6 6 6];
            repoint=[4.4,4.4];
    end
    %% Load reference PS and PF data
    
    if sum(i_func == local_function) >0
        load (strcat([fname,'_globalPSPF']))
        load (strcat([fname,'_localPSPF']));
        PS = global_PS;
        PF = global_PF;
    else
        load  (strcat([fname,'_Reference_PSPF_data']));
        load  (strcat([fname,'_Reference_PSPF_data']));
    end
    %% Initialize the population size and the maximum evaluations
    popsize=200;
    Max_fevs=200*100;
    Max_Gen=fix(Max_fevs/popsize);
    Indicator=[];
    
    for j=1:runtimes
        %% Search the PSs using NCDE
        fprintf('Running test function: %s, times = %d \n', fname,j);
        [ps,pf]=MMOP_FSCD(fname,xl,xu,n_obj,popsize,Max_Gen);
        % Indicators
        hyp=Hypervolume_calculation(pf,repoint);
        IGDx=IGD_calculation(ps,PS);
        IGDf=IGD_calculation(pf,PF);
        CR=CR_calculation(ps,PS);
        PSP=CR/IGDx;%
        Indicator.MMOP_FSCD(j,:)=[1./PSP,1./hyp,IGDx,IGDf];
        %Indicator.MMOP_FSCD(j,:)=[1./PSP,hyp,IGDx,IGDf];
        fprintf('Running test function: %s \n %d times rPSP=%f \n', fname,j,1./PSP);
        PSdata.MMOP_FSCD{j}=ps;
        PFdata.MMOP_FSCD{j}=pf;%
        clear ps pf hyp IGDx IGDf CR PSP
        
        %fprintf('Running test function: %s \n %d times \n', fname,j);
    end
    
    %% Choose one PS with the median indicators
    % Choose one PS for MMOP_FSCD
    choosen_In=Indicator.MMOP_FSCD(:,1);% Choose PS according to PSP (the fifth indicator) value
    median_index=find(choosen_In==median(choosen_In));
    choose_ps.MMOP_FSCD = PSdata.MMOP_FSCD{median_index};
    choose_pf.MMOP_FSCD  = PFdata.MMOP_FSCD{median_index};
    clear choosen_In median_index
    
    %% Calculate mean and std of the indicators
    % NCDE
    Indicator.MMOP_FSCD(runtimes+1,:)=min(Indicator.MMOP_FSCD(1:runtimes,:)); %the minimum is the best
    Indicator.MMOP_FSCD(runtimes+2,:)=max(Indicator.MMOP_FSCD(1:runtimes,:)); %the maximum is the worst
    Indicator.MMOP_FSCD(runtimes+3,:)=mean(Indicator.MMOP_FSCD(1:runtimes,:));
    Indicator.MMOP_FSCD(runtimes+4,:)=median(Indicator.MMOP_FSCD(1:runtimes,:));
    Indicator.MMOP_FSCD(runtimes+5,:)=std(Indicator.MMOP_FSCD(1:runtimes,:));
    % Generate Table data in the report
    Table.MMOP_FSCD.rPSP(i_func,:)=(Indicator.MMOP_FSCD(:,1))';%Talbe II data
    Table.MMOP_FSCD.rHV(i_func,:)=(Indicator.MMOP_FSCD(:,2))';%Talbe III data
    Table.MMOP_FSCD.IGDX(i_func,:)=(Indicator.MMOP_FSCD(:,3))';%Talbe IV data
    Table.MMOP_FSCD.IGDF(i_func,:)=(Indicator.MMOP_FSCD(:,4))';%Talbe V data
    
    %% save resultdata
    %save(strcat([fname,'PSPF_indicator_global_data_MMOP_FSCD']),'PSdata','PFdata','Indicator');
    clear PSdata PFdata Indicator
    %save(strcat([fname,'ChoosenPSPFdata']),'choose_ps','choose_pf');

           %% Plot figure
%            if size(choose_ps.MMOP_FSCD,2)==2
%                figure
%                plot(choose_ps.MMOP_FSCD(:,1),choose_ps.MMOP_FSCD(:,2),'o');
%                hold on;
%                plot(PS(:,1),PS(:,2),'r+');
%                legend 'Obtained PS' 'True PS'
           
%            elseif size(choose_ps.MMOP_FSCD,2)==3
%                figure
%                plot3(choose_ps.MMOP_FSCD(:,1),choose_ps.MMOP_FSCD(:,2),choose_ps.MMOP_FSCD(:,3),'o');
%                hold on;
%                plot3(PS(:,1),PS(:,2),PS(:,3),'r+');
%                legend 'Obtained PS' 'True PS'
              
%            end
%            clear choose_ps
end
save Table_global_MMOP_FSCD Table

