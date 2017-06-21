
function [boot_rates, signif_thresholds, signif_cell_inds,signif_cell_totals]= bootstrap_rate_stim_cells (cell, type, threshold,FPS)
%this fxn identifies selective cells based on the difference in rate
%between 2 behavioral states (ie open-closed arms EPM/EZM, center-periphery
%OFT, novelob-oppfield novel object task etc)
%input you need the behavior and cell data, the type of behavior (EPM, EZM,
%OFT, NOVEL), and a signif threshold ('all'= 1SD & 2SD, '1SD', or '2SD'
%only)
%outputs you the rates generated from the bootstrap, the signif thresholds
%per cell, the cell indices to be used in other analyses, and the total
%cell #s

    [cellrates_behstates,states,beh_length]=calc_stim_rates_forboot(cell,type,FPS);
    %pre-allocate bootstrap output sizes
    k=1:size(cell,2);
    cond1_events(k,1000)=zeros;
    cond2_events(k,1000)=zeros;
for i=1:1000%number of shuffle iterations
             cellshuffle=shake(cell,1); %make sure the 'shake' fxn is in your matlab path, shuffles the cell event data
             cond1_events_temp = cellshuffle(states{1},:); 
             %makes a temp matrix w the shuffled cell data w only the condition 1 rows
             cond1_events(:,i)=sum(cond1_events_temp);
             %sums the # of events per cell from that temp matrix
             
             %same as above but for condition2
            cond2_events_temp = cellshuffle(states{2},:);
             cond2_events(:,i)=sum(cond2_events_temp);

end

cond1_rate=cond1_events./beh_length(1); %converts events to rate
cond2_rate=cond2_events./beh_length(2); 

boot_rates={cond1_rate,cond2_rate};
%outputs a cell array w all of the bootstrap rates from above for future
%reference

signif_thresh_cond2_2SD= quantile(cond2_rate,0.95,2); %2SD quantile from bootstrap rate
signif_thresh_cond2_1SD= quantile(cond2_rate,0.68,2); %1SD " "

signif_thresh_cond1_2SD= quantile(cond1_rate,0.95,2);
signif_thresh_cond1_1SD= quantile(cond1_rate,0.68,2);

signif_thresholds=horzcat(signif_thresh_cond1_1SD, signif_thresh_cond2_1SD, signif_thresh_cond1_2SD,signif_thresh_cond2_2SD); 
%outputs the signif threshold levels per cell for future reference

cond2_selective_cells_2SD=cellrates_behstates(:,2)>signif_thresh_cond2_2SD; 
%logical for if cell's diff rate (cond2) exceeds signif threshold
cond2_selective_cells_1SD=cellrates_behstates(:,2)>signif_thresh_cond2_1SD;

cond1_selective_cells_2SD=cellrates_behstates(:,1)>signif_thresh_cond1_2SD;
%logical for if cell's diff rate (cond1) exceeds signif threshold
cond1_selective_cells_1SD=cellrates_behstates(:,1)>signif_thresh_cond1_1SD;

%%%gets cell indices in each condition
cond2_selective_cells2SD=find(cond2_selective_cells_2SD>0);
cond2_selective_cells1SD=find(cond2_selective_cells_1SD>0);

cond1_selective_cells2SD=find(cond1_selective_cells_2SD>0);
cond1_selective_cells1SD=find(cond1_selective_cells_1SD>0);

nonselective_cells2SD=find(cond2_selective_cells_2SD<1 & cond1_selective_cells_2SD<1);
nonselective_cells1SD=find(cond2_selective_cells_1SD<1 & cond1_selective_cells_1SD<1);
%%%

selective_cell_indices_1SD={nonselective_cells1SD, cond1_selective_cells1SD, cond2_selective_cells1SD};
selective_cell_indices_2SD={nonselective_cells2SD, cond1_selective_cells2SD, cond2_selective_cells2SD};

selective_cell_totals_1SD=horzcat(length(nonselective_cells1SD), length(cond1_selective_cells1SD), length(cond2_selective_cells1SD));
selective_cell_totals_2SD=horzcat(length(nonselective_cells2SD), length(cond1_selective_cells2SD), length(cond2_selective_cells2SD));

if strcmpi(threshold,'all')==1
    signif_cell_inds=horzcat(selective_cell_indices_1SD,selective_cell_indices_2SD);
    signif_cell_totals=horzcat(selective_cell_totals_1SD,selective_cell_totals_2SD);
elseif strcmpi(type,'1SD')==1
    signif_cell_inds=selective_cell_indices_1SD;
    signif_cell_totals=selective_cell_totals_1SD;
elseif strcmpi(type,'2SD')==1
    signif_cell_inds=selective_cell_indices_2SD;
    signif_cell_totals=selective_cell_totals_2SD;
end