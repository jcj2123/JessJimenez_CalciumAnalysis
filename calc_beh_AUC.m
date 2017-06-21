
function[aveAUC_behstates,numbtransients_behstates]=calc_beh_AUC(behavior, cell_AUC, cell_events, beh_type)
%outputs rates and events per cell per behavior state specified
%calls on the 'define_beh_states' function to define behavior columns for
%diff imaging sessions
%inputs behavior data, cell data, and specify what imaging session as
%'type'
%'type': (but see more behavior column detail in 'define_beh_states'
%'EPM'=elevated plus
%'EPMC'=EPM w center and open compartments counted as one
%'EZM'= elevated zero maze
%'OFT'= open field test
%'NOVEL'= novel object task

[states,numb_beh,~]= define_beh_states(behavior,beh_type);
%uses 'define_beh_states' function to define behavior columns to analyze
%outputs the behavior indices in a cell array "states", how many behavior
%states are extracted from that imaging session (numb_beh), and the length of time the
%animal spent in that behavior state to use for rate calculations
%(beh_length)

%specificy if you want all the rates from 'define_beh_states' function, or
%only rates from the 2 columns that are used in the bootstrap script
%'all' gives all rates
%'bootstrap' outputs the 2 rate columns used for bootstrap plus a
%difference rate column 3
cell_events=cell_events>0;
k=1:size(cell_AUC,2);
aveAUC_behstates(k,numb_beh)=zeros;
totalAUC_behstates(k,numb_beh)=zeros;
numbtransients_behstates(k,numb_beh)=zeros;
for b=1:numb_beh
    AUC_bycell=sum(cell_AUC(states{b},:),1);
    %loops through each behavioral state in cell array to sum events per cell in that state
    numbtransients=sum(cell_events(states{b},:),1);
    totalAUC_behstates(k,b)=AUC_bycell';
    numbtransients_behstates(k,b)=numbtransients';
    %cells in rows, # of events in each behavior, each behavior in diff column
end

for b=1:numb_beh
    aveAUC_behstates(k,b)=totalAUC_behstates(:,b)./numbtransients_behstates(k,b);
    %calculates rate by dividing number of events by the length of that
    %behavioral state
end
end