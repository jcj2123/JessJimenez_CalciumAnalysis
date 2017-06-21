
function[cellrates_behstates,states,beh_length]=calc_beh_rates_forboot(behavior, cell, type)
%outputs rates and events per cell per behavior state specified
%calls on the 'define_beh_states' function to define behavior columns for
%diff imaging sessions
%inputs behavior data, cell data, and specify what imaging session as
%'type'
%'type': (but see more behavior column detail in 'define_beh_states'
%'EPM'=elevated plus
%'EZM'= elevated zero maze
%'OFT'= open field test
%'NOVEL'= novel object task
[states,numb_beh,beh_length]= define_beh_states(behavior,type);
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

k=1:size(cell,2);
cellrates_behstates(k,4)=zeros;
cellevents_behstates(k,4)=zeros;
for b=1:2
    events_bycell=sum(cell(states{b},:),1);
    %loops through each behavioral state in cell array to sum events per cell in that state
    cellevents_behstates(k,b)=events_bycell';
    %cells in rows, # of events in each behavior, each behavior in diff column
end

for b=1:numb_beh
    cellrates_behstates(k,b)=cellevents_behstates(:,b)./beh_length(b);
    %calculates rate by dividing number of events by the length of that
    %behavioral state
end
cellrates_behstates(k,3)=cellrates_behstates(k,2)-cellrates_behstates(k,1);
%rate differences column2-column1, to be used in bootstrap script later
%cond 2 selective
cellrates_behstates(k,4)=cellrates_behstates(k,1)-cellrates_behstates(k,2);
%rate differences column1-column2, to be used in bootstrap script later
%cond1selective
end