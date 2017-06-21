function [skaggs_info]=calc_skaggs_info(cell_data,qualifiedcells,occupancy_map_raw,placemaps_raw)
    cell=cell_data(:,qualifiedcells);
    placemaps=placemaps_raw;
    
    ave_rate=(sum(cell,1)./size(cell,1))';
    ratebins_byave(size(placemaps,1),size(placemaps,2),size(cell,2))=zeros;
    skaggs_info_bins(size(placemaps,1),size(placemaps,2),size(cell,2))=zeros;
    skaggs_info1(1,size(placemaps,2),size(cell,2))=zeros;
     skaggs_info(size(cell,2),1)=zeros;
    for k=1:size(cell,2);
    
   %lambdai/lambda(firing rate per bin divided by ave rate of cell)
   ratebins_byave(:,:,k)=placemaps(:,:,k)./ave_rate(k);

   %(lambdai/lambda)*log2(lambdai/lambda)
   ratebins_byave_log2=ratebins_byave.*(log2(ratebins_byave));
   
   prob_bins=occupancy_map_raw./size(cell,2);
   
   skaggs_info_bins(:,:,k)=prob_bins.*ratebins_byave_log2(:,:,k);
   skaggs_info1(1,:,k)=sum(skaggs_info_bins(:,:,k));
   skaggs_info(k,1)=sum(skaggs_info1(:,:,k),2);
   
    end
   
   
   
   