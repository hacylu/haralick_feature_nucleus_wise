function [featvals names]= Lharalick_img_nuclei_wise(img,mask,Num_level)
% 
% img: org tissue image.
% mask: nuclei mask
% Num_level:8
% 
%  feats = haralick_img_nuclei_wise(img,mask,Num_level)
% cheng adapted haojia's function here

%function feats = haralick(img,fg_bg_mask,graylevels,hws,dist,fg_bg)
%fg_bg_mask is a mask of the foreground as ones (fg_bg =1) or background as ones (fg_bg =0)
%Returns 13 Haralick features in a structure containing names and features
%Uses MEX file graycomtx_3d for co-occurence matrix calculations
%James Monaco
% FEATS is a struct: img3 are the features, names are the names
% Options: 
% graylevels : 64,128,256
%hws :1,2,3
%dist: 1,2
% IMG needs to be uint16.

mask_size=size(mask);
img_size=size(img);
if mask_size~=img_size
    error('the size of image and mask must be same')
end


% ttime=cputime;
feats.names = {'contrast_energy','contrast_inverse_moment', ...
    'contrast_ave','contrast_var','contrast_entropy', ... 
    'intensity_ave','intensity_variance','intensity_entropy', ...
    'entropy,','energy', ...
    'correlation', ...
    'information_measure1','information_measure2'};

cc=bwconncomp(mask);
mask_label = bwlabel(mask);

% z coords of window, removing z coords outside of image
for nuclei_ind=1:cc.NumObjects
    roi_ind=cc.PixelIdxList{nuclei_ind};
    [ROI_row,ROI_colm]=ind2sub(mask_size,roi_ind);
    Gray_lim(1)=min(img(roi_ind));
    Gray_lim(2)=max(img(roi_ind));
    org_nuclei=img(min(ROI_row):max(ROI_row),min(ROI_colm):max(ROI_colm));
    mask_label_nuclei=mask_label(min(ROI_row):max(ROI_row),min(ROI_colm):max(ROI_colm));
    roi_mask=logical(mask_label_nuclei==nuclei_ind);
   
   
    
 %%   transpose roi_ind_nuclei to coordinate with the tansformeation within function graycomatrix_myversion
    
 
    SGLD=graycomatrix_myversion(org_nuclei,roi_mask,'offset',[0 1],'GrayLimits',Gray_lim,'NumLevels',Num_level, 'Symmetric', true);
	            SGLD = double(SGLD);
                SGLD=SGLD/sum(sum(SGLD));
	            %%% Calculate Statistics %%%%
	            [pi,pj,p] = find(SGLD);  
                
	            if length(p) <= 1
	                continue;
	            end
	            p = p/sum(p); pi = pi-1; pj = pj-1;

	            %marginal of x
	            px_all = sum(SGLD,2);                                       
	            [pxi,junk,px] = find(px_all>(eps^2));
	            px = px/sum(px); pxi = pxi-1;
% 	            if length(px) <=3
% 	                continue;
% 	            end             

	            %marginal of y
	            py_all = sum(SGLD,1)';                                       
	            [pyi,junk,py] = find(py_all>(eps^2)); 
	            py = py/sum(py); pyi = pyi-1;


	            %%% Calculate Contrast Features %%%%
	            all_contrast = abs(pi-pj);
	            [sorted_contrast,sind] = sort(all_contrast);
	            ind = [find(diff(sorted_contrast)); length(all_contrast)];
	            contrast = sorted_contrast(ind);
	            pcontrast = cumsum(p(sind));
	            pcontrast = diff([0; pcontrast(ind)]);

	            contrast_energy = sum( contrast.^2 .* pcontrast );                  
	            contrast_inverse_moment = sum( (1./(1+contrast.^2)) .* pcontrast );      
	            contrast_ave = sum( contrast .* pcontrast ); 
	            contrast_var = sum( (contrast - contrast_ave).^2 .* pcontrast );
	            contrast_entropy = -sum( pcontrast.*log(pcontrast) );

	            %%% Calculate Intensity Features %%%%
	            all_intensity = (pi+pj)/2;
	            [sorted_intensity,sind] = sort(all_intensity);
	            ind = [find(diff(sorted_intensity)); length(all_intensity)];
	            intensity = sorted_intensity(ind);
	            pintensity = cumsum(p(sind));
	            pintensity = diff([0; pintensity(ind)]);

	            intensity_ave = sum( intensity .* pintensity );                     
	            intensity_variance = sum( (intensity-intensity_ave).^2 .* pintensity ); 
	            intensity_entropy = -sum( pintensity.*log(pintensity) ); 

	            %%% Calculate Probability Features %%%%
	            entropy = -sum( p.*log(p) );                     
	            energy = sum( p.*p );                     

	            %%% Calculate Correlation Features %%%%
	            mu_x = sum(pxi.*px);                                 
	            sigma_x = sqrt(sum( (pxi-mu_x).^2 .* px ));            
	            mu_y = sum(pyi.*py);                                  
	            sigma_y = sqrt(sum( (pyi-mu_y).^2 .* py));

	            if sigma_x==0 || sigma_y==0
	                warning('Zero standard deviation.');
	            else
	                correlation = sum( (pi-mu_x).*(pj-mu_y).* p ) / (sigma_x*sigma_y); 
	            end

	            %%% Calculate Information Features %%%%    
	            [px_grid,py_grid] = meshgrid(px,py);
	            [log_px_grid,log_py_grid] = meshgrid(log(px),log(py)); 
	            h1 = -sum( p .* log(px_all(pj+1).*py_all(pi+1)) );
	            h2 = -sum( px_grid(:).*py_grid(:).*(log_px_grid(:)+log_py_grid(:)) );  
	            hx = -sum( px.*log(px) );
	            hy = -sum( py.*log(py) );

	            information_measure1 = (entropy-h1)/max(hx,hy);
	            information_measure2 = sqrt(1-exp(-2*(h2-entropy)));

	            for k = 1:length(feats.names)
	                feats.img3(k,nuclei_ind) = eval(feats.names{k});
	            end

end
 
            

% 
% 
%         if mod(nuclei_ind,5)==0, % print a dot for every row. change number to skip rows
%             dt=toc;
%             if dt<0.5,
%                 fprintf('-');
%             else
%                 fprintf('.');pause(0.04);    % to facilitate a break
%             end
%             tic
%         end
%         milestones=round(linspace(1,cc.NumObjects,11));
%         milestone_percents=0:10:100;
%         if any(milestones==nuclei_ind),
%             fprintf('%d%%',milestone_percents(milestones==nuclei_ind));pause(0.04);
%         end
%                         
%             
%                 
% if any(isinf(feats.img3(:))) || any(isnan(feats.img3(:))) || any(~isreal(feats.img3(:)))
%     ninfs=sum(isinf(feats.img3(:)));
%     nnans=sum(isnan(feats.img3(:)));
%     fprintf('There are some Infs (%d) or NaNs (%d) in the feature image.\n',ninfs,nnans)
% end

feats.img3 = shiftdim(feats.img3,1); 
% fprintf('  Time: %d\n',cputime-ttime);
%% take statistics here
%% statistics: mean, median, standard deviation, range, kurtosis, skewness , across bounds for each haralick feature
% check if the current nuclear morpholgy has no haralick features
all_feat=feats.img3;
if ~isempty(all_feat) && size(all_feat,1)>3
    
    
    featvals=[];
    featvals=[mean(all_feat) median(all_feat) std(all_feat)  range(all_feat) kurtosis(all_feat) skewness(all_feat)];
else
    featvals=[ zeros(1,78)];
end
%% feature names organization
count = 1;
modifier = [{'mean'} {'median'} {'std'} {'range'} {'kurtosis'} {'skewness'} ];
%the format is like this: morphofeaturenames_haralickfeatuename_statistcname
feature_names=feats.names;

for mi = 1:numel(modifier)
    for j = 1:numel(feature_names)
        names{count} = ['HarralickNuWise-' modifier{mi} '(' feature_names{j} ')'  ];
        count = count + 1;
    end
end


