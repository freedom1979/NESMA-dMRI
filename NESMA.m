% nonlocal estimation of multispectral magnitudes (NESMA) 
% referred to the paper "A simple and fast adaptive nonlocal multispectral filtering algorithm for
% efficient noise reduction in magnetic resonance imaging"


Subjects = {'sub-CC00069XX12_ses-26300','sub-CC00099AN18_ses-34200','sub-CC00117XX10_ses-38200',...
            'sub-CC00122XX07_ses-42000','sub-CC00126XX11_ses-43100','sub-CC00138XX15_ses-46200',...
            'sub-CC00162XX06_ses-53600','sub-CC00164XX08_ses-54000','sub-CC00168XX12_ses-55700',...
            'sub-CC00170XX06_ses-56100'};
        
Num_of_Subs = length(Subjects);
for iii = 1:Num_of_Subs
    tic
    
    dwi_nii = load_untouch_nii([Subjects{iii} '_dwi.nii']);
    brain_mask_nii = load_untouch_nii([Subjects{iii} '_brainmask.nii']);
    
    dwi_data = dwi_nii.img;
    brain_mask = brain_mask_nii.img;
    [M,N,L,S_spectra] = size(dwi_data);  data_size = [M,N,L];
    
    sigma = zeros(S_spectra,1);
    for i_Sp = 1:S_spectra
        T = squeeze(dwi_data(:,:,:,i_Sp));
        Tt = T(:);
        pd = fitdist(Tt, 'Rayleigh');
        sigma(i_Sp,1) = std(pd);
    end
    
    % sigma = 0.56;
    % R = 50;
    w = [4,4,4];
    seau_t = 0.04;
    dwi_data_filtered = zeros(M,N,L,S_spectra);
    
    
    for i_L = 1:L
        for i_N = 1:N
            for i_M = 1:M
                if brain_mask(i_M,i_N,i_L) ~= 0
                    v = [i_M,i_N,i_L];
                    [xc,yc,zc] = candidate_voxel(v,w,data_size);
                    [m,~] = size(xc);
                    
                    d = zeros(m,1);
                    for i_m = 1:m
                        dc = zeros(S_spectra,1);
                        dc2 = zeros(S_spectra,1);
                        for i_spectra = 1:S_spectra
                            si = dwi_data(i_M,i_N,i_L,i_spectra);
                            sj = dwi_data(xc(i_m),yc(i_m),zc(i_m),i_spectra);
                            dc(i_spectra,1) = (si-sj)^2;
                            dc2(i_spectra,1) = si^2;
                        end
                        d(i_m,1) = sum(dc)/sum(dc2);
                    end
                    [Dd,Ii] = sort(d); % ascending
                    
                    R_Dd = length(Dd);
                    R = floor(R_Dd*seau_t);
                    
                    for i_spectra = 1:S_spectra
                        SJ = zeros(R,1);
                        for i_R = 1:R
                            sj = dwi_data(xc(Ii(i_R)),yc(Ii(i_R)),zc(Ii(i_R)),i_spectra);
                            SJ(i_R,1) = sj^2;
                        end
                        T = sum(SJ)/R-2*sigma(i_spectra)^2;
                        dwi_data_filtered(i_M,i_N,i_L,i_spectra) = sqrt(max(T,0));
                    end
                    
                    
                end
            end
        end
    end
    
    dwi_nii.img = dwi_data_filtered;
    save_untouch_nii(dwi_nii, [Subjects{iii} 'dwi_filtered.nii']);
    
    toc
end

