function [fus_DS] = GLP_MTF_HPM_DS(MS,PAN,u,sensor,ratio)
%GLP_MTF_HPM_DS  GLP-MTF-HPM-DS
%   Inputs:
%       MS:      MS image upsampled at PAN scale;
%       PAN:     PAN image;
%       u:       Parameters selected for different scenarios.
%       sensor:  String for type of sensor (e.g. 'WV2','IKONOS');
%       ratio:   Scale ratio between MS and PAN. Pre-condition: Integer value.
%
%   Outputs:
%       fus_DS:  Pansharpened image.
%
%   Reference:
%       [1]P. Wang, H. Yao, C. Li, G. Zhang, and H. Leung, ¡°Multiresolution analysis 
%       based on dual-scale regression for pansharpening,¡± IEEE Trans. Geosci. 
%       Remote Sens., vol. 60, pp. 1¨C19, 2022.
%       [2]G. Vivone, R. Restaino, and J. Chanussot, ¡°A regression-based highpass 
%       modulation pansharpening approach,¡± IEEE Trans. Geosci. Remote Sens., 
%       vol. 56, no. 2, pp. 984¨C996, Feb. 2018.
%       [3]G. Vivone, R. Restaino, and J. Chanussot, ¡°Full scale regression-based
%       injection coefficients for panchromatic sharpening,¡± IEEE Trans. Image
%       Process., vol. 27, no. 7, pp. 3418¨C3431, 2018.

h = genMTF(ratio, sensor, size(MS,3));

P_L = zeros(size(MS));
for i = 1:1:size(MS,3)
    pl = imfilter(PAN,real(h(:,:,i)),'replicate'); 
    t = imresize(pl,1/ratio,'nearest');
    P_L(:,:,i) = interp23tap(t,ratio);
end

fus_DS = zeros(size(MS));
for j = 1 : size(MS,3)
    A = MS(:,:,j);
    P = PAN;
    PL = P_L(:,:,j);
    
    %%%% Regression coefficients
	C_A_P = cov(A(:),P(:));
    C_A_PL = cov(A(:),PL(:));
    C_P_PL = cov(P(:),PL(:));
    g = u * (C_A_P(1,2) / C_P_PL(1,2)) + (1 - u) * (C_A_PL(1,2) / C_P_PL(1,2));
    
    cb = mean(A(:))/g - mean(P(:));
        
    %%% Fusion rule
    fus_DS(:,:,j) = A .* (P + cb) ./ (PL + cb + eps);
end

end

