function S = getImageEdges(I,useGPU)

if nargin < 2
    useGPU = false;
end

if useGPU
    S = getImageEdges_GPU(I);
else
    
    I = double(I);

    % 2nd derivative
    % S2 = abs(imfilter(I,fspecial('log',14,2),'symmetric'));

    % 1st derivative
    [Gx,Gy] = smoothGradient(I,1);
    S1 = sqrt(Gx.^2 + Gy.^2);

    % Normalize and sum
    % S = S1/prctile_fast(S1(S1>0),80) + S2/prctile_fast(S2(S2>0),80);
    S = S1;
    S = S - min(S(:));
    S = uint8(255*S/max(S(:)));

    % Smooth out with a closing.
    S = imclose(S,strel('disk',3));

end

end


function S = getImageEdges_GPU(I)

try
    I = single(I);
    I_d = gpuArray(I);

    % 2nd derivative
%     S2_d = abs(imfilter(I_d,fspecial('log',14,2),'symmetric'));
%     S2 = gather(S2_d);
%     clear S2_d
%     
%     S2 = S2/prctile_fast(S2(S2>0),80);
    
    % 1st derivative
    [Gx_d,Gy_d] = smoothGradient(I_d,1);
    S1 = gather(sqrt(Gx_d.^2 + Gy_d.^2));
    clear I_d Gx_d Gy_d
    
%     S1 = S1/prctile_fast(S1(S1>0),80);
    

    % Normalize and sum
    S = S1;% + S2;

    S = S - min(S(:));
    S8 = uint8(255*S/max(S(:)));
    
    
    % Smooth out with a closing.
    S8_d = gpuArray(S8);
    S_d = imclose(S8_d,strel('disk',3));
    clear S8_d

    S = gather(S_d);
    clear S_d

catch ME
    reset(gpuDevice())
    gpuDevice([])
    rethrow(ME)
end
end

% function v = prctile_fast(x,p)
% % PERCTILE_FAST  An approximate percentile function with no error checking
% %
% % p - percentile [0,100]
% 
% % James Kapaldo
% % 2016-10-28
% x = double(x(:));
% x = sort(x,1);
% i = (1:numel(x))';
% 
% idx = round((p/100)*numel(x));
% 
% v = nakeinterp1(i,x,idx);
% 
% end