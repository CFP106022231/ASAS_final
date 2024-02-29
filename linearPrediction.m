
function [excitat, Aall] = linearPrediction(x, framelen, p, fs)
    [y,fs1] = audioread(x);
    y = y(:,1);
    y = resample(y,fs,fs1);
    Aall = [];
    %%
    L = framelen*fs;

    sw.emphasis = 1; % default = 1
    sw.visOnly = 0; % Visualize the signal only. default = 0.
    sw.datavisual = 0;
    numFrames = floor(length(y)/L);
    
    excitat = zeros(size(y));
    e_n = zeros(p+L,1);

    if sw.emphasis == 1
        y_emph = filter([1 -0.95],1,y); 
                    %[PARAM] -0.95 may be tuned anywhere from 0.9 to 0.99
    else
        y_emph = y;
    end
    
    %% Linear prediction and smooth synthesis of the estimated source e_n
    win = ones(L,1); % Rectangular window.
    
    % get excitation excitat
    for kk = 1:numFrames % frame index
        
        ind = (kk-1)*L+1:kk*L; % pointing to a frame of length L
        ywin = y_emph(ind).*win;
        %------------------- lpc term:A ,Aall:存所有frame的A -------------------
        A = lpc(ywin,p);
        Aall = [Aall; A]; % vertical concat
        %-------------------
        if kk == 1
            e_n(p+1:end) = filter(A,[1],ywin);
        else
            ywin_extended = y_emph((kk-1)*L+1-p:kk*L);
            e_n = filter(A,[1],ywin_extended);
        end
        excitat(ind) = e_n(p+1:end);
%         excitat2 = filter(1,[1 -0.95],excitat);

    end


end