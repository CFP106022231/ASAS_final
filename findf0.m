function [linAmp, f] = findf0(x, framelen, p, fs)
    [y,fs1] = audioread(x);
    y = y(:,1);
    fs = 16000;
    y = resample(y,fs,fs1);
    y_emph = filter([1 -0.95],1,y); 
%     y_emph = y;
    framelen = 0.032; % second. [INVESTIGATE]
    p = 20; % linear prediction order. [INVESTIGATE]
    
    L = framelen*fs;
    numFrames = floor(length(y)/L);
    Nfreqs = 2^nextpow2(2*L-1)/2; % Num points for plotting the inverse filter response
    df = fs/2/Nfreqs;
    ff = 0:df:fs/2-df;
    win = hann(L);
    f = [];
    linAmp = [];
    miss = zeros(length(y),1);
    for kk = 1:numFrames % frame index
        ind = (kk-1)*L+1:kk*L; % pointing to a frame of length L
        ywin = y_emph(ind).*win;
        Y = fft(ywin,2*Nfreqs);
        Ymag = 20*log10(abs(Y(1:Nfreqs)));
    
        [mag, indP] = findpeaks(Ymag);%,'MinPeakProminence',31,'Annotate','extents');
    
%         yi = interp1((indP-1)*df, mag, ff, "linear");
        [peak, ind2] = findpeaks(mag,'MinPeakProminence',15,'Annotate','extents');
        f0ind = find(Ymag==peak(1),1,'first');
    
        linAmp = [linAmp 10^(peak(1)/20)];
%         miss(ind) = linAmp*sin(2*pi*(f0ind-1)*df*(1:512)/fs);
        f = [f (f0ind-1)*df];
    
    %     plot((f0ind-1)*df, peak(1),'x'); hold on;
    %     plot(ff, yi);hold on;
    %     plot(ff,Ymag);hold off;
    %     set(gca,'xlim',[0 fs/2],'ylim',[-30, 25]);
    %     xlabel('Hz');
    %     drawnow;
    %     pause;
    end