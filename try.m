%% lpc parameter
clear; close all;
p = 20;
fs = 16000;
framelen = 0.032;

%% 音樂與人聲input

INST_DIR = './input/';
VOICE_DIR = './input/';
outdir = './output/';
instr = 'flu.wav';
voice = 'break.m4a';
INST_FILE = instr;
VOICE_FILE = voice;
sw.plot = 1; % visualization

[E1, A1] = linearPrediction([INST_DIR INST_FILE], framelen, p, fs);  %鋼琴
[E2, A2] = linearPrediction([VOICE_DIR VOICE_FILE], framelen, p, fs);  % 人聲

% [Eme, Ame] = linearPrediction('2minute.m4a', 0.032, 20);
% [Ehow, Ahow] = linearPrediction('how.mp4', 0.032, 20);
% [Emor, Amor] = linearPrediction('motorcycle.wav', 0.032, 20);

% noise = wgn(58142,1,-50); % 白噪音：測試filter的聲音，結果會像氣音

input = 0.5*E1;
[voice, fs_v]= audioread([VOICE_DIR VOICE_FILE]); % import for plotting

if fs_v ~= fs
    voice = resample(voice, fs, fs_v);
end

%% parameter

L = framelen*fs;
numFrames = floor(length(input)/L);
win = ones(L,1);
cross = zeros(size(input));
e_n = zeros(p+L,1);
Nfreqs = 2^nextpow2(2*L-1)/2; % Num points for plotting the inverse filter response
df = fs/2/Nfreqs;
ff = 0:df:fs/2-df;

videoFWriter  = vision.VideoFileWriter('./output/flubreakfast.avi','FileFormat','AVI','AudioInputPort',true,'FrameRate', 1/framelen);

for kk=1:numFrames
    ind = (kk-1)*L+1:kk*L;
    pwin = input(ind).*win;
    vwin = voice(ind).*win;
    v_emph = filter([1 -0.95],1,vwin);
    if kk == 1
        e_n(p+1:end) = filter(1,A2(kk,:),pwin);
    else
        pwinext = input((kk-1)*L+1-p:kk*L); % 往前偷看p個點
        e_n = filter(1, A2(kk,:), pwinext);

    end
    
    cross(ind) = e_n(p+1:end);
    
    % visualisation
    if sw.plot % && (mod(kk, 2)==0)

%         h = figure;
        subplot(3, 1, 1);
        % 畫 input, A1
        Ysc= fft(pwin, 2*Nfreqs); % spectrum of source e_n
        Asc = lpc(pwin, p); % spectrum contour of e_n
        %[H,W] = freqz([0 -A(2:end)], [1], Nfreqs);
        [Hsc,Wsc] = freqz([1],Asc, Nfreqs); % freq resp of filter
        
        Hmagsc = 20*log10(abs(Hsc));
        Ymagsc = 20*log10(abs(Ysc(1:Nfreqs)));
        Hmaxsc = max(Hmagsc);
        offsetsc = max(Hmagsc) - max(Ymagsc);
        plot(ff,Hmagsc); hold on;
        plot(ff,Ymagsc+offsetsc,'r'); hold off;
        set(gca,'xlim',[0 fs/2],'ylim',[Hmaxsc-50, Hmaxsc+5]);
        xlabel('Hz')
        title("source e_n");

        subplot(3, 1, 2);
        % 畫 voice, A2
        Yflt= fft(v_emph, 2*Nfreqs); % spectrum of filter
        %[H,W] = freqz([0 -A(2:end)], [1], Nfreqs);
        [Hflt,Wflt] = freqz([1],A2(kk, :), Nfreqs); % freq resp of filter
        
        Hmagflt= 20*log10(abs(Hflt));
        Ymagflt= 20*log10(abs(Yflt(1:Nfreqs)));
        Hmaxflt= max(Hmagflt);
        offsetflt= max(Hmagflt) - max(Ymagflt);
        plot(ff,Hmagflt); hold on;
        plot(ff,Ymagflt+offsetflt,'r'); hold off;
        set(gca,'xlim',[0 fs/2],'ylim',[Hmaxflt-50, Hmaxflt+5]);
        xlabel('Hz')
        title("'vocal' filter");
        
        subplot(3, 1, 3);
        % 畫 cross, A2
        Y = fft(cross(ind),2*Nfreqs); % spectrum of cross-synth sound
        %[H,W] = freqz([0 -A(2:end)], [1], Nfreqs);
        [H,W] = freqz([1],A2(kk, :), Nfreqs); % preq resp of filter
        
        Hmag = 20*log10(abs(H));
        Ymag = 20*log10(abs(Y(1:Nfreqs))); 
        Hmax = max(Hmag);
        offset = max(Hmag) - max(Ymag);
        plot(ff,Hmag); hold on;
        plot(ff,Ymag+offset,'r'); hold off;
        set(gca,'xlim',[0 fs/2],'ylim',[Hmax-50, Hmax+5]);
        xlabel('Hz');
        title("cross-synthesised");
        
%         F = getframe(h);
%         step(videoFWriter, F.cdata, cross(ind)); 
%         drawnow;
          pause;
    end
    
end
% cross_de = filter(1,[1 -0.95],cross);


release(videoFWriter);
