

% ofdm decoder
clear all
close all
clc
plotcrap = 1;  % set to one to see debug plots

FILTERSPEC = '.wav';
TITLE = 'Pick a Recorded OFDM IQ wav file';
FILE = 'C:\sdr#\sdrsharp-x86-noskin(1)\';

[FILENAME, PATHNAME, FILTERINDEX] = uigetfile(FILTERSPEC, TITLE, FILE);
xstring = [PATHNAME FILENAME];

[x,fswav] = audioread(xstring);

% tx bandwidth

fs = 8e3;  % sample rate of ofdm signal 
% pulse at start for sync 
N = 8000;  %number of samples for chirp sync pulse
t = [0:(N-1)]/fs;
f0 = 300;
f1 = 3000;
T = N / fs;
c = (f1 - f0)/T;
x1 = (c/2)*t.^2 + f0*t;
pream = sin(2*pi*x1);
Nguard = 512;   % time between pulse and preamble fft data

Nchar = 1024;  % make power of 2
bits = Nchar * 8;  % total number of bits

Gbins = 512;  % guard fft bins
G2 = round(Gbins/2);  % guard fft bins on each side
Nfft2 = 1024;
bitsperfft = Nfft2 - Gbins;
Nfft = (2*Nfft2)  - 1 ; 
Nblocks = ceil(bits / bitsperfft );

bitsmatrix = zeros(Nblocks,Nfft);



%                    xp            xd           guard
sigtime = Nguard + 2*Nfft + 2*Nfft*Nblocks + 2*Nguard;






% just get one channel of audio
x = x(:,1);
x = x.';






if fswav ~= fs
    n = gcd(fs,fswav);
    p = fs / n;
    q = fswav / n;
    x = resample(x,p,q);   % resample file to expected sample rate
end


hpw = conj(pream(end:-1:1));
xdet = filter(hpw,1,x);

if plotcrap
    figure(100)
    plot(abs(xdet))
    title('chirp detection')
end

[maxv, maxi] = max(abs(xdet)); % maxi will be at end of chirp
xstart = maxi + Nguard  ;
xend = xstart + sigtime ;
x = x(xstart:xend);

% if plotcrap
%     figure(101)
%     xplot = abs(x) / max(abs(x));
%     xplot = 20 * log10(xplot);
%     yplot = length(xdet);
%     yplot = zeros(1,yplot);
%     yplot(maxi) = -40;
%     plot(xplot)
%     hold on
%     plot(yplot)
%     hold off
% end
z = round(0.95*Nfft);


a = z + 1 ;
b = a + Nfft   - 1 ;
c = (2*Nfft) - a;
c = round(c);

a1 =   2*Nfft + z + 1 ;
b2 = a1 + Nfft  - 1 ;


for k = 1:Nblocks
    xp = x(a:b);
    xd = x(a1:b2);
    Xp = fft(xp);
    Xd = fft(xd);
    if plotcrap 
        length(Xp)
        length(Xd)
    end
    Xp = Xp( 1:floor(end/2) );  % 
    Xd = Xd( 1:floor(end/2) );

    if plotcrap 
        figure(456)
        plot(abs(Xp))
        hold on
        plot(abs(Xd))
        hold off
        title('fft of pilot and data')
        legend('pilot fft', 'data fft')
        pause(1)
    end
    det = angle(Xp ./ Xd);          % ratio to get delta phase diff
    det = abs(det);
    if plotcrap
        figure(123)
        plot(det)
        title('phase vector')
        pause(1)
    end
    n = length(det);
    bitsmatrix(k,1:n) = det;
    a = a + c + z + 1;
    b = a + Nfft - 1;
    a1 = a1 + c + z + 1;
    b2 = a1 + Nfft - 1;
end





bitsmatrix = bitsmatrix(:,G2+1:end);
bitsmatrix = bitsmatrix(:,1:bitsperfft );
bits = reshape(bitsmatrix',1,[]);

thres = pi/2;

det = bits;






det(det <= thres) = 1;
det(det > thres) = 0;

% if plotcrap
%     figure(105)
%     plot(det)
%     title('decoded bits vector')
% end


det = reshape(det,8,[]);
det = det';


charmes = zeros(1,Nchar);
w = [7:-1:0];
w = 2.^w;


for k = 1:Nchar
    x = det(k,:);
    x = x .* w;
    x = sum(x);
    charmes(k) = x;
end



%     clc

   xstr =  char(charmes)

   msgbox(xstr,'replace')
















