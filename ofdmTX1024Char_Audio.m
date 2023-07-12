
clc;
clear all
close all

fs = 8e3;  % sample rate of ofdm signal 
% create chirp signal for sync 
N = 8000;  %number of samples for chirp sync pulse
t = [0:(N-1)]/fs;
f0 = 300;
f1 = 3000;
T = N/fs;
c = (f1 - f0)/T;
x = (c/2)*t.^2 + f0*t;
pream = sin(2*pi*x);
pream = pream / max(pream);






% samples from sync pulse to ofdm data
guardN = 512;

% time vector 
td = 5;
tdsamples = round(td * fs);
tdvec = zeros(1,tdsamples);


Nchar = 1024;  % make power of 2,  number of char for text
bits = Nchar * 8;

Gbins = 512;  % guard fft bins 
G2 = round(Gbins/2);  % guard fft bins on each side of spectrum
fftsize = 1024;
bitsperfft = fftsize - Gbins;

Nfftblocks = ceil ( bits / bitsperfft);



% get message from txt file
fid = fopen('C:\SDR_Work\textfile.txt');
mtxt = fread(fid);
fclose(fid);
mtxt = mtxt';
if length(mtxt) >= Nchar
    mtxt = mtxt(1:Nchar);
    disp('text file greater char limit')
else
    z = length(mtxt);
    z = Nchar - z;
    mtxt = [mtxt zeros(1,z)];
    disp('textfile under char limit')
end
mbits = dec2bin(mtxt,8);
mbits = reshape(mbits',1,[]);
mbits = mbits(1:bits);
txbits = zeros(1,bits);
for k = 1:bits
    if mbits(k) == '1'
        txbits(k) = 1;
    else
        txbits(k) = -1;
    end
end

txbits = reshape(txbits,bitsperfft,[]);
txbits = txbits';
[nblocks,nbits] = size(txbits);




% make pilot vector
xp = randn(1,fftsize);
xp = xp / max(xp);
xp = xp * 13;
xp = 2 * pi * xp;


% pilot vector with  guard on each side
xp = exp(1i*xp);
xp(1:G2) = 0;
xp(end-(G2-1):end) = 0;

% make data matrix bpsk,  1 = 1,  0 = -1

xd = [zeros(nblocks,G2)   txbits    zeros(nblocks,G2) ];  % guard on each side




for k = 1:nblocks

    if k == 1
    %relate data vector to pilot vector
        xd(1,:) = xp .* xd(1,:);  % e^a * e^b =  e^(a + b)
    else
        xd(k,:) = xd(k-1,:) .* xd(k,:);
    end
end


    



% make data be real: conj and symetric 
xp1 = conj( xp(end:-1:2)  );
xp = [xp1 xp];

xd1 = conj( xd(:,end:-1:2) );
xd = [xd1 xd];
[nblock,fftsamples] = size(xd);

% pilot time vector
xp = ifftshift(xp);
xp = ifft(xp);
xp = xp / max(abs(xp));

% data time matrix
xd = ifftshift(xd,2);
xd = ifft(xd,[],2);

xdtx = zeros(nblock,2*fftsamples );

% repeat data blocks twice for cyclic prefix
for k = 1:nblock

    xdtx(k,:) = [ xd(k,:) xd(k,:) ];
    
end

xdtx = reshape(xdtx.',1,[]);   %  1xN
xdtx = xdtx / max(abs(xdtx));



% put the whole movie together
xz = zeros(1,guardN);
xt = [tdvec pream xz xp xp xdtx tdvec];

% ensure bandwidth, reduce out of band signals
hlpf = fir1(64,0.85);
xt = filter(hlpf,1,xt);
xt = xt / max(abs(xt));  % 1xN

figure()
plot(real(xt))
hold on
plot(imag(xt))
hold off


datafile = [ real(xt).' ];



 audiowrite('c:\SDR_Work\ofdmtest1024Char8khzAudio.wav',datafile,fs,'BitsPerSample',16);
% 













