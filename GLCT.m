function [Spec,tfr,tt,ff] = GLCT(x,N,fs,hlength);
%      General Linear Chirplet Transform.
%	x      : Signal.
%	N      : Number of calculated LCT with different chirp rate.
%	fs     : Sample Frequency .
%	hlength: Length of window function.

%	tfr    : Time-Frequency Representation.
%	tt     : Time.
%      ff     : Frequency.

[xrow,xcol] = size(x);
if (nargin < 3),
    error('At least 3 parameter is required');
end;
Siglength=xrow;

if (nargin < 4),
    hlength=floor(Siglength/4);
end;

hlength=hlength+1-rem(hlength,2);
h = tftb_window(hlength);

%t=1:xrow;
%[trow,tcol] = size(t);

[hrow,hcol]=size(h); Lh=(hrow-1)/2;

h=h/norm(h);

Frelength=round(Siglength/2);

slope=(-pi/2+pi/(N+1)):pi/(N+1):(pi/2-pi/(N+1));
Spec=zeros(N,Frelength,Siglength);%%%

t=Siglength/fs;

for i=1:N
    [Spec(i,:,:)] = LCT(x,tan(slope(i))*(fs/2)/t,fs,h);
end

GLCTspec=zeros(Frelength,Siglength);

for i=1:Frelength
    for j=1:Siglength
        absSpec=abs(Spec(:,i,j));
        index=find(max(absSpec)==absSpec);
        GLCTspec(i,j)=Spec(index(1),i,j);
    end
end

tfr=GLCTspec;

if (nargout==3),
    tt=(1:Siglength)/fs;
elseif (nargout==4),
    tt=(1:Siglength)/fs;
    ff=0:fs/Siglength:(fs/2)-(fs/2)/Siglength;
end;
