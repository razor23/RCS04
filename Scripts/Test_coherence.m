%% Coherence check %%

clear all
close all

x=[0:0.2:5000];
y=sin(x);

figure; plot(x,y);

noise=randn(1,length(x));
y1=y+noise;

noise=randn(1,length(x));
y2=y+noise;

figure; 
subplot(2,1,1);
plot(x,y1)
subplot(2,1,2);
plot(x,y2)

%%
fs=250;
[cxy,f]=mscohere(y1,y2,fs,fs/2,[1:100],fs);
figure; 
plot(f,cxy,'b'); hold on;

%% Cut out chunks of data %%

for g=1:100
    % Pick random number %
    rnd=round(abs(rand(1))*(length(y1)-51));
    
    y1(rnd:rnd+50)=[];
    y2(rnd:rnd+50)=[];
end

fs=250;
[cxy,f]=mscohere(y1,y2,fs,fs/2,[1:100],fs);
plot(f,cxy,'r');
legend('Normal data','Chunks cut');

%%

yabs=abs(y1);

figure; plot(yabs)

sabs3=std(yabs)*3;

%% Artifact rejection - longer indexing on either side

ind=yabs>sabs3;

figure; 
sb1=subplot(2,1,1); plot(yabs); 
sb2=subplot(2,1,2); plot(ind);
linkaxes([sb1,sb2],'x')

indsm=smooth(ind,125);
figure; plot(indsm)

ind2=indsm>0.1;

figure; plot(ind2);


