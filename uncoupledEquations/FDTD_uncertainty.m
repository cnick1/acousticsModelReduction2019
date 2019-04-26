%this is a version 3 fdtd code
%This has in buildings and attempts to create a buffer zone to absorb
%energy traveling into the far field
%clear

%griding out 2-D problem
numOb=[100 100];
ratSpace=[2 2];
num=[numOb(1)*ratSpace(1) numOb(2)*ratSpace(2)];%number of elements in x and y
num_ele=prod(num);
dxy=[0.5 0.5];%dimensions of grid
[pos,vec]=buildRectMesh(num,dxy);
subset=repmat([1:numOb]*ratSpace(1),numOb(2),1)+repmat(([0:numOb(2)-1]*ratSpace(1)*num(1)),numOb(1),1)';
subset=reshape(subset.',prod(numOb),1);

%setting up constants
ro=1.2;%density of air
T=20%deg Celcius
c=331.3*sqrt(1+(T/273.15));%speed of sound

%setting up loop to run through time
dt=0.0005;%timestep in seconds
numSkip=4;%number of time steps to skip before saving
numTimeSave=240;
max_t=numTimeSave*numSkip*dt;%length of simulation in time
numStep=numSkip*numTimeSave;
con=dt*ro*c*c;

%defining an absorptive boundary  
numRings=12;minfact=0.8;
fact=[minfact:(1-minfact)/numRings:1-(1e-9)];
fact=cos(acos(minfact):-(acos(minfact))/numRings:0.00001);
F=ones(num_ele,1);
for i=1:numRings
boundB=[i:num(1)-i]+num(1)*(i-1);
boundR=([i:num(2)-i]*num(1))+(1-i);
boundT=num_ele-([i:num(1)-i]+num(1)*(i-1)+1);  
boundL=([i:num(2)-i]*num(1))+i;
F([boundR,boundL,boundT,boundB])=fact(i);
end

%setting up initial conditions
vx=zeros(num_ele+1,2);
vy=zeros(num_ele+1,2);
p=zeros(num_ele+1,2);
pOb=zeros(prod(numOb),numTimeSave);

%defining input signal/pulse
pul=zeros(numStep,1);

% pul(1:30)=hanning(30);
% numAvg=numTimeSave;

% freF=50;
% pul=sin(2*pi*freF*dt*[1:numStep]);
% numAvg=round(2/mean(freF)/dt/numSkip);

% freF=[40 60];
% [B,A] = BUTTER(4,freF/(0.5/dt));
% pul(1)=1;
% pul=filter(B,A,pul);
% numAvg=numTimeSave;

%locating source
fx=60;
%fx=[9:numOb(1)-8];
fy=31;
posf=fx+(fy-1)*num(1);

%making walls/buildings
wallBin=[];
nx=[60 90];ny=[106 121];
wallBin=[repmat([nx(1):nx(2)],1,diff(ny)+1)+reshape(repmat([ny(1):ny(2)]*num(1),diff(nx)+1,1),1,(diff(nx)+1)*(diff(ny)+1))];
nx=[109 124];ny=[40 80];
wallBin=[wallBin,repmat([nx(1):nx(2)],1,diff(ny)+1)+reshape(repmat([ny(1):ny(2)]*num(1),diff(nx)+1,1),1,(diff(nx)+1)*(diff(ny)+1))];
nx=[115 134];ny=[106 121];
wallBin=[wallBin,repmat([nx(1):nx(2)],1,diff(ny)+1)+reshape(repmat([ny(1):ny(2)]*num(1),diff(nx)+1,1),1,(diff(nx)+1)*(diff(ny)+1))];
nx=[60 90];ny=[136 151];
wallBin=[wallBin,repmat([nx(1):nx(2)],1,diff(ny)+1)+reshape(repmat([ny(1):ny(2)]*num(1),diff(nx)+1,1),1,(diff(nx)+1)*(diff(ny)+1))];
nx=[115 134];ny=[136 151];
wallBin=[wallBin,repmat([nx(1):nx(2)],1,diff(ny)+1)+reshape(repmat([ny(1):ny(2)]*num(1),diff(nx)+1,1),1,(diff(nx)+1)*(diff(ny)+1))];
wallBinVx=union(wallBin,wallBin-1);
wallBinVy=union(wallBin,wallBin-num(1));

%starting the loop in time
for n=1:numTimeSave
    n
    for u=1:numSkip
    vx(1:num_ele,2)=F.*(vx(1:num_ele,1)-dt/ro/dxy(1)*(p(vec(:,1),1)-p(1:num_ele,1)));
    vy(1:num_ele,2)=F.*(vy(1:num_ele,1)-dt/ro/dxy(2)*(p(vec(:,3),1)-p(1:num_ele,1)));
    %buildings/walls 
    vx(wallBinVx,2)=0;
    vy(wallBinVy,2)=0;
    
    p(1:num_ele,2)=F.*(p(1:num_ele,1)- con*(vx(1:num_ele,2)-vx(vec(:,2),2))/dxy(1) - con*(vy(1:num_ele,2)-vy(vec(:,4),2))/dxy(2));
    p(wallBin,2)=0;
    
    %forcing
    %p(posf,2)=p(posf,2)+pul(u+(n-1)*numSkip); 
    vx(posf,2)=vx(posf,2)+pul(u+(n-1)*numSkip); 
    vx(vec(posf,2),2)=vx(vec(posf,2),2)-pul(u+(n-1)*numSkip);
     vy(posf,2)=vy(posf,2)+pul(u+(n-1)*numSkip); 
    vy(vec(posf,4),2)=vy(vec(posf,4),2)-pul(u+(n-1)*numSkip);
    
    %replacing old data
    p(:,1)=p(:,2);
    vx(:,1)=vx(:,2);
    vy(:,1)=vy(:,2);
    end
   pOb(:,n)= p(subset,1);
end
    
%plotting
% figure(1)
% for i=1:numTimeSave;
%     surf([1:numOb(1)]*dxy(1)*ratSpace(1),[1:numOb(2)]*dxy(2)*ratSpace(2),real(reshape(pOb(:,i),numOb(1),numOb(2)).'));
%     view(45,70);
%     axis([0 numOb(1)*dxy(1)*ratSpace(1) 0 numOb(2)*dxy(2)*ratSpace(2) -30 30])
%     CAXIS([-30 30])
%     xlabel('x')
%     ylabel('y')
%     drawnow
%     M(i) = getframe;
% end

figure(3)
posx=[dxy(1)*ratSpace(1):dxy(1)*ratSpace(1):numOb(1)*dxy(1)*ratSpace(1)];
posy=[dxy(2)*ratSpace(2):dxy(2)*ratSpace(2):numOb(2)*dxy(2)*ratSpace(2)];
contourf(posx,posy,10*log10(0.00001+reshape(mean(abs(pOb(:,numTimeSave-numAvg+1:end).^2),2),numOb(1),numOb(2))).',[-60:2:40])
%contourf(posx,posy,10*log10(0.00001+reshape(pDiff,numOb(1),numOb(2))).',[-30:6:20])
%contourf(posx,posy,10*log10(0.00001+reshape(pSum/49,numOb(1),numOb(2))).',[-30:6:60])
%contourf(posx,posy,10*log10(0.00001+reshape(sum(pOb20.^2,2),numOb(1),numOb(2))).',[-30:6:60])
axis([numRings*dxy(1) (num(1)-numRings)*dxy(1) numRings*dxy(2) (num(2)-numRings)*dxy(2)])
xlabel('x(m)')
ylabel('y(m)')
%title('Level difference due to 20^0C change or 3.6% increase in c - 50Hz - 6dB contours')
%title('Level 20^0C - 50Hz - 6dB contours')
title('Level averaged 0^0C-30^0C - 50Hz - 6dB contours')
hold off
colorbar

