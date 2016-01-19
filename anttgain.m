clc;
clear;
close all;
format long;
warning off;
%---Definition of constants and initialization---
PI = 4.0*atan(1.0); % pi
E = 120.0*PI; % eta ie intrinsic impedance
THETA = PI/180.0; % theta max value
UMAX = 0.0; % maximum radiation intensity
PRAD = 0.0; % power radiated
TOL = 1.0E-6; % tolerance
f1=57e9; % start frequency
f2=64e9; % stop frequency
freq=f1:0.1e9:f2; % frequency range
c=3e8; % speed of light
n=120*pi; % intrinsic impedance
k=2*pi; % 2 * pi
%**************************************************************************
%***************************MONOPOLE ANTENNA ELEMENT*************************
VSWR=zeros(length(freq),1); % initializing VSWR array for storing results
gain=zeros(length(freq),1);% initializing gain array for storing results for sweep response
S11=zeros(length(freq),1); % initliatizing S11 array for storing results for sweep response
for i=1:length(freq)
lamda=c/freq(i); % wavelength
L=lamda/4; % quarter wavelength monopole antenna element
r=0.2*L; % radius of monopole antenna element
A = L*PI*134; % Maximum Effective apperture area of monopole antenna element
I = 1;
while(I <= 180) % Evaluating for each theta value
XI = I*PI/180.0; % Instantaneous theta value
if(XI ~= PI)
U = ((cos(A*cos(XI))-cos(A))/sin(XI))^2*(E/(8.0*PI^2)); % Instantaneous value of Radiation intensity
if(U > UMAX)
UMAX = U; % Maximum radiation intensity
end
end
UA = U*sin(XI)*THETA*2.0*PI; % Radiation intensity a factor of theta
PRAD = PRAD+UA; % total power radiated
I = I+1;
end
RR = 2.0*PRAD; % Radiation resistance due to maximum current
if(A ~= PI)
RIN = RR/(sin(A))^2; % Input Resistance
end
Xm=30*(2*sinint(k*L)+cos(k*L)*(2*sinint(k*L)-sinint(2*k*L))- ... % Reactance due to maximum current
sin(k*L)*(2*cosint(k*L)-cosint(2*k*L)-cosint(2*k*r^2/L)));
Xin=Xm/(sin(k*L/2))^2; % Input reactance
Zin=complex(RR,Xm); % Input Impedance
gam=abs((Zin-50)/(Zin+50)); % reflection Coefficient
Vsw=(1+gam)/(1-gam); % VSWR value for each frequency in the sweep
VSWR(i)=Vsw; % VSWR value for entire sweep
s1=-20*log10((Vsw+1)/(Vsw-1)); % Return loss value for each frequency in the sweep
S11(i)=s1; % Return Loss value for entire sweep
end
plot(freq/1e9,VSWR) % Plot of frequency vs VSWR
grid('on')
title('VSWR Vs Frequency Plot of Vertical Monopole Antenna');
xlabel('Freq in GHz');
ylabel('VSWR')
figure
plot(freq/1e9,S11) % Plot of frequency vs Return Loss
grid('on')
title('Return Loss Vs Frequency Plot of Vertical Monopole Antenna');
xlabel('Freq in GHz');
ylabel('S11')
%**************************************************************************
%**************************************************************************
UMAX = 0.0; % Maximum Radiation Intensity
PRAD = 0.0; % Power Radiated
for i=1:length(freq)
lamda=c/freq(i); % Wavelength
L=lamda/4; % quarter wavelength monopole element
r=0.1*lamda; % radius of monopole antenna element
A = L*PI; % effective apperture
I = 1;
while(I <= 180)
XI = I*PI/180.0;
if(XI ~= PI)
U = ((cos(A*cos(XI))-cos(A))/sin(XI))^2*(E/(8.0*PI^2)); % Radiation Intensity
if(U > UMAX)
UMAX = U; % Maximum Radiation Intensity
end
end
UA = U*sin(XI)*THETA*2.0*PI; % Radiation intesity factor of theta
PRAD = PRAD+UA; % total power radiated
I = I+1;
end
D = (4.0*PI*UMAX)/PRAD; % directivity
DDB = 10.0*log10(D); % in dBs
gain(i)=DDB; % antenna gain
UMAX = 0.0;
PRAD = 0.0;
end
figure
plot(freq/1e9,gain) % frequency vs antenna gain plot
grid('on')
title('Gain Vs Frequency Plot of Vertical Monopole Antenna');
xlabel('Freq in GHz');
ylabel('Gain (dBs)')
%**************************************************************************
%**************************************************************************
%---Calculation of elevation far-field patterns in 1 degree increments---
T = zeros(180,1);
ET = zeros(180,1);
EdB = zeros(180,1);
x = 1;
while(x<=180)
T(x) = x-0.99;
ET(x) = (cos(PI*L*cos(T(x)*THETA))-cos(PI*L))/sin(T(x)*THETA); % Intensity in theta direction
x = x+1;
end
ET = abs(ET); % absolute value
ETmax = max(abs(ET)); % maximum intensity
EdB = 20*log10(abs(ET)/ETmax); % gain in dbs
T=T'; EdB=EdB';
EdB=[EdB fliplr(EdB)];
% expanding U to span entire space
EdB=EdB';
U=zeros(360,180);
for k=1:180
U(:,k)=EdB(:,1); % gain in phi orientation
end
T=[T T+180];
az=-180:1:180; % azimuth coverage
el=-90:1:90; % elevation coverage
fmax=60e9;
freqvector=[59 60].*1e9; % frequency vector
UV=[U U EdB];
UU=UV(1:181,1:361); % gain in 3D space
%**************************************************************************
%**************************************************************************
% Antenna Custom Element using the phased array toolbox, here azimuth and
% elevation range is defined along with the gain values calculated along
% theta and phi planes in 3D
sAnt=phased.CustomAntennaElement('FrequencyVector',freqvector,'AzimuthAngles',az,'ElevationAngles',el,'RadiationPattern',UU);
figure
% plotting 3D gain response of single vertical antenna element in dBs
plotResponse(sAnt,fmax,'Format','Polar','RespCut','3D');
title('3D RESPONSE OF MONOPOLE','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D response of single vertical monopole antenna element in azimuth plane in dBs
figure
subplot(1,2,1)
plotResponse(sAnt,fmax,'Format','Polar','Respcut','az','Cutangle',90);
title('2D RESPONSE OF MONOPOLE IN AZIMUTH dBs','FontSize',10,'FontWeight','bold','Color','b');
subplot(1,2,2)
%plotting 2D response of single vertical monopole antenna element in elevation plane in dBs
plotResponse(sAnt,fmax,'Format','Polar','Respcut','el');
title('2D RESPONSE OF MONOPOLE IN ELEVATION dBs','FontSize',10,'FontWeight','bold','Color','b');
%**************************************************************************
%**************************************************************************
figure
%plotting 3D magnitude response of single vertical monopole antenna element
%(dimensionless)
plotResponse(sAnt,fmax,'RespCut','3D','Format','Polar','Unit','mag');
title('3D GAIN RESPONSE OF MONOPOLE ANTENNA (DIMENSIONLESS)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D magnitude response of single vertical monopole antenna element
%in azimuth plane(dimensionless)
figure
subplot(1,2,1)
plotResponse(sAnt,fmax,'RespCut','Az','Format','Polar','Unit','mag');
title('2D RESPONSE OF MONOPOLE IN AZIMUTH (DIMENSIONLESS)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D magnitude response of single vertical monopole antenna element
%in elevation plane(dimensionless)
subplot(1,2,2)
plotResponse(sAnt,fmax,'RespCut','El','Format','Polar','Unit','mag');
title('2D RESPONSE OF MONOPOLE IN ELEVATION (DIMENSIONLESS)','FontSize',10,'FontWeight','bold','Color','b');
figure
plotResponse(sAnt,fmax,'RespCut','El');
%**************************************************************************
%***************1X4 LINEAR ARRAY*************************************
fmax=freqvector(end);
c=3e8;
lamda =c/fmax;
%defining uniform linear array of 1x4 vertical monopole antenna elements
%designed earlier with the spacing of 0.5 lambda
sArray= phased.ULA('Element',sAnt,'NumElements',4,'ElementSpacing',0.5* lamda);
figure
%plotting 3d response of 1x4 antenna array in dBs
plotResponse(sArray,fmax,c,'RespCut','3D','Format','Polar');
title('3D GAIN RESPONSE OF 1X4 ELEMENTS ARRAY (dBs)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2d azimuth response of 1x4 antenna array in dBs
figure
subplot(1,2,1)
plotResponse(sArray,fmax,c,'RespCut','Az','Format','Polar');
title('2D RESPONSE OF 1x4 ELEMENTS ARRAY IN AZIMUTH PLANE (dBs)','FontSize',10,'FontWeight','bold','Color','b');
subplot(1,2,2)
%plotting 2d elevation response of 1x4 antenna array in dBs
plotResponse(sArray,fmax,c,'RespCut','El','Format','Polar');
title('2D RESPONSE OF 1x4 ELEMENTS ARRAY IN ELEVATION PLANE (dBs)','FontSize',10,'FontWeight','bold','Color','b');
figure
%view 1x4 linear array
viewArray(sArray,'ShowIndex','All');
%**************************************************************************
%**************************************************************************
figure
%plotting 3D magnitude response of 1x4 antenna array
plotResponse(sArray,fmax,c,'RespCut','3D','Format','Polar','Unit','mag');
title('3D GAIN RESPONSE OF 1X4 ELEMENTS ARRAY (DIMENSIONLESS)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D magnitude response of 1x4 antenna array in azimuth plane
figure
subplot(1,2,1)
plotResponse(sArray,fmax,c,'RespCut','Az','Format','Polar','Unit','mag');
title('2D RESPONSE OF 1x4 ELEMENTS ARRAY IN AZIMUTH PLANE (DIMENSIONLESS)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D magnitude response of 1x4 antenna array in elevation plane
subplot(1,2,2)
plotResponse(sArray,fmax,c,'RespCut','El','Format','Polar','Unit','mag');
title('2D RESPONSE OF 1x4 ELEMENTS ARRAY IN ELEVATION PLANE (DIMENSIONLESS)','FontSize',10,'FontWeight','bold','Color','b');
%**************************************************************************
%**************************************************************************
%plotting 2D gain response of 1x4 antenna array on rectangular plot for
%beam realization
figure
plotResponse(sArray,fmax,c)
title('2D RESPONSE OF 1x4 ELEMENTS ARRAY(dBs)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D gain response of 1x4 antenna array in azimuth plane on rectangular plot for
%beam realization
figure
subplot(1,2,1)
plotResponse(sArray,fmax,c,'RespCut','Az');
title('2D RESPONSE OF 1x4 ELEMENTS ARRAY IN AZIMUTH PLANE(dBs)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D gain response of 1x4 antenna array in elevation plane on rectangular plot for
%beam realization
subplot(1,2,2)
plotResponse(sArray,fmax,c,'RespCut','El');
title('2D RESPONSE OF 1x4 ELEMENTS ARRAY IN ELEVATION PLANE(dBs)','FontSize',10,'FontWeight','bold','Color','b');
%**************************************************************************
%***********************4X4 RECTANGULAR ARRAY*******************************
%defining 4x4 rectangular array of monopole antenna elements placed at a
%distance of lambda/2
sArray = phased.URA('Element',sAnt,'Size',[4 4],'ElementSpacing',0.5* lamda);
ura=[4 4];
%plotting 3D gain response of antenna array in dBs
figure
plotResponse(sArray,fmax,c,'RespCut','3D','Format','Polar');
title('3D GAIN PLOT OF 4x4 ELEMENT ARRAY (dBs)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D gain response of 4x4 antenna array in azimuth plane in dBs
figure
subplot (1,2,1)
plotResponse(sArray,fmax,c,'RespCut','Az','Format','Polar');
title('2D RESPONSE OF 4x4 ELEMENT ARRAY IN AZIMUTH (dBs)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D gain response of 4x4 antenna array in elevation plane in dBs
subplot (1,2,2)
plotResponse(sArray,fmax,c,'RespCut','El','Format','Polar');
title('2D RESPONSE OF 4X4 ELEMENT ARRAY IN ELEVATION PLANE (dBs)','FontSize',10,'FontWeight','bold','Color','b');
% visualization of array geometry
figure
viewArray(sArray,'ShowIndex','All');
%**************************************************************************
%**************************************************************************
figure
%plotting 3D magnitude response of 4x4 antenna array
plotResponse(sArray,fmax,c,'RespCut','3D','Format','Polar','Unit','mag');
title('3D GAIN RESPONSE OF 4X4 ELEMENTS ARRAY (DIMENSIONLESS)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D magnitude response of 4x4 antenna array in azimuth plane
figure
subplot(1,2,1)
plotResponse(sArray,fmax,c,'RespCut','Az','Format','Polar','Unit','mag');
title('2D RESPONSE OF 4x4 ELEMENTS ARRAY IN AZIMUTH PLANE (DIMENSIONLESS)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D magnitude response of 4x4 antenna array in elevation plane
subplot(1,2,2)
plotResponse(sArray,fmax,c,'RespCut','El','Format','Polar','Unit','mag');
title('2D RESPONSE OF 4x4 ELEMENTS ARRAY IN ELEVATION PLANE (DIMENSIONLESS)','FontSize',10,'FontWeight','bold','Color','b');
%**************************************************************************
%**************************************************************************
%plotting 2D gain response of 4x4 antenna array on rectangular plot for
%beam realization
figure
plotResponse(sArray,fmax,c)
title('2D RESPONSE OF 4x4 ELEMENTS ARRAY(dBs)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D gain response of 4x4 antenna array in azimuth plane on rectangular plot for
%beam realization
figure
subplot(1,2,1)
plotResponse(sArray,fmax,c,'RespCut','Az');
title('2D RESPONSE OF 4x4 ELEMENTS ARRAY IN AZIMUTH PLANE(dBs)','FontSize',10,'FontWeight','bold','Color','b');
%plotting 2D gain response of 4x4 antenna array in elevation plane on rectangular plot for
%beam realization
subplot(1,2,2)
plotResponse(sArray,fmax,c,'RespCut','El');
title('2D RESPONSE OF 4x4 ELEMENTS ARRAY IN ELEVATION PLANE(dBs)','FontSize',10,'FontWeight','bold','Color','b');
%**************************************************************************
%****************BEAM STEERING PART**********************************
r= input('Start Phi Angle = ');
q= input('End Phi Angle = ');
u= input('No. of Scans = ');
scanPhi = r:q;
global hScope;
for x=1:u
hScan = phased.SteeringVector('SensorArray',sArray,...
'PropagationSpeed',c);
hResponse = phased.ArrayResponse('WeightsInputPort',true,...
'SensorArray',sArray);
scanTheta = zeros(1,numel(scanPhi));
scanAngles = [scanPhi;scanTheta];
weights = step(hScan,fmax,scanAngles); % Calculate Weights for steering
az = -180:180;
el = zeros(1,numel(az));
ang_pairs = [az;el];
if isempty(hScope)
hScope = ArrayResponseDemo2DPolarScope; % Initialize scope
end
for t = 1:length(scanTheta) % Calculate response
w = weights(:,t);
temp_out = abs(step(hResponse,fmax,ang_pairs,w));
temp_out = temp_out/max(temp_out);
% temp_out=20*log10 (temp_out);
pause(0.1);
step(hScope,temp_out); % Display on scope
%disp(w);
m{t}=w;
end
%disp(m)
end
