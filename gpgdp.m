a = 'D:\Dropbox\sim\output\';
b = 'data 7.txt';
c = strcat(a,b);
pars =importUData(c);
%[frame     tension     Rwallposx   lwally pos  frame   lwallosc  gamma    sigma    wallv]
frames      = data(:,1);
normF       = data(:,2);
Rwallposx   = data(:,3);
Lwallposy   = data(:,4);
LWallForce  = data(:,6);
gamma       = data(:,7);
sigma       = data(:,8);
wallv       = data(:,9);
%mas waitime woarena hoarena wall_percent walfrq num_body amp_max b k mot L

spf = .05;
p = pi/180;
mass = pars(1);
sf = pars(2);
w_oarena = pars(3);
h_oarena=pars(4);
compLen= pars(5);
wallfreq=pars(6);
num_body= pars(7);
amp = pars(8);
amp2=pars(8)/wallfreq;
b = pars(9);
k = pars(10);
if length(pars)>= 11
mot= pars(11);
L=pars(12);
else
    mot = '?';
    L='?';
end
% Rwallposx= abs(max(Rwallposx)-Rwallposx);
% RWallposx = w_oarena-Rwallposx;
%  str$(1/wall_freq)+"*cos("+str$(WALL_FREQ)+"*(t-"+str$(waitTime)+"))*(t>="+str$(waitTime)+")"
% v= 1/wallfreq*cos(wallfreq*(t-waittime))

%stress = sigma = F/A
%strain = gamma = deltax/h (for me h = (ArenaWidth at time=t))

arg = (wallfreq*(frames(sf:end)-(sf))*spf); 
analytic = zeros(length(LWallForce),1);
analytic(1:sf-1)=LWallForce(1:sf-1);
analytic(sf:end)=(LWallForce(sf:end)+mass*(amp*wallfreq/wallfreq)*sin(arg));

%these should all be area of the blob not of arena and shat
sigma = analytic./(h_oarena*w_oarena);
gamma = Lwallposy./Rwallposx;
gammadot = diff(gamma)/(spf*length(frames));
gprime = sigma./gamma;
gdprime= sigma(2:end)./gammadot;
strainrate=gammadot;
strain=gamma;
stress= sigma;


figure(1)
hold on;
dvdt=(diff(wallv)./(spf));
dvdt=(LWallForce(2:end)-dvdt.*mass);

plot(frames,analytic,'b')
plot(frames(2:end),dvdt,'r')
plot(frames,LWallForce,'g')

leg = {'analytic subtraction -f','dv/dt*m subtraction - f','f'};
legend(leg,0);
title('lwall tension vs. frames');
xlabel('frames');
ylabel('force');

figure(2)
hold all;
plot(strain(sf:end),stress(sf:end),'r');
plot(strainrate(sf:end),stress(sf:end-1),'b');
% plot(wallv(2:end),gdprime);

xlabel('Strain, Strain rate [1/s]', 'FontSize', 14);
ylabel('Stress [Pa]', 'FontSize', 14);
title('Lissajous plots','FontSize', 14);
leg = {'Stress vs. Strain','Stress vs. Strain Rate'};
legend(leg,0);

figure(4)
hold all;
plot(strain(sf+1:(sf+1/spf*wallfreq)),stress(sf+1:(sf+1/spf*wallfreq)),'r');
plot(strainrate(sf-1:sf-1+1/spf*wallfreq),stress(sf:sf+1/spf*wallfreq),'b');
% plot(wallv(2:end),gdprime);

xlabel('Strain, Strain rate [1/s]', 'FontSize', 14);
ylabel('Stress [Pa]', 'FontSize', 14);
title('Lissajous plots for 1 period', 'FontSize', 14);
leg = {'Stress vs. Strain','Stress vs. Strain Rate'};
legend(leg,0);


figure(3)
hold on;
dvdtmean = mean(abs(dvdt));
amean=mean(abs(analytic(sf+1:end)));
% amean = mean(find(abs(analytic(sf:end)));
% plot(frames,analytic,'b')
m = max(abs(dvdt(sf+1:end)));
m = max(abs(analytic(sf+2:end)));

% plot(frames(1:end),analytic,'r')
% plot(frames,Lwallposy/(max(Lwallposy))*amean,'g')
plot(frames,stress,'r')
plot(frames,strain,'g')
smoothed = smooth(frames,stress,0.1,'loess');
plot(frames,smoothed,'b')
leg = {'stress','strain','smoothed strain'};
legend(leg,0);
title(strcat('b=',num2str(b),' compAmt=',num2str(compLen),' \omega='...
    ,num2str(wallfreq), ' n=', num2str(num_body), ' m=',num2str(mass),...
    ' A=', num2str(amp2,3),' k=',num2str(k), ' mots=', num2str(mot),...
    ' L=', num2str(L)),...
    'FontSize', 12);
xlabel('frames');
% ylabel('force/position(normalized)');
ylabel('stress[Pa]/strain');
