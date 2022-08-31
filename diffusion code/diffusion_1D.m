function [Fo_all,Fo_new,timesteps,t_Fo,x,a,b]=diffusion_1D(dx,X,nr,years,a,b)

T=1200+273.15; %temperature in K
fug=10^(8.912-(25160/T)); %log fugacity in bars
fug_pa=fug.*100000; %fugacity in Pa
P=60*1e6; %pressure in Pa
R=8.31446; %gas constant in J/(mole K)

%calculate dx for profile

x=(0:dx:X*dx-dx)';
Fo_init=ones(X,1).*0.9;
Fo_low=0.81;

%calculate diffusivity for stability criterion

DcFo_i=(a.*Fo_low+b).*((10^-9.21).*((fug_pa/(10^-7)).^(1/6)).*(10.^(3.*(0.9-(Fo_low)))).*exp(-(201000+(P-1e5).*7e-6)/(R.*T)))*(1e12); %diffusivity along c-axis in um2/s
DaFo_i=DcFo_i./6; %diffusivity along a-axis in um2/s
DbFo_i=DaFo_i; %diffusivity along b-axis in um2/s
%DtFo_i=(DaFo_i.*(cos(deg2rad(alpha))^2))+(DbFo_i.*(cos(deg2rad(beta))^2))+(DcFo_i.*(cos(deg2rad(gamma))^2)); %diffusivity along traverse in um2/s

%calculate dt to fit stability criterion R<0.5
r=0.4;
dt=r.*((dx.^2)/max(DcFo_i));

%define duration of diffusion model

duration=years*365*24*60*60; %in sec

timesteps=ceil(duration./dt);
r_steps=timesteps./nr; %frequence des recharges
v=r_steps:r_steps:timesteps; %linear recharge rate
%v=sort(rand(nr,1)).*timesteps; 
v(nr)=timesteps;

%redefine concentration that will be renewed at every step
Fo_new=ones(X,nr).*Fo_init;

DcFo_new=ones(X,nr).*DcFo_i;

%Fo_all=zeros(X,timesteps);

t_Fo=zeros(1,timesteps);


%perform calculations, solve continuity equation
%Fo content
curr_nr=1;
curr_v=1;

%display a waitbar
wb = waitbar(0,'Diffusing...');
tic; %start the waitbar stopwatch

for j=1:timesteps
    for m=1:curr_nr
        
        Fo_old(:,m)=Fo_new(:,m);
        DcFo_old(:,m)=(a.*Fo_old(:,m)+b).*((10^-9.21).*((fug_pa/(10^-7)).^(1/6)).*(10.^(3.*(0.9-(Fo_old(:,m))))).*exp(-(201000+(P-1e5).*7e-6)/(R.*T)))*(1e12); %diffusivity along c-axis in um2/s
        DbFo_old=DcFo_old./6;
        DaFo_old=DcFo_old./6;
        
        for i=2:(X-1)
            Fo_new(i,m)=Fo_old(i,m) + (dt./dx.^2).*(((DcFo_old(i+1,m)-DcFo_old(i,m)).*(Fo_old(i+1,m)-Fo_old(i,m)))+DcFo_old(i,m).*(Fo_old(i+1,m)-2.*Fo_old(i,m)+Fo_old(i-1,m)));
        end
        
        Fo_new(1,m)=0.8; %boundary condition for first element
        
        Fo_new(X,m)= Fo_new(X-1,m); % boundary condition for last element
        Fo_all(j,m)=Fo_new(X,m);
        
    end
    t_Fo(j)=(dt.*j)./(3600.*24);
    
    if j>v(curr_v)
        curr_nr=curr_nr+1;
        curr_v=curr_v+1;
    end
    t=toc;
    Perc=double(j)./double(timesteps);
    Trem=t/Perc-t;
    Hrs=floor(Trem/3600);
    Min=floor((Trem-Hrs.*3600)/60);
    waitbar(Perc,wb,[sprintf('%0.1f',Perc*100) '%,'...
        sprintf('%03.0f', Hrs) ':'...
        sprintf('%02.0f', Min) ':'...
        sprintf('%02.0f',rem(Trem,60)) 'remaining'])
end

close (wb)


   