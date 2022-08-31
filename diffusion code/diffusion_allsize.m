clear all

%input starting conditions
X=50;
nr=1;
years=0.1;

[Fo_all1,Fo_new1,timesteps1,t_Fo1,x1,a1,b1]=diffusion_1D(2.5,X,nr,years,-0.134,13.432);
[Fo_all2,Fo_new2,timesteps2,t_Fo2,x2,a2,b2]=diffusion_1D(5,X,nr,years,-0.1353,13.564);
[Fo_all3,Fo_new3,timesteps3,t_Fo3,x3,a3,b3]=diffusion_1D(10,X,nr,years,-0.1402,14.038);


n1=2.*nr;
Fo_all1=reshape(ones(n1,1).*Fo_all1(timesteps1,:),n1.*nr,1);
n2=2.*nr;
Fo_all2=reshape(ones(n2,1).*Fo_all2(timesteps2,:),n2.*nr,1);
n3=2.*nr;
Fo_all3=reshape(ones(n3,1).*Fo_all3(timesteps3,:),n3.*nr,1);


Fo_low=0.8;
Fo_high=0.9;
n_pop=200;
Fo_dist=normrnd(0.81, 0.01,[n_pop,1]);
Fo_pop=ones(X,n_pop);


for i=1:n_pop
    Fo_init=ones(X,1).*Fo_dist(i);
    Fo_pop(:,i)=Fo_init;
end
[Fo_all4,Fo_new4,timesteps4,t_Fo4,x4,a4,b4]=diffusion_1D_pop(10,X,Fo_pop,n_pop,years,-0.1402,14.038);

Fo_all4=Fo_all4(timesteps4,:);

% Fo_all=[Fo_all1
%     Fo_all2
%     Fo_all3
%     Fo_all4'];

edges=[0.8:0.005:0.9];

hist1=histc(Fo_all1,edges);
hist2=histc(Fo_all2,edges);
hist3=histc(Fo_all3,edges);
hist4=histc(Fo_all4,edges);
hist4=hist4';
hist_all=[hist1 hist2 hist3 hist4];

figure(1)

% bar(edges,hist_all,'stacked'); hold on
h = bar(edges,hist_all,'stacked');
set(h, {'DisplayName'}, {'250um','500um','1mm','pre-existing ol'}')
legend() 

