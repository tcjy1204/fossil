%% coal Stated Policies cost: carbon
clear;clc;
tic
path=('F:\paper\science\data\coal.xlsx');
path1=('F:\paper\science\data\futureproduction.xlsx');
pathtcs=('F:\paper\science\results\tax\coal\STEPS.xlsx');
coal_d0=xlsread(path,'23','d2:i1266');% emission/yr,coal output/yr, and workforce of each coal mine
latlon=xlsread(path,'23','b2:c1266');
procode=xlsread(path,'23','aj2:aj1266');
l=size(coal_d0,1);
t=ones(l,1);% 1 for sites with output
coal_d=[coal_d0,t,coal_d0(:,2),latlon,procode];%11 column
coal_d2=coal_d;
coal_d3=coal_d;
coal_ew=xlsread(path1,'coal','k5:k31');% future coal production
coal_ew2=coal_ew-910.7473617;% exclude coal that are not covered by inventory
coal_r=-(coal_ew2-sum(coal_d(:,2)));%reduction in inventory needed
per=2024:1:2050;% period span
lper=size(per,2);
pei=xlsread(path,'23','m2:m1266');%emission inensity of production tCO2e/t
% sorted by output
o_nc_o=zeros(l,lper);% output reduction by natural closure
o_pf_o=zeros(l,lper);% output reduction by carbon tax-driven closure
o_fo_o=zeros(l,lper);% output reduction by directly forced closure
e_nc_o=zeros(l,lper);% emission reduction by natural closure
e_pf_o=zeros(l,lper);% emission reduction by carbon tax-driven closure
e_fo_o=zeros(l,lper);% emission reduction by directly forced closure
op=xlsread(path,'23','j2:j1266');%original profit without carbon tax
cart=xlsread(path1,'carbontax','j5:k31');%carbon tax 2030-2050
c=xlsread(path,'23','k2:k1266');%original cost without carbon tax and employment
pf=zeros(l,lper);% profit with carbon tax for production emission
spf=zeros(l,lper);% profit with carbon tax for all emission
tpf=zeros(l,lper);
coalemo=zeros(l,lper);
coalomo=zeros(l,lper);
coallmo=zeros(l,lper);
coalpfmo=zeros(l,lper);
coalomob=zeros(l,lper);
lato=zeros(l,lper);
lono=zeros(l,lper);
procodeo=zeros(l,lper);
for i=1:lper % 24-50
    pf(:,i)=op-pei.*cart(i,1);
    spf(:,i)=coal_d(:,5).*cart(i,1);%social cost considering all emissions
    tpf(:,i)=spf(:,i).*coal_d(:,2);% total profit M CNY
    coal_d=[coal_d2, tpf(:,i)];% 12 col
    for j=1:l
        if coal_d(j,6)<=i% no resources
            o_nc_o(j,i)=coal_d2(j,2);
            e_nc_o(j,i)=coal_d2(j,1);
            coal_d(j,7)=0;
            coal_d(j,1)=0;
            coal_d(j,2)=0;
            coal_d(j,3)=0;
            coal_d(j,12)=0;
        end
        if pf(j,i)<=0 && coal_d(j,6)>i% stranded by carbon tax
            o_pf_o(j,i)=coal_d2(j,2);
            e_pf_o(j,i)=coal_d2(j,1);
            coal_d(j,7)=0;
            coal_d(j,1)=0;
            coal_d(j,2)=0;
            coal_d(j,3)=0;
            coal_d(j,12)=0;
        end
    end
    coal_d_so=sortrows(coal_d,[7 8 -4]);% sorted by status, output and depth
    coalemo(:,i)=coal_d_so(:,1);% emission matrix sorted by status, output and depth
    coalomo(:,i)=coal_d_so(:,2);% output matrix sorted by status, output and depth
    coallmo(:,i)=coal_d_so(:,3);% labor matrix sorted by status, output and depth
    coalpfmo(:,i)=coal_d_so(:,12);% profit matrix sorted by status, output and depth
    lato(:,i)=coal_d_so(:,9);
    lono(:,i)=coal_d_so(:,10);
    procodeo(:,i)=coal_d_so(:,11);
    coalemo0=coalemo;
    coalomo0=coalomo;
    coallmo0=coallmo;
    coalpfmo0=coalpfmo;
    coalomob(:,i)=coal_d_so(:,8);
    coal_r(coal_r>sum(coalomob(:,i)))=sum(coalomob(:,i));% reduction cannot beyond inventory
    for j=1:l
        if sum(coalomob(1:j,i))>=coal_r(i) && sum(coalomob(1:(j-1),i))<coal_r(i)
            coalemo(1:(j-1),i)=0;
            coalomo(1:(j-1),i)=0;
            coallmo(1:(j-1),i)=0;
            coalpfmo(1:(j-1),i)=0;
            coalomo(j,i)=sum(coalomob(1:j,i))-coal_r(i);
            coalemo(j,i)=coalemo0(j,i)./coalomob(j,i).*(sum(coalomob(1:j,i))-coal_r(i));
            coallmo(j,i)=coallmo0(j,i)./coalomob(j,i).*(sum(coalomob(1:j,i))-coal_r(i));
            coalpfmo(j,i)=coalpfmo0(j,i)./coalomob(j,i).*(sum(coalomob(1:j,i))-coal_r(i));           
        end
    end
    for j=1:l
        if coalemo(j,i)==0 && coal_d_so(j,7)==1
            o_fo_o(j,i)=coalomob(j,i);
            e_fo_o(j,i)=coalemo0(j,i);
        end
        if coalemo(j,i)~=0 && coalemo(j,i)<coalemo0(j,i)
            o_fo_o(j,i)=(coal_r(i)-sum(coalomob(1:j-1,i)));
            e_fo_o(j,i)=coalemo0(j,i)./coalomob(j,i).*o_fo_o(j,i);
        end
    end
    %
end
coal_so_ye=sum(coalemo,1)/1000;% yearly emissionGtCO2e
coal_so_ye_cum=cumsum(coal_so_ye);% cumulative yearly emission
coal_so_yo=sum(coalomo,1);% yearly output Mt





% sorted by emission intensity
%for carbon tax
o_nc_e=zeros(l,lper);% output reduction by natural closure
o_pf_e=zeros(l,lper);% output reduction by carbon tax-driven closure
o_fo_e=zeros(l,lper);% output reduction by directly forced closure
e_nc_e=zeros(l,lper);% emission reduction by natural closure
e_pf_e=zeros(l,lper);% emission reduction by carbon tax-driven closure
e_fo_e=zeros(l,lper);% emission reduction by directly forced closure
op=xlsread(path,'23','j2:j1266');%original cost without carbon tax
cart=xlsread(path1,'carbontax','j5:k31');%carbon tax 2030-2050
c=xlsread(path,'23','k2:k1266');%original cost without carbon tax and employment
pf=zeros(l,lper);
spf=zeros(l,lper);% profit with carbon tax for all emission
tpf=zeros(l,lper);
%
coaleme=zeros(l,lper);
coalome=zeros(l,lper);
coallme=zeros(l,lper);
coalpfme=zeros(l,lper);
coalomeb=zeros(l,lper);
late=zeros(l,lper);
lone=zeros(l,lper);
procodee=zeros(l,lper);
coal_d=[coal_d0,t,coal_d0(:,2),latlon,procode];
for i=1:lper % 24-50
    %for carbon tax
    pf(:,i)=op-pei.*cart(i,1);
    spf(:,i)=coal_d(:,5).*cart(i,1);%social cost considering all emissions
    tpf(:,i)=spf(:,i).*coal_d(:,2);% total profit M CNY
    coal_d=[coal_d2, tpf(:,i)];% 12 col
    %
    for j=1:l
        if coal_d(j,6)<=i% no resources
            %for carbon tax
            o_nc_e(j,i)=coal_d2(j,2);
            e_nc_e(j,i)=coal_d2(j,1);
            coal_d(j,7)=0;
            coal_d(j,1)=0;
            coal_d(j,2)=0;
            coal_d(j,3)=0;
            coal_d(j,12)=0;
        end
        %for carbon tax
        if pf(j,i)<=0 && coal_d(j,6)>i% stranded by carbon tax
            o_pf_e(j,i)=coal_d2(j,2);
            e_pf_e(j,i)=coal_d2(j,1);
            coal_d(j,7)=0;
            coal_d(j,1)=0;
            coal_d(j,2)=0;
            coal_d(j,3)=0;
            coal_d(j,12)=0;           
        end
    end
    coal_d_se=sortrows(coal_d,[7 -5 -4]);% sorted by status, emission intensity and depth
    coaleme(:,i)=coal_d_se(:,1);% emission matrix sorted by status, output and depth
    coalome(:,i)=coal_d_se(:,2);% output matrix sorted by status, output and depth
    coallme(:,i)=coal_d_se(:,3);% labor matrix sorted by status, output and depth
    coalpfme(:,i)=coal_d_se(:,12);% profit matrix sorted by status, output and depth
    late(:,i)=coal_d_se(:,9);
    lone(:,i)=coal_d_se(:,10);
    procodee(:,i)=coal_d_se(:,11);
    coaleme0=coaleme;
    coalome0=coalome;
    coallme0=coallme;
    coalpfme0=coalpfme;
    coalomeb(:,i)=coal_d_se(:,8);
    coal_r(coal_r>sum(coalomeb(:,i)))=sum(coalomeb(:,i));% reduction cannot beyond inventory
    for j=1:l
        if sum(coalomeb(1:j,i))>=coal_r(i) && sum(coalomeb(1:(j-1),i))<coal_r(i)
            coaleme(1:(j-1),i)=0;
            coalome(1:(j-1),i)=0;
            coallme(1:(j-1),i)=0;
            coalpfme(1:(j-1),i)=0;
            coalome(j,i)=sum(coalomeb(1:j,i))-coal_r(i);
            coaleme(j,i)=coaleme0(j,i)./coalomeb(j,i).*(sum(coalomeb(1:j,i))-coal_r(i));
            coallme(j,i)=coallme0(j,i)./coalomeb(j,i).*(sum(coalomeb(1:j,i))-coal_r(i));
            coalpfme(j,i)=coalpfme0(j,i)./coalomeb(j,i).*(sum(coalomeb(1:j,i))-coal_r(i));           
        end
    end
    %for carbon tax
    for j=1:l
        if coaleme(j,i)==0 && coal_d_se(j,7)==1
            o_fo_e(j,i)=coalomeb(j,i);
            e_fo_e(j,i)=coaleme0(j,i);
        end
        if coaleme(j,i)~=0 && coaleme(j,i)<coaleme0(j,i)
            o_fo_e(j,i)=(coal_r(i)-sum(coalomeb(1:j-1,i)));
            e_fo_e(j,i)=coaleme0(j,i)./coalomeb(j,i).*o_fo_e(j,i);
        end
    end
    %
end
coal_se_ye=sum(coaleme,1)/1000;% yearly emission GtCO2e
coal_se_ye_cum=cumsum(coal_se_ye);% cumulative yearly emission
coal_se_yo=sum(coalome,1);% yearly output Mt

x=2024:1:2050;
plot(x,coal_so_ye,'-','color',[235, 23, 23]./255,'Markersize',6,'LineWidth',5)
hold on
plot(x,coal_se_ye,'-o','color',[235, 23, 23]./255,'Markersize',6,'LineWidth',5)
ylabel('Coal (GtCO_2e)','FontSize',50);
xlim([2024 2050]);
ylim([0 10]);
xticks([2025,2030,2035,2040,2045,2050]);
set(gca,'FontSize',100);

diffoe_cs=coal_so_ye_cum(1,size(coal_so_ye_cum,2))-coal_se_ye_cum(1,size(coal_se_ye_cum,2));
% output
xlswrite(pathtcs,coalomo,'coalomo','a1');
xlswrite(pathtcs,coalemo,'coalemo','a1');
xlswrite(pathtcs,coallmo,'coallmo','a1');
xlswrite(pathtcs,coalpfmo,'coalpfmo','a1');
xlswrite(pathtcs,coalome,'coalome','a1');
xlswrite(pathtcs,coaleme,'coaleme','a1');
xlswrite(pathtcs,coallme,'coallme','a1');
xlswrite(pathtcs,coalpfme,'coalpfme','a1');
xlswrite(pathtcs,coal_so_yo,'coal_so_yo','a1');
xlswrite(pathtcs,coal_so_ye,'coal_so_ye','a1');
xlswrite(pathtcs,coal_se_yo,'coal_se_yo','a1');
xlswrite(pathtcs,coal_se_ye,'coal_se_ye','a1');
xlswrite(pathtcs,o_nc_o,'o_nc_o','a1');
xlswrite(pathtcs,o_pf_o,'o_pf_o','a1');
xlswrite(pathtcs,o_fo_o,'o_fo_o','a1');
xlswrite(pathtcs,o_nc_e,'o_nc_e','a1');
xlswrite(pathtcs,o_pf_e,'o_pf_e','a1');
xlswrite(pathtcs,o_fo_e,'o_fo_e','a1');
xlswrite(pathtcs,e_nc_o,'e_nc_o','a1');
xlswrite(pathtcs,e_pf_o,'e_pf_o','a1');
xlswrite(pathtcs,e_fo_o,'e_fo_o','a1');
xlswrite(pathtcs,e_nc_e,'e_nc_e','a1');
xlswrite(pathtcs,e_pf_e,'e_pf_e','a1');
xlswrite(pathtcs,e_fo_e,'e_fo_e','a1');
xlswrite(pathtcs,lato,'lato','a1');
xlswrite(pathtcs,lono,'lono','a1');
xlswrite(pathtcs,procodeo,'procodeo','a1');
xlswrite(pathtcs,late,'late','a1');
xlswrite(pathtcs,lone,'lone','a1');
xlswrite(pathtcs,procodee,'procodee','a1');

xlswrite(pathtcs,coal_r','overview','b11');

pceo=[zeros(31,27),(1:1:31)'];
for i=1:27
    for j=1:l
        for g=1:31
            if procodeo(j,i)==pceo(g,28)
                pceo(g,i)=pceo(g,i)+coalemo(j,i);
            end
        end
    end
end
tpceo=sum(pceo(:,1:27),2);
xlswrite(pathtcs,tpceo,'tpceo','a1');

pcee=[zeros(31,27),(1:1:31)'];
for i=1:27
    for j=1:l
        for g=1:31
            if procodee(j,i)==pcee(g,28)
                pcee(g,i)=pcee(g,i)+coaleme(j,i);
            end
        end
    end
end
tpcee=sum(pcee(:,1:27),2);
xlswrite(pathtcs,tpcee,'tpcee','a1');


%composition coal_so_yo
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(o_nc_o);sum(o_pf_o);sum(o_fo_o)]';
area(x,y);
ylabel0=['STEPS' newline 'output reduction-Coal (Bt)'];
ylabel(ylabel0)
xlim([2024 2050]);
legend({'Natural','Carbon tax','Forced'},'Location','northwest')
set(gca,'FontSize',50);
print('F:\paper\science\figure\output reduction composition tax steps coal.tiff', '-dtiff', ['-r' num2str(dpi)]);

%composition coal_so_ye
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(e_nc_o);sum(e_pf_o);sum(e_fo_o)]'./1000;
area(x,y);
legend({'Natural','Carbon tax','Forced'},'Location','northwest')
ylabel0=['STEPS-PSO' newline 'emission reduction-Coal (Mt)'];
ylabel(ylabel0)
xlim([2024 2050]);
set(gca,'FontSize',50);
print('F:\paper\science\figure\emission reduction composition tax steps Scale coal.tiff', '-dtiff', ['-r' num2str(dpi)]);

%composition coal_se_ye
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(e_nc_e);sum(e_pf_e);sum(e_fo_e)]'./1000;
area(x,y);
legend({'Natural','Carbon tax','Forced'},'Location','northwest')
ylabel0=['STEPS-PHEI' newline 'emission reduction-Coal (Mt)'];
ylabel(ylabel0)
xlim([2024 2050]);
set(gca,'FontSize',50);
print('F:\paper\science\figure\emission reduction composition tax steps Emission coal.tiff', '-dtiff', ['-r' num2str(dpi)]);

%emission intensity tco2e/t
eiso=round(coal_so_ye,6).*1000./round(coal_so_yo,6);
eise=round(coal_se_ye,6).*1000./round(coal_se_yo,6);


%emission intensity tco2e/t
eiso=coal_so_ye.*1000./coal_so_yo;
eise=coal_se_ye.*1000./coal_se_yo;



f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
x=2024:1:2050;
plot(x,coal_so_ye,'-','color',[235, 23, 23]./255,'Markersize',6,'LineWidth',5)
hold on
plot(x,coal_se_ye,'-o','color',[235, 23, 23]./255,'Markersize',6,'LineWidth',5)
hold on






%% coal Announced Pledges
path=('F:\paper\science\data\coal.xlsx');
path1=('F:\paper\science\data\futureproduction.xlsx');
pathtca=('F:\paper\science\results\tax\coal\APS.xlsx');
coal_d0=xlsread(path,'23','d2:i1266');% emission/yr,coal output/yr, and workforce of each coal mine
l=size(coal_d0,1);
t=ones(l,1);% 1 for sites with output
coal_d=[coal_d0,t,coal_d0(:,2),latlon,procode];
coal_d2=coal_d;
coal_d3=coal_d;
coal_ew=xlsread(path1,'coal','l5:l31');% future coal production
coal_ew2=coal_ew-910.7473617;% exclude coal that are not covered by inventory
coal_r=-(coal_ew2-sum(coal_d(:,2)));%reduction in inventory needed
per=2024:1:2050;% period span
lper=size(per,2);
%for carbon tax
pei=xlsread(path,'23','m2:m1266');%emission inensity of production tCO2e/t
%
% sorted by output
%for carbon tax
o_nc_o=zeros(l,lper);% output reduction by natural closure
o_pf_o=zeros(l,lper);% output reduction by carbon tax-driven closure
o_fo_o=zeros(l,lper);% output reduction by directly forced closure
e_nc_o=zeros(l,lper);% emission reduction by natural closure
e_pf_o=zeros(l,lper);% emission reduction by carbon tax-driven closure
e_fo_o=zeros(l,lper);% emission reduction by directly forced closure
op=xlsread(path,'23','j2:j1266');%original profit without carbon tax
cart=xlsread(path1,'carbontax','j5:k31');%carbon tax 2030-2050
c=xlsread(path,'23','k2:k1266');%original cost without carbon tax and employment
pf=zeros(l,lper);% output reduction by directly forced closure
spf=zeros(l,lper);% profit with carbon tax for all emission
tpf=zeros(l,lper);
%
coalemo=zeros(l,lper);
coalomo=zeros(l,lper);
coallmo=zeros(l,lper);
coalpfmo=zeros(l,lper);
coalomob=zeros(l,lper);
lato=zeros(l,lper);
lono=zeros(l,lper);
procodeo=zeros(l,lper);

for i=1:lper % 24-50
    %for carbon tax
    pf(:,i)=op-pei.*cart(i,2);
    spf(:,i)=coal_d(:,5).*cart(i,2);%social cost considering all emissions
    tpf(:,i)=spf(:,i).*coal_d(:,2);% total profit M CNY
    coal_d=[coal_d2, tpf(:,i)];% 12 col
    %
    for j=1:l
        if coal_d(j,6)<=i% no resources
            %for carbon tax
            o_nc_o(j,i)=coal_d2(j,2);
            e_nc_o(j,i)=coal_d2(j,1);
            %
            coal_d(j,7)=0;
            coal_d(j,1)=0;
            coal_d(j,2)=0;
            coal_d(j,3)=0;
            coal_d(j,12)=0;           
        end
        %for carbon tax
        if pf(j,i)<=0 && coal_d(j,6)>i% stranded by carbon tax
            o_pf_o(j,i)=coal_d2(j,2);
            e_pf_o(j,i)=coal_d2(j,1);
            coal_d(j,7)=0;
            coal_d(j,1)=0;
            coal_d(j,2)=0;
            coal_d(j,3)=0;
            coal_d(j,12)=0;           
        end
        %
    end
    coal_d_so=sortrows(coal_d,[7 8 -4]);% sorted by status, output and depth
    coalemo(:,i)=coal_d_so(:,1);% emission matrix sorted by status, output and depth
    coalomo(:,i)=coal_d_so(:,2);% output matrix sorted by status, output and depth
    coallmo(:,i)=coal_d_so(:,3);% labor matrix sorted by status, output and depth
    coalpfmo(:,i)=coal_d_so(:,12);% profit matrix sorted by status, output and depth
    lato(:,i)=coal_d_so(:,9);
    lono(:,i)=coal_d_so(:,10);
    procodeo(:,i)=coal_d_so(:,11);
    coalemo0=coalemo;
    coalomo0=coalomo;
    coallmo0=coallmo;
    coalpfmo0=coalpfmo;
    coalomob(:,i)=coal_d_so(:,8);
    coal_r(coal_r>sum(coalomob(:,i)))=sum(coalomob(:,i));% reduction cannot beyond inventory
    for j=1:l
        if sum(coalomob(1:j,i))>=coal_r(i) && sum(coalomob(1:(j-1),i))<coal_r(i)
            coalemo(1:(j-1),i)=0;
            coalomo(1:(j-1),i)=0;
            coallmo(1:(j-1),i)=0;
            coalpfmo(1:(j-1),i)=0;
            coalomo(j,i)=sum(coalomob(1:j,i))-coal_r(i);
            coalemo(j,i)=coalemo0(j,i)./coalomob(j,i).*(sum(coalomob(1:j,i))-coal_r(i));
            coallmo(j,i)=coallmo0(j,i)./coalomob(j,i).*(sum(coalomob(1:j,i))-coal_r(i));
            coalpfmo(j,i)=coalpfmo0(j,i)./coalomob(j,i).*(sum(coalomob(1:j,i))-coal_r(i));           
        end
    end
    %for carbon tax
    for j=1:l
        if coalemo(j,i)==0 && coal_d_so(j,7)==1
            o_fo_o(j,i)=coalomob(j,i);
            e_fo_o(j,i)=coalemo0(j,i);
        end
        if coalemo(j,i)~=0 && coalemo(j,i)<coalemo0(j,i)
            o_fo_o(j,i)=(coal_r(i)-sum(coalomob(1:j-1,i)));
            e_fo_o(j,i)=coalemo0(j,i)./coalomob(j,i).*o_fo_o(j,i);
        end
    end
    %
end
coal_so_ye=sum(coalemo,1)/1000;% yearly emissionGtCO2e
coal_so_ye_cum=cumsum(coal_so_ye);% cumulative yearly emission
coal_so_yo=sum(coalomo,1);% yearly output Mt





% sorted by emission intensity
%for carbon tax
o_nc_e=zeros(l,lper);% output reduction by natural closure
o_pf_e=zeros(l,lper);% output reduction by carbon tax-driven closure
o_fo_e=zeros(l,lper);% output reduction by directly forced closure
e_nc_e=zeros(l,lper);% emission reduction by natural closure
e_pf_e=zeros(l,lper);% emission reduction by carbon tax-driven closure
e_fo_e=zeros(l,lper);% emission reduction by directly forced closure
op=xlsread(path,'23','j2:j1266');%original cost without carbon tax
cart=xlsread(path1,'carbontax','j5:k31');%carbon tax 2030-2050
c=xlsread(path,'23','k2:k1266');%original cost without carbon tax and employment
pf=zeros(l,lper);
spf=zeros(l,lper);% profit with carbon tax for all emission
tpf=zeros(l,lper);
%
coaleme=zeros(l,lper);
coalome=zeros(l,lper);
coallme=zeros(l,lper);
coalpfme=zeros(l,lper);
coalomeb=zeros(l,lper);
late=zeros(l,lper);
lone=zeros(l,lper);
procodee=zeros(l,lper);
coal_d=[coal_d0,t,coal_d0(:,2),latlon,procode];
for i=1:lper % 24-50
    %for carbon tax
    pf(:,i)=op-pei.*cart(i,2);
    spf(:,i)=coal_d(:,5).*cart(i,2);%social cost considering all emissions
    tpf(:,i)=spf(:,i).*coal_d(:,2);% total profit M CNY
    coal_d=[coal_d2, tpf(:,i)];% 12 col
    %
    for j=1:l
        if coal_d(j,6)<=i% no resources
            %for carbon tax
            o_nc_e(j,i)=coal_d2(j,2);
            e_nc_e(j,i)=coal_d2(j,1);
            %
            coal_d(j,7)=0;
            coal_d(j,1)=0;
            coal_d(j,2)=0;
            coal_d(j,3)=0;
            coal_d(j,12)=0;           
        end
        %for carbon tax
        if pf(j,i)<=0 && coal_d(j,6)>i% stranded by carbon tax
            o_pf_e(j,i)=coal_d2(j,2);
            e_pf_e(j,i)=coal_d2(j,1);
            coal_d(j,7)=0;
            coal_d(j,1)=0;
            coal_d(j,2)=0;
            coal_d(j,3)=0;  
            coal_d(j,12)=0;           
        end
        %
    end   
    coal_d_se=sortrows(coal_d,[7 -5 -4]);% sorted by status, emission intensity and depth
    coaleme(:,i)=coal_d_se(:,1);% emission matrix sorted by status, output and depth
    coalome(:,i)=coal_d_se(:,2);% output matrix sorted by status, output and depth
    coallme(:,i)=coal_d_se(:,3);% labor matrix sorted by status, output and depth
    coalpfme(:,i)=coal_d_se(:,12);% profit matrix sorted by status, output and depth
    late(:,i)=coal_d_se(:,9);
    lone(:,i)=coal_d_se(:,10);
    procodee(:,i)=coal_d_se(:,11);
    coaleme0=coaleme;
    coalome0=coalome;
    coallme0=coallme;    
    coalpfme0=coalpfme;
    coalomeb(:,i)=coal_d_se(:,8);
    coal_r(coal_r>sum(coalomeb(:,i)))=sum(coalomeb(:,i));% reduction cannot beyond inventory
    for j=1:l
        if sum(coalomeb(1:j,i))>=coal_r(i) && sum(coalomeb(1:(j-1),i))<coal_r(i)
            coaleme(1:(j-1),i)=0;
            coalome(1:(j-1),i)=0;
            coallme(1:(j-1),i)=0;     
            coalpfme(1:(j-1),i)=0;
            coalome(j,i)=sum(coalomeb(1:j,i))-coal_r(i);
            coaleme(j,i)=coaleme0(j,i)./coalomeb(j,i).*(sum(coalomeb(1:j,i))-coal_r(i));
            coallme(j,i)=coallme0(j,i)./coalomeb(j,i).*(sum(coalomeb(1:j,i))-coal_r(i));  
            coalpfme(j,i)=coalpfme0(j,i)./coalomeb(j,i).*(sum(coalomeb(1:j,i))-coal_r(i));           
        end
    end
    %for carbon tax
    for j=1:l
        if coaleme(j,i)==0 && coal_d_se(j,7)==1
            o_fo_e(j,i)=coalomeb(j,i);
            e_fo_e(j,i)=coaleme0(j,i);
        end
        if coaleme(j,i)~=0 && coaleme(j,i)<coaleme0(j,i)
            o_fo_e(j,i)=(coal_r(i)-sum(coalomeb(1:j-1,i)));
            e_fo_e(j,i)=coaleme0(j,i)./coalomeb(j,i).*o_fo_e(j,i);
        end
    end
    %
end
coal_se_ye=sum(coaleme,1)/1000;% yearly emission GtCO2e
coal_se_ye_cum=cumsum(coal_se_ye);% cumulative yearly emission
coal_se_yo=sum(coalome,1);% yearly output Mt


plot(x,coal_so_ye,'-','color',[6, 107, 82]./255,'Markersize',6,'LineWidth',5)
hold on
plot(x,coal_se_ye,'-o','color',[6, 107, 82]./255,'Markersize',6,'LineWidth',5)
ylabel('Coal (GtCO_2e)','FontSize',50);
xlim([2024 2050]);
ylim([0 10]);
xticks([2025,2030,2035,2040,2045,2050]);
lgd = legend('STEPS-small','STEPS-high','APS-small','APS-high','Location','northeast','Orientation','vertical');
set(gca,'FontSize',100);
diffoe_cs=coal_so_ye_cum(1,size(coal_so_ye_cum,2))-coal_se_ye_cum(1,size(coal_se_ye_cum,2));
% output
xlswrite(pathtca,coalomo,'coalomo','a1');
xlswrite(pathtca,coalemo,'coalemo','a1');
xlswrite(pathtca,coallmo,'coallmo','a1');
xlswrite(pathtca,coalpfmo,'coalpfmo','a1');
xlswrite(pathtca,coalome,'coalome','a1');
xlswrite(pathtca,coaleme,'coaleme','a1');
xlswrite(pathtca,coallme,'coallme','a1');
xlswrite(pathtca,coalpfme,'coalpfme','a1');
xlswrite(pathtca,coal_so_yo,'coal_so_yo','a1');
xlswrite(pathtca,coal_so_ye,'coal_so_ye','a1');
xlswrite(pathtca,coal_se_yo,'coal_se_yo','a1');
xlswrite(pathtca,coal_se_ye,'coal_se_ye','a1');
xlswrite(pathtca,o_nc_o,'o_nc_o','a1');
xlswrite(pathtca,o_pf_o,'o_pf_o','a1');
xlswrite(pathtca,o_fo_o,'o_fo_o','a1');
xlswrite(pathtca,o_nc_e,'o_nc_e','a1');
xlswrite(pathtca,o_pf_e,'o_pf_e','a1');
xlswrite(pathtca,o_fo_e,'o_fo_e','a1');

xlswrite(pathtca,e_nc_o,'e_nc_o','a1');
xlswrite(pathtca,e_pf_o,'e_pf_o','a1');
xlswrite(pathtca,e_fo_o,'e_fo_o','a1');
xlswrite(pathtca,e_nc_e,'e_nc_e','a1');
xlswrite(pathtca,e_pf_e,'e_pf_e','a1');
xlswrite(pathtca,e_fo_e,'e_fo_e','a1');

xlswrite(pathtca,lato,'lato','a1');
xlswrite(pathtca,lono,'lono','a1');
xlswrite(pathtca,procodeo,'procodeo','a1');
xlswrite(pathtca,late,'late','a1');
xlswrite(pathtca,lone,'lone','a1');
xlswrite(pathtca,procodee,'procodee','a1');

xlswrite(pathtca,coal_r','overview','b11');

pceo=[zeros(31,27),(1:1:31)'];
for i=1:27
    for j=1:l
        for g=1:31
            if procodeo(j,i)==pceo(g,28)
                pceo(g,i)=pceo(g,i)+coalemo(j,i);
            end
        end
    end
end
tpceo=sum(pceo(:,1:27),2);
xlswrite(pathtca,tpceo,'tpceo','a1');

pcee=[zeros(31,27),(1:1:31)'];
for i=1:27
    for j=1:l
        for g=1:31
            if procodee(j,i)==pcee(g,28)
                pcee(g,i)=pcee(g,i)+coaleme(j,i);
            end
        end
    end
end
tpcee=sum(pcee(:,1:27),2);
xlswrite(pathtca,tpcee,'tpcee','a1');

%composition coal_so_yo
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(o_nc_o);sum(o_pf_o);sum(o_fo_o)]';
area(x,y);
ylabel0=['APS' newline 'output reduction-Coal (Bt)'];
ylabel(ylabel0)
xlim([2024 2050]);
legend({'Natural','Carbon tax','Forced'},'Location','northwest')
set(gca,'FontSize',50);
print('F:\paper\science\figure\output reduction composition tax aps coal.tiff', '-dtiff', ['-r' num2str(dpi)]);

%composition coal_so_ye
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(e_nc_o);sum(e_pf_o);sum(e_fo_o)]'./1000;
area(x,y);
legend({'Natural','Carbon tax','Forced'},'Location','northwest')
ylabel0=['APS-PSO' newline 'emission reduction-Coal (Mt)'];
ylabel(ylabel0)
xlim([2024 2050]);
set(gca,'FontSize',50);
print('F:\paper\science\figure\emission reduction composition tax aps Scale coal.tiff', '-dtiff', ['-r' num2str(dpi)]);

%composition coal_se_ye
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(e_nc_e);sum(e_pf_e);sum(e_fo_e)]'./1000;
area(x,y);
legend({'Natural','Carbon tax','Forced'},'Location','northwest')
ylabel0=['APS-PHEI' newline 'emission reduction-Coal (Mt)'];
ylabel(ylabel0)
xlim([2024 2050]);
set(gca,'FontSize',50);
print('F:\paper\science\figure\emission reduction composition tax aps Emission coal.tiff', '-dtiff', ['-r' num2str(dpi)]);


%emission intensity tco2e/t
eiao=round(coal_so_ye,6).*1000./round(coal_so_yo,6);
eiae=round(coal_se_ye,6).*1000./round(coal_se_yo,6);

f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
x=2024:1:2050;
plot(x,eiso,'-','color',[235, 23, 23]./255,'Markersize',6,'LineWidth',5)
hold on
plot(x,eise,'-o','color',[235, 23, 23]./255,'Markersize',6,'LineWidth',5)
hold on
plot(x,eiao,'-','color',[6, 107, 82]./255,'Markersize',6,'LineWidth',5)
hold on
plot(x,eiae,'-o','color',[6, 107, 82]./255,'Markersize',6,'LineWidth',5)
ylabel('Emission intensity-Coal-tax (tCO_2e/t)','FontSize',45);
xlim([2024 2050]);
ylim([1.5 2.5]);
xticks([2025,2030,2035,2040,2045,2050]);
lgd = legend('STEPS-PSO','STEPS-PHEI','APS-PSO','APS-PHEI','Location','southwest','Orientation','vertical');
set(gca,'FontSize',45);


%% province cumulative emissions
%steps-small high
path=('F:\paper\science\data\coal.xlsx');
procum0=xlsread(path,'provincecode','aj2:am24');% Mt
[~,province,~]=xlsread(path,'provincecode','ap2:ap24');
procum=procum0./1000;%Gt
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
val=[procum(:,1),procum(:,2)];
sorted_val=sortrows(val,1);
bar(sorted_val);
xticks(1:23);
xtickangle(-90);
xticklabels(province);
ylabel(['Cumulative emissions' newline '(GtCO_2e)'],'FontSize',60);
set(gca,'FontSize',70);
lgd = legend('STEPS-small','STEPS-high','Location','northwest','Orientation','vertical');
print('F:\paper\science\figure\provincial cumulative emissions steps.tiff', '-dtiff', ['-r' num2str(dpi)]);



%aps-small high
path=('F:\paper\science\data\coal.xlsx');
procum0=xlsread(path,'provincecode','aj2:am24');% Mt
[~,province,~]=xlsread(path,'provincecode','ap2:ap24');
procum=procum0./1000;%Gt
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
val=[procum(:,3),procum(:,4)];
sorted_val=sortrows(val,1);
bar(sorted_val);
xticks(1:23);
xtickangle(-90);
xticklabels(province);
ylabel(['Cumulative emissions' newline '(GtCO_2e)'],'FontSize',60);
set(gca,'FontSize',70);
lgd = legend('APS-small','APS-high','Location','northwest','Orientation','vertical');
print('F:\paper\science\figure\provincial cumulative emissions aps.tiff', '-dtiff', ['-r' num2str(dpi)]);

%% output and emission intensity of coal
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
plot(coal_d0(:,2),coal_d0(:,5),'ok','Markersize',15,'LineWidth',2)
xlim([0 35]);
ylim([1.5 3.5]);
set(gca,'FontSize',110);
xlabel('Output of coal mine (t)','FontSize',110);
ylabel('tCO_2e/t','FontSize',110);
print('F:\paper\science\figure\coal output emission intensity.tiff', '-dtiff', ['-r' num2str(dpi)]);


close all
toc



% tax
%% oil Stated Policies
clear;clc;
tic
path=('F:\paper\science\data\oil.xlsx');
path1=('F:\paper\science\data\futureproduction.xlsx');
pathtcs=('F:\paper\science\results\tax\oil\STEPS.xlsx');
oil_d0=xlsread(path,'23','d2:i842');% emission/yr,oil output/yr, and workforce of each oil mine
latlon=xlsread(path,'23','b2:c842');
procode=xlsread(path,'23','jr2:jr842');
procode(isnan(procode))=32;
l=size(oil_d0,1);
t=ones(l,1);% 1 for sites with output
oil_d=[oil_d0,t,oil_d0(:,2),latlon,procode];
oil_d2=oil_d;
oil_d3=oil_d;
oil_ew=xlsread(path1,'oil','k5:k31');% future oil production
oil_ew2=oil_ew-24.29;% exclude oil that are not covered by inventory
oil_r=-(oil_ew2-sum(oil_d(:,2)));%reduction in inventory needed
per=2024:1:2050;% period span
lper=size(per,2);
%for carbon tax
pei=xlsread(path,'23','m2:m842');%emission inensity of production tCO2e/t
%
% sorted by output
%for carbon tax
o_nc_o=zeros(l,lper);% output reduction by natural closure
o_pf_o=zeros(l,lper);% output reduction by carbon tax-driven closure
o_fo_o=zeros(l,lper);% output reduction by directly forced closure
e_nc_o=zeros(l,lper);% emission reduction by natural closure
e_pf_o=zeros(l,lper);% emission reduction by carbon tax-driven closure
e_fo_o=zeros(l,lper);% emission reduction by directly forced closure
op=xlsread(path,'23','j2:j842');%original profit without carbon tax
cart=xlsread(path1,'carbontax','j5:k31');%carbon tax 2030-2050
c=xlsread(path,'23','k2:k842');%original cost without carbon tax and employment
pf=zeros(l,lper);% output reduction by directly forced closure
%
spf=zeros(l,lper);% profit with carbon tax for all emission
tpf=zeros(l,lper);
oilemo=zeros(l,lper);
oilomo=zeros(l,lper);
oillmo=zeros(l,lper);
oilpfmo=zeros(l,lper);
oilomob=zeros(l,lper);
lato=zeros(l,lper);
lono=zeros(l,lper);
procodeo=zeros(l,lper);

for i=1:lper % 24-50
    %for carbon tax
    pf(:,i)=op-pei.*cart(i,1);
    spf(:,i)=c+oil_d(:,5).*cart(i,1);% cost considering all emissions
    tpf(:,i)=spf(:,i).*oil_d(:,2);% total cost M CNY
    oil_d=[oil_d2, tpf(:,i)];% 12 col
    for j=1:l
        if oil_d(j,6)<=i% no resources
            %for carbon tax
            o_nc_o(j,i)=oil_d2(j,2);
            e_nc_o(j,i)=oil_d2(j,1);
            %
            oil_d(j,7)=0;
            oil_d(j,1)=0;
            oil_d(j,2)=0;
            oil_d(j,3)=0;
            oil_d(j,12)=0;          
        end
        %for carbon tax
        if pf(j,i)<=0 && oil_d(j,6)>i% stranded by carbon tax
            o_pf_o(j,i)=oil_d2(j,2);
            e_pf_o(j,i)=oil_d2(j,1);
            oil_d(j,7)=0;
            oil_d(j,1)=0;
            oil_d(j,2)=0;
            oil_d(j,3)=0;
            oil_d(j,12)=0;                    
        end
        %
    end
    oil_d_so=sortrows(oil_d,[7 8 -4]);% sorted by status, output and cost
    oilemo(:,i)=oil_d_so(:,1);% emission matrix sorted by status, output and cost
    oilomo(:,i)=oil_d_so(:,2);% output matrix sorted by status, output and cost
    oillmo(:,i)=oil_d_so(:,3);% output matrix sorted by status, output and cost
    oilpfmo(:,i)=oil_d_so(:,12);% profit matrix sorted by status, output and depth
    lato(:,i)=oil_d_so(:,9);
    lono(:,i)=oil_d_so(:,10);
    procodeo(:,i)=oil_d_so(:,11);
    oilemo0=oilemo;
    oilomo0=oilomo;
    oillmo0=oillmo;
    oilpfmo0=oilpfmo;    
    oilomob(:,i)=oil_d_so(:,8);
    oil_r(oil_r>sum(oilomob(:,i)))=sum(oilomob(:,i));% reduction cannot beyond inventory
    for j=1:l
        if sum(oilomob(1:j,i))>=oil_r(i) && sum(oilomob(1:(j-1),i))<oil_r(i) && oilomo0(j,i)==1
            oilemo(1:(j-1),i)=0;
            oilomo(1:(j-1),i)=0;
            oillmo(1:(j-1),i)=0;
            oilpfmo(1:(j-1),i)=0;            
            oilomo(j,i)=sum(oilomob(1:j,i))-oil_r(i);
            oilemo(j,i)=oilemo0(j,i)./oilomob(j,i).*(sum(oilomob(1:j,i))-oil_r(i));
            oillmo(j,i)=oillmo0(j,i)./oilomob(j,i).*(sum(oilomob(1:j,i))-oil_r(i));       
            oilpfmo(j,i)=oilpfmo0(j,i)./oilomob(j,i).*(sum(oilomob(1:j,i))-oil_r(i));                  
        end
    end
    %for carbon tax
    for j=1:l
        if oilemo(j,i)==0 && oil_d_so(j,7)==1
            o_fo_o(j,i)=oilomob(j,i);
            e_fo_o(j,i)=oilemo0(j,i);
        end
        if oilemo(j,i)~=0 && oilemo(j,i)<oilemo0(j,i)
            o_fo_o(j,i)=(oil_r(i)-sum(oilomob(1:j-1,i)));
            e_fo_o(j,i)=oilemo0(j,i)./oilomob(j,i).*o_fo_o(j,i);
        end
    end
    %
end
oil_so_ye=sum(oilemo,1)/1000;% yearly emissionGtCO2e
oil_so_ye_cum=cumsum(oil_so_ye);% cumulative yearly emission
oil_so_yo=sum(oilomo,1);% yearly output Mt





% sorted by emission intensity
%for carbon tax
o_nc_e=zeros(l,lper);% output reduction by natural closure
o_pf_e=zeros(l,lper);% output reduction by carbon tax-driven closure
o_fo_e=zeros(l,lper);% output reduction by directly forced closure
e_nc_e=zeros(l,lper);% emission reduction by natural closure
e_pf_e=zeros(l,lper);% emission reduction by carbon tax-driven closure
e_fo_e=zeros(l,lper);% emission reduction by directly forced closure
op=xlsread(path,'23','j2:j842');%original cost without carbon tax
cart=xlsread(path1,'carbontax','j5:k31');%carbon tax 2030-2050
c=xlsread(path,'23','k2:k842');%original cost without carbon tax and employment
pf=zeros(l,lper);% output reduction by directly forced closure
spf=zeros(l,lper);% profit with carbon tax for all emission
tpf=zeros(l,lper);
%
oileme=zeros(l,lper);
oilome=zeros(l,lper);
oillme=zeros(l,lper);
oilpfme=zeros(l,lper);
oilomeb=zeros(l,lper);
late=zeros(l,lper);
lone=zeros(l,lper);
procodee=zeros(l,lper);
oil_d=[oil_d0,t,oil_d0(:,2),latlon,procode];
for i=1:lper % 24-50
    %for carbon tax
    pf(:,i)=op-pei.*cart(i,1);
    spf(:,i)=c+oil_d(:,5).*cart(i,1);%social cost considering all emissions
    tpf(:,i)=spf(:,i).*oil_d(:,2);% total profit M CNY
    oil_d=[oil_d2, tpf(:,i)];% 12 col
    for j=1:l
        if oil_d(j,6)<=i% no resources
            %for carbon tax
            o_nc_e(j,i)=oil_d2(j,2);
            e_nc_e(j,i)=oil_d2(j,1);
            %
            oil_d(j,7)=0;
            oil_d(j,1)=0;
            oil_d(j,2)=0;
            oil_d(j,3)=0;   
            oil_d(j,12)=0;              
        end
        %for carbon tax
        if pf(j,i)<=0 && oil_d(j,6)>i% stranded by carbon tax
            o_pf_e(j,i)=oil_d2(j,2);
            e_pf_e(j,i)=oil_d2(j,1);
            oil_d(j,7)=0;
            oil_d(j,1)=0;
            oil_d(j,2)=0;
            oil_d(j,3)=0;  
            oil_d(j,12)=0;                          
        end
        %
    end
    
    oil_d_se=sortrows(oil_d,[7 -5 -4]);% sorted by status, emission intensity and cost
    oileme(:,i)=oil_d_se(:,1);% emission matrix sorted by status, output and cost
    oilome(:,i)=oil_d_se(:,2);% output matrix sorted by status, output and cost
    oillme(:,i)=oil_d_se(:,3);% output matrix sorted by status, output and cost 
    oilpfme(:,i)=oil_d_se(:,12);% profit matrix sorted by status, output and depth
    late(:,i)=oil_d_se(:,9);
    lone(:,i)=oil_d_se(:,10);
    procodee(:,i)=oil_d_se(:,11);
    oileme0=oileme;
    oilome0=oilome;
    oillme0=oillme; 
    oilpfme0=oilpfme;     
    oilomeb(:,i)=oil_d_se(:,8);
    oil_r(oil_r>sum(oilomeb(:,i)))=sum(oilomeb(:,i));% reduction cannot beyond inventory
    for j=1:l
        if sum(oilomeb(1:j,i))>=oil_r(i) && sum(oilomeb(1:(j-1),i))<oil_r(i) && oilome0(j,i)==1
            oileme(1:(j-1),i)=0;
            oilome(1:(j-1),i)=0;
            oillme(1:(j-1),i)=0;  
            oilpfme(1:(j-1),i)=0;              
            oilome(j,i)=sum(oilomeb(1:j,i))-oil_r(i);
            oileme(j,i)=oileme0(j,i)./oilomeb(j,i).*(sum(oilomeb(1:j,i))-oil_r(i));
            oillme(j,i)=oillme0(j,i)./oilomeb(j,i).*(sum(oilomeb(1:j,i))-oil_r(i));  
            oilpfme(j,i)=oilpfme0(j,i)./oilomeb(j,i).*(sum(oilomeb(1:j,i))-oil_r(i));              
        end
    end
    %for carbon tax
    for j=1:l
        if oileme(j,i)==0 && oil_d_se(j,7)==1
            o_fo_e(j,i)=oilomeb(j,i);
            e_fo_e(j,i)=oileme0(j,i);
        end
        if oileme(j,i)~=0 && oileme(j,i)<oileme0(j,i)
            o_fo_e(j,i)=(oil_r(i)-sum(oilomeb(1:j-1,i)));
            e_fo_e(j,i)=oileme0(j,i)./oilomeb(j,i).*o_fo_e(j,i);
        end
    end
    %
end
oil_se_ye=sum(oileme,1)/1000;% yearly emission GtCO2e
oil_se_ye_cum=cumsum(oil_se_ye);% cumulative yearly emission
oil_se_yo=sum(oilome,1);% yearly output Mt

f=figure;
x=2024:1:2050;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi=300;
plot(x,oil_so_ye,'-','color',[235, 23, 23]./255,'Markersize',6,'LineWidth',5)
hold on
plot(x,oil_se_ye,'-','color',[235, 23, 23]./255,'Markersize',6,'LineWidth',5)
ylabel('Emissions-oil (GtCO_2e)','FontSize',50);
xlim([2024 2050]);
ylim([0 1]);
xticks([2025,2030,2035,2040,2045,2050]);
set(gca,'FontSize',80);
print('F:\paper\science\figure\Emissions tax steps oil.tiff', '-dtiff', ['-r' num2str(dpi)]);
diffoe_cs=oil_so_ye_cum(1,size(oil_so_ye_cum,2))-oil_se_ye_cum(1,size(oil_se_ye_cum,2));
% output
xlswrite(pathtcs,oilomo,'oilomo','a1');
xlswrite(pathtcs,oilemo,'oilemo','a1');
xlswrite(pathtcs,oillmo,'oillmo','a1');
xlswrite(pathtcs,oilpfmo,'oilpfmo','a1');
xlswrite(pathtcs,oilome,'oilome','a1');
xlswrite(pathtcs,oileme,'oileme','a1');
xlswrite(pathtcs,oillme,'oillme','a1');
xlswrite(pathtcs,oilpfme,'oilpfme','a1');
xlswrite(pathtcs,oil_so_ye,'oil_so_ye','a1');
xlswrite(pathtcs,oil_se_ye,'oil_se_ye','a1');
xlswrite(pathtcs,o_nc_o,'o_nc_o','a1');
xlswrite(pathtcs,o_pf_o,'o_pf_o','a1');
xlswrite(pathtcs,o_fo_o,'o_fo_o','a1');
xlswrite(pathtcs,o_nc_e,'o_nc_e','a1');
xlswrite(pathtcs,o_pf_e,'o_pf_e','a1');
xlswrite(pathtcs,o_fo_e,'o_fo_e','a1');

xlswrite(pathtcs,e_nc_o,'e_nc_o','a1');
xlswrite(pathtcs,e_pf_o,'e_pf_o','a1');
xlswrite(pathtcs,e_fo_o,'e_fo_o','a1');
xlswrite(pathtcs,e_nc_e,'e_nc_e','a1');
xlswrite(pathtcs,e_pf_e,'e_pf_e','a1');
xlswrite(pathtcs,e_fo_e,'e_fo_e','a1');

xlswrite(pathtcs,oil_r','overview','a11');

xlswrite(pathtcs,oil_so_yo,'oil_so_yo','a1');
xlswrite(pathtcs,oil_se_yo,'oil_se_yo','a1');

xlswrite(pathtcs,lato,'lato','a1');
xlswrite(pathtcs,lono,'lono','a1');
xlswrite(pathtcs,procodeo,'procodeo','a1');
xlswrite(pathtcs,late,'late','a1');
xlswrite(pathtcs,lone,'lone','a1');
xlswrite(pathtcs,procodee,'procodee','a1');

pceo=[zeros(31,27),(1:1:31)'];
for i=1:27
    for j=1:l
        for g=1:31
            if procodeo(j,i)==pceo(g,28)
                pceo(g,i)=pceo(g,i)+oilemo(j,i);
            end
        end
    end
end
tpceo=sum(pceo(:,1:27),2);
xlswrite(pathtcs,tpceo,'tpceo','a1');

pcee=[zeros(31,27),(1:1:31)'];
for i=1:27
    for j=1:l
        for g=1:31
            if procodee(j,i)==pcee(g,28)
                pcee(g,i)=pcee(g,i)+oileme(j,i);
            end
        end
    end
end
tpcee=sum(pcee(:,1:27),2);
xlswrite(pathtcs,tpcee,'tpcee','a1');

%composition oil_so_yo
f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(o_nc_o);sum(o_fo_o)]';
area(x,y);
ylabel0=['STEPS' newline 'output reduction-oil-tax (Bt)'];
ylabel(ylabel0)
xlim([2024 2050]);
legend({'Natural'},'Location','northwest')
set(gca,'FontSize',50);
print('F:\paper\science\figure\output reduction composition tax steps oil.tiff', '-dtiff', ['-r' num2str(dpi)]);

%composition oil_so_ye
f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(e_nc_o);sum(e_fo_o)]'./1000;
area(x,y);
legend({'Natural'},'Location','northwest')
ylabel0=['STEPS' newline 'emission reduction-oil-tax (Mt)'];
ylabel(ylabel0)
xlim([2024 2050]);
set(gca,'FontSize',50);
print('F:\paper\science\figure\emission reduction composition tax steps oil.tiff', '-dtiff', ['-r' num2str(dpi)]);


%emission intensity tco2e/t
eiso=round(oil_so_ye,6).*1000./round(oil_so_yo,6);
eise=round(oil_se_ye,6).*1000./round(oil_se_yo,6);




%% oil Announced Pledges
path=('F:\paper\science\data\oil.xlsx');
path1=('F:\paper\science\data\futureproduction.xlsx');
pathtca=('F:\paper\science\results\tax\oil\APS.xlsx');
oil_d0=xlsread(path,'23','d2:i842');% emission/yr,oil output/yr, and workforce of each oil mine
l=size(oil_d0,1);
t=ones(l,1);% 1 for sites with output
oil_d=[oil_d0,t,oil_d0(:,2),latlon,procode];
oil_d2=oil_d;
oil_d3=oil_d;
oil_ew=xlsread(path1,'oil','l5:l31');% future oil production
oil_ew2=oil_ew-24.29;% exclude oil that are not covered by inventory
oil_r=-(oil_ew2-sum(oil_d(:,2)));%reduction in inventory needed
per=2024:1:2050;% period span
lper=size(per,2);
%for carbon tax
pei=xlsread(path,'23','m2:m842');%emission inensity of production tCO2e/t
%
% sorted by output
%for carbon tax
o_nc_o=zeros(l,lper);% output reduction by natural closure
o_pf_o=zeros(l,lper);% output reduction by carbon tax-driven closure
o_fo_o=zeros(l,lper);% output reduction by directly forced closure
e_nc_o=zeros(l,lper);% emission reduction by natural closure
e_pf_o=zeros(l,lper);% emission reduction by carbon tax-driven closure
e_fo_o=zeros(l,lper);% emission reduction by directly forced closure
op=xlsread(path,'23','j2:j842');%original profit without carbon tax
cart=xlsread(path1,'carbontax','j5:k31');%carbon tax 2030-2050
c=xlsread(path,'23','k2:k842');%original cost without carbon tax and employment
pf=zeros(l,lper);% output reduction by directly forced closure
spf=zeros(l,lper);% profit with carbon tax for all emission
tpf=zeros(l,lper);
oilemo=zeros(l,lper);
oilomo=zeros(l,lper);
oillmo=zeros(l,lper);
oilpfmo=zeros(l,lper);
oilomob=zeros(l,lper);
lato=zeros(l,lper);
lono=zeros(l,lper);
procodeo=zeros(l,lper);

for i=1:lper % 24-50
    %for carbon tax
    pf(:,i)=op-pei.*cart(i,2);
    spf(:,i)=c+oil_d(:,5).*cart(i,2);%social cost considering all emissions
    tpf(:,i)=spf(:,i).*oil_d(:,2);% total profit M CNY
    oil_d=[oil_d2, tpf(:,i)];% 12 col
    for j=1:l
        if oil_d(j,6)<=i% no resources
            %for carbon tax
            o_nc_o(j,i)=oil_d2(j,2);
            e_nc_o(j,i)=oil_d2(j,1);
            %
            oil_d(j,7)=0;
            oil_d(j,1)=0;
            oil_d(j,2)=0;
            oil_d(j,3)=0;  
            oil_d(j,12)=0;             
        end
        %for carbon tax
        if pf(j,i)<=0 && oil_d(j,6)>i% stranded by carbon tax
            o_pf_o(j,i)=oil_d2(j,2);
            e_pf_o(j,i)=oil_d2(j,1);
            oil_d(j,7)=0;
            oil_d(j,1)=0;
            oil_d(j,2)=0;
            oil_d(j,3)=0;     
            oil_d(j,12)=0;                        
        end
        %
    end
    oil_d_so=sortrows(oil_d,[7 8 -4]);% sorted by status, output and cost
    oilemo(:,i)=oil_d_so(:,1);% emission matrix sorted by status, output and cost
    oilomo(:,i)=oil_d_so(:,2);% output matrix sorted by status, output and cost
    oillmo(:,i)=oil_d_so(:,3);% output matrix sorted by status, output and cost    
    oilpfmo(:,i)=oil_d_so(:,12);% output matrix sorted by status, output and cost        
    lato(:,i)=oil_d_so(:,9);
    lono(:,i)=oil_d_so(:,10);
    procodeo(:,i)=oil_d_so(:,11);
    oilemo0=oilemo;
    oilomo0=oilomo;
    oillmo0=oillmo; 
    oilpfmo0=oilpfmo;     
    oilomob(:,i)=oil_d_so(:,8);
    oil_r(oil_r>sum(oilomob(:,i)))=sum(oilomob(:,i));% reduction cannot beyond inventory
    for j=1:l
        if sum(oilomob(1:j,i))>=oil_r(i) && sum(oilomob(1:(j-1),i))<oil_r(i) && oilomo0(j,i)==1
            oilemo(1:(j-1),i)=0;
            oilomo(1:(j-1),i)=0;
            oillmo(1:(j-1),i)=0;  
            oilpfmo(1:(j-1),i)=0;              
            oilomo(j,i)=sum(oilomob(1:j,i))-oil_r(i);
            oilemo(j,i)=oilemo0(j,i)./oilomob(j,i).*(sum(oilomob(1:j,i))-oil_r(i));
            oilpfmo(j,i)=oilpfmo0(j,i)./oilomob(j,i).*(sum(oilomob(1:j,i))-oil_r(i));              
        end
    end
    %for carbon tax
    for j=1:l
        if oilemo(j,i)==0 && oil_d_so(j,7)==1
            o_fo_o(j,i)=oilomob(j,i);
            e_fo_o(j,i)=oilemo0(j,i);
        end
        if oilemo(j,i)~=0 && oilemo(j,i)<oilemo0(j,i)
            o_fo_o(j,i)=(oil_r(i)-sum(oilomob(1:j-1,i)));
            e_fo_o(j,i)=oilemo0(j,i)./oilomob(j,i).*o_fo_o(j,i);
        end
    end
    %
end
oil_so_ye=sum(oilemo,1)/1000;% yearly emissionGtCO2e
oil_so_ye_cum=cumsum(oil_so_ye);% cumulative yearly emission
oil_so_yo=sum(oilomo,1);% yearly output Mt





% sorted by emission intensity
%for carbon tax
o_nc_e=zeros(l,lper);% output reduction by natural closure
o_pf_e=zeros(l,lper);% output reduction by carbon tax-driven closure
o_fo_e=zeros(l,lper);% output reduction by directly forced closure
e_nc_e=zeros(l,lper);% emission reduction by natural closure
e_pf_e=zeros(l,lper);% emission reduction by carbon tax-driven closure
e_fo_e=zeros(l,lper);% emission reduction by directly forced closure
op=xlsread(path,'23','j2:j842');%original cost without carbon tax
cart=xlsread(path1,'carbontax','j5:k31');%carbon tax 2030-2050
c=xlsread(path,'23','k2:k842');%original cost without carbon tax and employment
pf=zeros(l,lper);% output reduction by directly forced closure
spf=zeros(l,lper);% profit with carbon tax for all emission
tpf=zeros(l,lper);
%
oileme=zeros(l,lper);
oilome=zeros(l,lper);
oillme=zeros(l,lper);
oilpfme=zeros(l,lper);
oilomeb=zeros(l,lper);
late=zeros(l,lper);
lone=zeros(l,lper);
procodee=zeros(l,lper);
oil_d=[oil_d0,t,oil_d0(:,2),latlon,procode];
for i=1:lper % 24-50
    %for carbon tax
    pf(:,i)=op-pei.*cart(i,2);
    spf(:,i)=c+oil_d(:,5).*cart(i,2);%social cost considering all emissions
    tpf(:,i)=spf(:,i).*oil_d(:,2);% total profit M CNY
    oil_d=[oil_d2, tpf(:,i)];% 12 col
    %
    for j=1:l
        if oil_d(j,6)<=i% no resources
            %for carbon tax
            o_nc_e(j,i)=oil_d2(j,2);
            e_nc_e(j,i)=oil_d2(j,1);
            %
            oil_d(j,7)=0;
            oil_d(j,1)=0;
            oil_d(j,2)=0;
            oil_d(j,3)=0; 
            oil_d(j,12)=0;             
        end
        %for carbon tax
        if pf(j,i)<=0 && oil_d(j,6)>i% stranded by carbon tax
            o_pf_e(j,i)=oil_d2(j,2);
            e_pf_e(j,i)=oil_d2(j,1);
            oil_d(j,7)=0;
            oil_d(j,1)=0;
            oil_d(j,2)=0;
            oil_d(j,3)=0;  
            oil_d(j,12)=0;                        
        end
        %
    end
    
    oil_d_se=sortrows(oil_d,[7 -5 -4]);% sorted by status, emission intensity and cost
    oileme(:,i)=oil_d_se(:,1);% emission matrix sorted by status, output and cost
    oilome(:,i)=oil_d_se(:,2);% output matrix sorted by status, output and cost
    oillme(:,i)=oil_d_se(:,3);% output matrix sorted by status, output and cost 
    oilpfme(:,i)=oil_d_se(:,12);% output matrix sorted by status, output and cost    
    late(:,i)=oil_d_se(:,9);
    lone(:,i)=oil_d_se(:,10);
    procodee(:,i)=oil_d_se(:,11);
    oileme0=oileme;
    oilome0=oilome;
    oillme0=oillme;  
    oilpfme0=oilpfme;      
    oilomeb(:,i)=oil_d_se(:,8);
    oil_r(oil_r>sum(oilomeb(:,i)))=sum(oilomeb(:,i));% reduction cannot beyond inventory
    for j=1:l
        if sum(oilomeb(1:j,i))>=oil_r(i) && sum(oilomeb(1:(j-1),i))<oil_r(i) && oilome0(j,i)==1
            oileme(1:(j-1),i)=0;
            oilome(1:(j-1),i)=0;
            oillme(1:(j-1),i)=0;   
            oilpfme(1:(j-1),i)=0;               
            oilome(j,i)=sum(oilomeb(1:j,i))-oil_r(i);
            oileme(j,i)=oileme0(j,i)./oilomeb(j,i).*(sum(oilomeb(1:j,i))-oil_r(i));
            oillme(j,i)=oillme0(j,i)./oilomeb(j,i).*(sum(oilomeb(1:j,i))-oil_r(i)); 
            oilpfme(j,i)=oilpfme0(j,i)./oilomeb(j,i).*(sum(oilomeb(1:j,i))-oil_r(i));             
        end
    end
    %for carbon tax
    for j=1:l
        if oileme(j,i)==0 && oil_d_se(j,7)==1
            o_fo_e(j,i)=oilomeb(j,i);
            e_fo_e(j,i)=oileme0(j,i);
        end
        if oileme(j,i)~=0 && oileme(j,i)<oileme0(j,i)
            o_fo_e(j,i)=(oil_r(i)-sum(oilomeb(1:j-1,i)));
            e_fo_e(j,i)=oileme0(j,i)./oilomeb(j,i).*o_fo_e(j,i);
        end
    end
    %
end
oil_se_ye=sum(oileme,1)/1000;% yearly emission GtCO2e
oil_se_ye_cum=cumsum(oil_se_ye);% cumulative yearly emission
oil_se_yo=sum(oilome,1);% yearly output Mt

f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
x=2024:1:2050;
plot(x,oil_so_ye,'-','color','k','Markersize',6,'LineWidth',5)
hold on
plot(x,oil_se_ye,'-','color','k','Markersize',6,'LineWidth',5)
ylabel('Oil (GtCO_2e)','FontSize',50);
xlim([2024 2050]);
ylim([0 1]);
xticks([2025,2030,2035,2040,2045,2050]);
set(gca,'FontSize',100);
diffoe_cs=oil_so_ye_cum(1,size(oil_so_ye_cum,2))-oil_se_ye_cum(1,size(oil_se_ye_cum,2));
% output
xlswrite(pathtca,oilomo,'oilomo','a1');
xlswrite(pathtca,oilemo,'oilemo','a1');
xlswrite(pathtca,oillmo,'oillmo','a1');
xlswrite(pathtca,oilpfmo,'oilpfmo','a1');
xlswrite(pathtca,oilome,'oilome','a1');
xlswrite(pathtca,oileme,'oileme','a1');
xlswrite(pathtca,oillme,'oillme','a1');
xlswrite(pathtca,oilpfme,'oilpfme','a1');
xlswrite(pathtca,oil_so_ye,'oil_so_ye','a1');
xlswrite(pathtca,oil_se_ye,'oil_se_ye','a1');
xlswrite(pathtca,o_nc_o,'o_nc_o','a1');
xlswrite(pathtca,o_pf_o,'o_pf_o','a1');
xlswrite(pathtca,o_fo_o,'o_fo_o','a1');
xlswrite(pathtca,o_nc_e,'o_nc_e','a1');
xlswrite(pathtca,o_pf_e,'o_pf_e','a1');
xlswrite(pathtca,o_fo_e,'o_fo_e','a1');

xlswrite(pathtca,e_nc_o,'e_nc_o','a1');
xlswrite(pathtca,e_pf_o,'e_pf_o','a1');
xlswrite(pathtca,e_fo_o,'e_fo_o','a1');
xlswrite(pathtca,e_nc_e,'e_nc_e','a1');
xlswrite(pathtca,e_pf_e,'e_pf_e','a1');
xlswrite(pathtca,e_fo_e,'e_fo_e','a1');

xlswrite(pathtca,oil_r','overview','a11');
xlswrite(pathtca,oil_so_yo,'oil_so_yo','a1');
xlswrite(pathtca,oil_se_yo,'oil_se_yo','a1');

xlswrite(pathtca,lato,'lato','a1');
xlswrite(pathtca,lono,'lono','a1');
xlswrite(pathtca,procodeo,'procodeo','a1');
xlswrite(pathtca,late,'late','a1');
xlswrite(pathtca,lone,'lone','a1');
xlswrite(pathtca,procodee,'procodee','a1');

pceo=[zeros(31,27),(1:1:31)'];
for i=1:27
    for j=1:l
        for g=1:31
            if procodeo(j,i)==pceo(g,28)
                pceo(g,i)=pceo(g,i)+oilemo(j,i);
            end
        end
    end
end
tpceo=sum(pceo(:,1:27),2);
xlswrite(pathtca,tpceo,'tpceo','a1');

pcee=[zeros(31,27),(1:1:31)'];
for i=1:27
    for j=1:l
        for g=1:31
            if procodee(j,i)==pcee(g,28)
                pcee(g,i)=pcee(g,i)+oileme(j,i);
            end
        end
    end
end
tpcee=sum(pcee(:,1:27),2);
xlswrite(pathtca,tpcee,'tpcee','a1');

%composition oil_so_yo
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(o_nc_o);sum(o_pf_o)]';
area(x,y);
ylabel0=['APS' newline 'output reduction-oil-tax (Bt)'];
ylabel(ylabel0)
xlim([2024 2050]);
legend({'Natural','Carbon tax'},'Location','northwest')
set(gca,'FontSize',50);
print('F:\paper\science\figure\output reduction composition tax APS oil.tiff', '-dtiff', ['-r' num2str(dpi)]);

%composition oil_so_ye
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(e_nc_o);sum(e_pf_o)]'./1000;
area(x,y);
legend({'Natural','Carbon tax'},'Location','northwest')
ylabel0=['APS' newline 'emission reduction-oil-tax (Mt)'];
ylabel(ylabel0)
xlim([2024 2050]);
set(gca,'FontSize',50);
print('F:\paper\science\figure\emission reduction composition tax APS Scale oil.tiff', '-dtiff', ['-r' num2str(dpi)]);

%composition oil_se_ye
f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(e_nc_e);sum(e_pf_e)]'./1000;
area(x,y);
legend({'Natural','Carbon tax'},'Location','northwest')
ylabel0=['APS' newline 'emission reduction-oil-tax (Mt)'];
ylabel(ylabel0)
xlim([2024 2050]);
set(gca,'FontSize',50);
print('F:\paper\science\figure\emission reduction composition tax APS Emission oil.tiff', '-dtiff', ['-r' num2str(dpi)]);


%emission intensity tco2e/t
eiao=round(oil_so_ye,6).*1000./round(oil_so_yo,6);
eiae=round(oil_se_ye,6).*1000./round(oil_se_yo,6);

f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
x=2024:1:2050;
plot(x,eiso,'-','color','k','Markersize',6,'LineWidth',5)
%hold on
%plot(x,eise,'-','color',[235, 23, 23]./255,'Markersize',6,'LineWidth',5)
%hold on
%plot(x,eiao,'-','color',[6, 107, 82]./255,'Markersize',6,'LineWidth',5)
%hold on
%plot(x,eiae,'-','color',[6, 107, 82]./255,'Markersize',6,'LineWidth',5)
ylabel('Emission intensity-oil-tax (tCO_2e/t)','FontSize',45);
xlim([2024 2050]);
ylim([3 4]);
xticks([2025,2030,2035,2040,2045,2050]);
%lgd = legend('STEPS-PSO','STEPS-PHEI','APS-PSO','APS-PHEI','Location','southwest','Orientation','vertical');
set(gca,'FontSize',45);
print('F:\paper\science\figure\emission reduction composition tax APS Emission oil.tiff', '-dtiff', ['-r' num2str(dpi)]);

%% output and emission intensity of oil
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
plot(oil_d0(:,2),oil_d0(:,5),'ok','Markersize',15,'LineWidth',2)
ylabel('tCO_2e/t','FontSize',110);
xlabel('Output of oil field (t)','FontSize',110);
xlim([0 15]);
ylim([3 4]);
set(gca,'FontSize',110);
print('F:\paper\science\figure\oil output emission intensity.tiff', '-dtiff', ['-r' num2str(dpi)]);


close all
toc


%tax
%% gas Stated Policies
clear;clc;
tic
path=('F:\paper\science\data\gas.xlsx');
path1=('F:\paper\science\data\futureproduction.xlsx');
pathtcs=('F:\paper\science\results\tax\gas\STEPS.xlsx');
gas_d0=xlsread(path,'23','d2:i588');% emission/yr,gas output/yr, and workforce of each gas mine
latlon=xlsread(path,'23','b2:c588');
procode=xlsread(path,'23','jr2:jr588');
procode(isnan(procode))=32;
l=size(gas_d0,1);
t=ones(l,1);% 1 for sites with output
gas_d=[gas_d0,t,gas_d0(:,2),latlon,procode];
gas_d2=gas_d;
gas_d3=gas_d;
gas_ew=xlsread(path1,'gas','k5:k31');% future gas production
gas_ew2=gas_ew-41.188;% exclude gas that are not covered by inventory
gas_r=-(gas_ew2-sum(gas_d(:,2)));%reduction in inventory needed
per=2024:1:2050;% period span
lper=size(per,2);
%for carbon tax
pei=xlsread(path,'23','m2:m588');%emission inensity of production tCO2e/t
%
% sorted by output
%for carbon tax
o_nc_o=zeros(l,lper);% output reduction by natural closure
o_pf_o=zeros(l,lper);% output reduction by carbon tax-driven closure
o_fo_o=zeros(l,lper);% output reduction by directly forced closure
e_nc_o=zeros(l,lper);% emission reduction by natural closure
e_pf_o=zeros(l,lper);% emission reduction by carbon tax-driven closure
e_fo_o=zeros(l,lper);% emission reduction by directly forced closure
op=xlsread(path,'23','j2:j588');%original profit without carbon tax
cart=xlsread(path1,'carbontax','j5:k31');%carbon tax 2030-2050
c=xlsread(path,'23','k2:k588');%original cost without carbon tax and employment
pf=zeros(l,lper);% output reduction by directly forced closure
spf=zeros(l,lper);% profit with carbon tax for all emission
tpf=zeros(l,lper);
%
gasemo=zeros(l,lper);
gasomo=zeros(l,lper);
gaslmo=zeros(l,lper);
gaspfmo=zeros(l,lper);
gasomob=zeros(l,lper);
lato=zeros(l,lper);
lono=zeros(l,lper);
procodeo=zeros(l,lper);
for i=1:lper % 24-50
    %for carbon tax
    pf(:,i)=op-pei.*cart(i,1);
    spf(:,i)=c+gas_d(:,5).*cart(i,1);%social cost considering all emissions
    tpf(:,i)=spf(:,i).*gas_d(:,2);% total cost B CNY
    gas_d=[gas_d2, tpf(:,i)];% 12 col
    %
    for j=1:l
        if gas_d(j,6)<=i% no resources
            %for carbon tax
            o_nc_o(j,i)=gas_d2(j,2);
            e_nc_o(j,i)=gas_d2(j,1);
            %
            gas_d(j,7)=0;
            gas_d(j,1)=0;
            gas_d(j,2)=0;
            gas_d(j,3)=0;    
            gas_d(j,12)=0;                
        end
        %for carbon tax
        if pf(j,i)<=0 && gas_d(j,6)>i% stranded by carbon tax
            o_pf_o(j,i)=gas_d2(j,2);
            e_pf_o(j,i)=gas_d2(j,1);
            gas_d(j,7)=0;
            gas_d(j,1)=0;
            gas_d(j,2)=0;
            gas_d(j,3)=0;    
            gas_d(j,12)=0;                            
        end
        %
    end
    gas_d_so=sortrows(gas_d,[7 8 -4]);% sorted by status, output and cost
    gasemo(:,i)=gas_d_so(:,1);% emission matrix sorted by status, output and cost
    gasomo(:,i)=gas_d_so(:,2);% output matrix sorted by status, output and cost
    gaslmo(:,i)=gas_d_so(:,3);% output matrix sorted by status, output and cost   
    gaspfmo(:,i)=gas_d_so(:,12);% output matrix sorted by status, output and cost      
    lato(:,i)=gas_d_so(:,9);
    lono(:,i)=gas_d_so(:,10);
    procodeo(:,i)=gas_d_so(:,11);
    gasemo0=gasemo;
    gasomo0=gasomo;
    gaslmo0=gaslmo;  
    gaspfmo0=gaspfmo;  
    gasomob(:,i)=gas_d_so(:,8);
    gas_r(gas_r>sum(gasomob(:,i)))=sum(gasomob(:,i));% reduction cannot beyond inventory
    for j=1:l
        if sum(gasomob(1:j,i))>=gas_r(i) && sum(gasomob(1:(j-1),i))<gas_r(i) && gasomo0(j,i)==1
            gasemo(1:(j-1),i)=0;
            gasomo(1:(j-1),i)=0;
            gaslmo(1:(j-1),i)=0;
            gaspfmo(1:(j-1),i)=0;           
            gasomo(j,i)=sum(gasomob(1:j,i))-gas_r(i);
            gasemo(j,i)=gasemo0(j,i)./gasomob(j,i).*(sum(gasomob(1:j,i))-gas_r(i));
            gaslmo(j,i)=gaslmo0(j,i)./gasomob(j,i).*(sum(gasomob(1:j,i))-gas_r(i));
            gaspfmo(j,i)=gaspfmo0(j,i)./gasomob(j,i).*(sum(gasomob(1:j,i))-gas_r(i));           
        end
    end
    %for carbon tax
    for j=1:l
        if gasemo(j,i)==0 && gas_d_so(j,7)==1
            o_fo_o(j,i)=gasomob(j,i);
            e_fo_o(j,i)=gasemo0(j,i);
        end
        if gasemo(j,i)~=0 && gasemo(j,i)<gasemo0(j,i)
            o_fo_o(j,i)=(gas_r(i)-sum(gasomob(1:j-1,i)));
            e_fo_o(j,i)=gasemo0(j,i)./gasomob(j,i).*o_fo_o(j,i);
        end
    end
    %
end
gas_so_ye=sum(gasemo,1)/1000;% yearly emissionGtCO2e
gas_so_ye_cum=cumsum(gas_so_ye);% cumulative yearly emission
gas_so_yo=sum(gasomo,1);% yearly output Mt





% sorted by emission intensity
%for carbon tax
o_nc_e=zeros(l,lper);% output reduction by natural closure
o_pf_e=zeros(l,lper);% output reduction by carbon tax-driven closure
o_fo_e=zeros(l,lper);% output reduction by directly forced closure
e_nc_e=zeros(l,lper);% emission reduction by natural closure
e_pf_e=zeros(l,lper);% emission reduction by carbon tax-driven closure
e_fo_e=zeros(l,lper);% emission reduction by directly forced closure
op=xlsread(path,'23','j2:j588');%original cost without carbon tax
cart=xlsread(path1,'carbontax','j5:k31');%carbon tax 2030-2050
c=xlsread(path,'23','k2:k588');%original cost without carbon tax and employment
pf=zeros(l,lper);% output reduction by directly forced closure
spf=zeros(l,lper);% profit with carbon tax for all emission
tpf=zeros(l,lper);
%
gaseme=zeros(l,lper);
gasome=zeros(l,lper);
gaslme=zeros(l,lper);
gaspfme=zeros(l,lper);
gasomeb=zeros(l,lper);
late=zeros(l,lper);
lone=zeros(l,lper);
procodee=zeros(l,lper);
gas_d=[gas_d0,t,gas_d0(:,2),latlon,procode];
for i=1:lper % 24-50
    %for carbon tax
    pf(:,i)=op-pei.*cart(i,1);
    spf(:,i)=c+gas_d(:,5).*cart(i,1);%social cost considering all emissions
    tpf(:,i)=spf(:,i).*gas_d(:,2);% total profit B CNY
    gas_d=[gas_d2, tpf(:,i)];% 12 col
    %
    for j=1:l
        if gas_d(j,6)<=i% no resources
            %for carbon tax
            o_nc_e(j,i)=gas_d2(j,2);
            e_nc_e(j,i)=gas_d2(j,1);
            %
            gas_d(j,7)=0;
            gas_d(j,1)=0;
            gas_d(j,2)=0;
            gas_d(j,3)=0;    
            gas_d(j,12)=0;               
        end
        %for carbon tax
        if pf(j,i)<=0 && gas_d(j,6)>i% stranded by carbon tax
            o_pf_e(j,i)=gas_d2(j,2);
            e_pf_e(j,i)=gas_d2(j,1);
            gas_d(j,7)=0;
            gas_d(j,1)=0;
            gas_d(j,2)=0;
            gas_d(j,3)=0;   
            gas_d(j,12)=0;                           
        end
        %
    end   
    gas_d_se=sortrows(gas_d,[7 -5 -4]);% sorted by status, emission intensity and cost
    gaseme(:,i)=gas_d_se(:,1);% emission matrix sorted by status, output and cost
    gasome(:,i)=gas_d_se(:,2);% output matrix sorted by status, output and cost
    gaslme(:,i)=gas_d_se(:,3);% output matrix sorted by status, output and cost   
    gaspfme(:,i)=gas_d_se(:,12);% output matrix sorted by status, output and cost      
    late(:,i)=gas_d_se(:,9);
    lone(:,i)=gas_d_se(:,10);
    procodee(:,i)=gas_d_se(:,11);
    gaseme0=gaseme;
    gasome0=gasome;
    gaslme0=gaslme;  
    gaspfme0=gaspfme;      
    gasomeb(:,i)=gas_d_se(:,8);
    gas_r(gas_r>sum(gasomeb(:,i)))=sum(gasomeb(:,i));% reduction cannot beyond inventory
    for j=1:l
        if sum(gasomeb(1:j,i))>=gas_r(i) && sum(gasomeb(1:(j-1),i))<gas_r(i) && gasome0(j,i)==1
            gaseme(1:(j-1),i)=0;
            gasome(1:(j-1),i)=0;
            gaslme(1:(j-1),i)=0;  
            gaspfme(1:(j-1),i)=0;             
            gasome(j,i)=sum(gasomeb(1:j,i))-gas_r(i);
            gaseme(j,i)=gaseme0(j,i)./gasomeb(j,i).*(sum(gasomeb(1:j,i))-gas_r(i));
            gaspfme(j,i)=gaspfme0(j,i)./gasomeb(j,i).*(sum(gasomeb(1:j,i))-gas_r(i));             
        end
    end
    %for carbon tax
    for j=1:l
        if gaseme(j,i)==0 && gas_d_se(j,7)==1
            o_fo_e(j,i)=gasomeb(j,i);
            e_fo_e(j,i)=gaseme0(j,i);
        end
        if gaseme(j,i)~=0 && gaseme(j,i)<gaseme0(j,i)
            o_fo_e(j,i)=(gas_r(i)-sum(gasomeb(1:j-1,i)));
            e_fo_e(j,i)=gaseme0(j,i)./gasomeb(j,i).*o_fo_e(j,i);
        end
    end
    %
end
gas_se_ye=sum(gaseme,1)/1000;% yearly emission GtCO2e
gas_se_ye_cum=cumsum(gas_se_ye);% cumulative yearly emission
gas_se_yo=sum(gasome,1);% yearly output Mt

f=figure;
x=2024:1:2050;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi=300;
plot(x,gas_so_ye,'-','color','k','Markersize',6,'LineWidth',5)
hold on
plot(x,gas_se_ye,'-','color','k','Markersize',6,'LineWidth',5)
ylabel('Gas (GtCO_2e)','FontSize',50);
xlim([2024 2050]);
ylim([0 1]);
xticks([2025,2030,2035,2040,2045,2050]);
set(gca,'FontSize',50);
print('F:\paper\science\figure\Emissions tax steps gas.tiff', '-dtiff', ['-r' num2str(dpi)]);
diffoe_cs=gas_so_ye_cum(1,size(gas_so_ye_cum,2))-gas_se_ye_cum(1,size(gas_se_ye_cum,2));
% output
xlswrite(pathtcs,gasomo,'gasomo','a1');
xlswrite(pathtcs,gasemo,'gasemo','a1');
xlswrite(pathtcs,gaslmo,'gaslmo','a1');
xlswrite(pathtcs,gaspfmo,'gaspfmo','a1');
xlswrite(pathtcs,gasome,'gasome','a1');
xlswrite(pathtcs,gaseme,'gaseme','a1');
xlswrite(pathtcs,gaslme,'gaslme','a1');
xlswrite(pathtcs,gaspfme,'gaspfme','a1');
xlswrite(pathtcs,gas_so_ye,'gas_so_ye','a1');
xlswrite(pathtcs,gas_se_ye,'gas_se_ye','a1');
xlswrite(pathtcs,o_nc_o,'o_nc_o','a1');
xlswrite(pathtcs,o_pf_o,'o_pf_o','a1');
xlswrite(pathtcs,o_fo_o,'o_fo_o','a1');
xlswrite(pathtcs,o_nc_e,'o_nc_e','a1');
xlswrite(pathtcs,o_pf_e,'o_pf_e','a1');
xlswrite(pathtcs,o_fo_e,'o_fo_e','a1');

xlswrite(pathtcs,e_nc_o,'e_nc_o','a1');
xlswrite(pathtcs,e_pf_o,'e_pf_o','a1');
xlswrite(pathtcs,e_fo_o,'e_fo_o','a1');
xlswrite(pathtcs,e_nc_e,'e_nc_e','a1');
xlswrite(pathtcs,e_pf_e,'e_pf_e','a1');
xlswrite(pathtcs,e_fo_e,'e_fo_e','a1');

xlswrite(pathtcs,gas_r','overview','a11');

xlswrite(pathtcs,gas_so_yo,'gas_so_yo','a1');
xlswrite(pathtcs,gas_se_yo,'gas_se_yo','a1');

xlswrite(pathtcs,lato,'lato','a1');
xlswrite(pathtcs,lono,'lono','a1');
xlswrite(pathtcs,procodeo,'procodeo','a1');
xlswrite(pathtcs,late,'late','a1');
xlswrite(pathtcs,lone,'lone','a1');
xlswrite(pathtcs,procodee,'procodee','a1');

pceo=[zeros(31,27),(1:1:31)'];
for i=1:27
    for j=1:l
        for g=1:31
            if procodeo(j,i)==pceo(g,28)
                pceo(g,i)=pceo(g,i)+gasemo(j,i);
            end
        end
    end
end
tpceo=sum(pceo(:,1:27),2);
xlswrite(pathtcs,tpceo,'tpceo','a1');

pcee=[zeros(31,27),(1:1:31)'];
for i=1:27
    for j=1:l
        for g=1:31
            if procodee(j,i)==pcee(g,28)
                pcee(g,i)=pcee(g,i)+gaseme(j,i);
            end
        end
    end
end
tpcee=sum(pcee(:,1:27),2);
xlswrite(pathtcs,tpcee,'tpcee','a1');

%composition gas_so_yo
f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(o_nc_o);sum(o_fo_o)]';
area(x,y);
ylabel0=['output reduction-gas-tax (Bm_3)'];
ylabel(ylabel0)
xlim([2024 2050]);
legend({'Natural'},'Location','northwest')
set(gca,'FontSize',50);
print('F:\paper\science\figure\output reduction composition tax steps gas.tiff', '-dtiff', ['-r' num2str(dpi)]);

%composition gas_so_ye
f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(e_nc_o);sum(e_fo_o)]'./1000;
area(x,y);
legend({'Natural'},'Location','northwest')
ylabel0=['emission reduction-gas-tax (Gt)'];
ylabel(ylabel0)
xlim([2024 2050]);
set(gca,'FontSize',50);
print('F:\paper\science\figure\emission reduction composition tax steps gas.tiff', '-dtiff', ['-r' num2str(dpi)]);


%emission intensity tco2e/t
eiso=round(gas_so_ye,6).*1000./round(gas_so_yo,6);
eise=round(gas_se_ye,6).*1000./round(gas_se_yo,6);



%% gas Announced Pledges
path=('F:\paper\science\data\gas.xlsx');
path1=('F:\paper\science\data\futureproduction.xlsx');
pathtca=('F:\paper\science\results\tax\gas\APS.xlsx');
gas_d0=xlsread(path,'23','d2:i588');% emission/yr,gas output/yr, and workforce of each gas mine
l=size(gas_d0,1);
t=ones(l,1);% 1 for sites with output
gas_d=[gas_d0,t,gas_d0(:,2),latlon,procode];
gas_d2=gas_d;
gas_d3=gas_d;
gas_ew=xlsread(path1,'gas','l5:l31');% future gas production
gas_ew2=gas_ew-41.188;% exclude gas that are not covered by inventory
gas_r=-(gas_ew2-sum(gas_d(:,2)));%reduction in inventory needed
per=2024:1:2050;% period span
lper=size(per,2);
%for carbon tax
pei=xlsread(path,'23','m2:m588');%emission inensity of production tCO2e/t
%
% sorted by output
%for carbon tax
o_nc_o=zeros(l,lper);% output reduction by natural closure
o_pf_o=zeros(l,lper);% output reduction by carbon tax-driven closure
o_fo_o=zeros(l,lper);% output reduction by directly forced closure
e_nc_o=zeros(l,lper);% emission reduction by natural closure
e_pf_o=zeros(l,lper);% emission reduction by carbon tax-driven closure
e_fo_o=zeros(l,lper);% emission reduction by directly forced closure
op=xlsread(path,'23','j2:j588');%original profit without carbon tax
cart=xlsread(path1,'carbontax','j5:k31');%carbon tax 2030-2050
c=xlsread(path,'23','k2:k588');%original cost without carbon tax and employment
pf=zeros(l,lper);% output reduction by directly forced closure
spf=zeros(l,lper);% profit with carbon tax for all emission
tpf=zeros(l,lper);
%
gasemo=zeros(l,lper);
gasomo=zeros(l,lper);
gaslmo=zeros(l,lper);
gaspfmo=zeros(l,lper);
gasomob=zeros(l,lper);
lato=zeros(l,lper);
lono=zeros(l,lper);
procodeo=zeros(l,lper);
for i=1:lper % 24-50
    %for carbon tax
    pf(:,i)=op-pei.*cart(i,2);
    spf(:,i)=c+gas_d(:,5).*cart(i,2);%social cost considering all emissions
    tpf(:,i)=spf(:,i).*gas_d(:,2);% total profit B CNY
    gas_d=[gas_d2, tpf(:,i)];% 12 col
    %
    for j=1:l
        if gas_d(j,6)<=i% no resources
            %for carbon tax
            o_nc_o(j,i)=gas_d2(j,2);
            e_nc_o(j,i)=gas_d2(j,1);
            %
            gas_d(j,7)=0;
            gas_d(j,1)=0;
            gas_d(j,2)=0;
            gas_d(j,3)=0;    
            gas_d(j,12)=0;    
           
        end
        %for carbon tax
        if pf(j,i)<=0 && gas_d(j,6)>i% stranded by carbon tax
            o_pf_o(j,i)=gas_d2(j,2);
            e_pf_o(j,i)=gas_d2(j,1);
            gas_d(j,7)=0;
            gas_d(j,1)=0;
            gas_d(j,2)=0;
            gas_d(j,3)=0;      
            gas_d(j,12)=0;                
        end
        %
    end
    gas_d_so=sortrows(gas_d,[7 8 -4]);% sorted by status, output and cost
    gasemo(:,i)=gas_d_so(:,1);% emission matrix sorted by status, output and cost
    gasomo(:,i)=gas_d_so(:,2);% output matrix sorted by status, output and cost
    gaslmo(:,i)=gas_d_so(:,3);% output matrix sorted by status, output and cost    
    gaspfmo(:,i)=gas_d_so(:,12);% output matrix sorted by status, output and cost       
    lato(:,i)=gas_d_so(:,9);
    lono(:,i)=gas_d_so(:,10);
    procodeo(:,i)=gas_d_so(:,11);
    gasemo0=gasemo;
    gasomo0=gasomo;
    gaslmo0=gaslmo;   
    gaspfmo0=gaspfmo;       
    gasomob(:,i)=gas_d_so(:,8);
    gas_r(gas_r>sum(gasomob(:,i)))=sum(gasomob(:,i));% reduction cannot beyond inventory
    for j=1:l
        if sum(gasomob(1:j,i))>=gas_r(i) && sum(gasomob(1:(j-1),i))<gas_r(i) && gasomo0(j,i)==1
            gasemo(1:(j-1),i)=0;
            gasomo(1:(j-1),i)=0;
            gaslmo(1:(j-1),i)=0;       
            gaspfmo(1:(j-1),i)=0;                  
            gasomo(j,i)=sum(gasomob(1:j,i))-gas_r(i);
            gasemo(j,i)=gasemo0(j,i)./gasomob(j,i).*(sum(gasomob(1:j,i))-gas_r(i));
            gaslmo(j,i)=gaslmo0(j,i)./gasomob(j,i).*(sum(gasomob(1:j,i))-gas_r(i));  
            gaspfmo(j,i)=gaspfmo0(j,i)./gasomob(j,i).*(sum(gasomob(1:j,i))-gas_r(i));             
        end
    end
    %for carbon tax
    for j=1:l
        if gasemo(j,i)==0 && gas_d_so(j,7)==1
            o_fo_o(j,i)=gasomob(j,i);
            e_fo_o(j,i)=gasemo0(j,i);
        end
        if gasemo(j,i)~=0 && gasemo(j,i)<gasemo0(j,i)
            o_fo_o(j,i)=(gas_r(i)-sum(gasomob(1:j-1,i)));
            e_fo_o(j,i)=gasemo0(j,i)./gasomob(j,i).*o_fo_o(j,i);
        end
    end
    %
end
gas_so_ye=sum(gasemo,1)/1000;% yearly emissionGtCO2e
gas_so_ye_cum=cumsum(gas_so_ye);% cumulative yearly emission
gas_so_yo=sum(gasomo,1);% yearly output Mt



% sorted by emission intensity
%for carbon tax
o_nc_e=zeros(l,lper);% output reduction by natural closure
o_pf_e=zeros(l,lper);% output reduction by carbon tax-driven closure
o_fo_e=zeros(l,lper);% output reduction by directly forced closure
e_nc_e=zeros(l,lper);% emission reduction by natural closure
e_pf_e=zeros(l,lper);% emission reduction by carbon tax-driven closure
e_fo_e=zeros(l,lper);% emission reduction by directly forced closure
op=xlsread(path,'23','j2:j588');%original cost without carbon tax
cart=xlsread(path1,'carbontax','j5:k31');%carbon tax 2030-2050
c=xlsread(path,'23','k2:k588');%original cost without carbon tax and employment
pf=zeros(l,lper);% output reduction by directly forced closure
spf=zeros(l,lper);% profit with carbon tax for all emission
tpf=zeros(l,lper);
%
gaseme=zeros(l,lper);
gasome=zeros(l,lper);
gaslme=zeros(l,lper);
gaspfme=zeros(l,lper);
gasomeb=zeros(l,lper);
late=zeros(l,lper);
lone=zeros(l,lper);
procodee=zeros(l,lper);
gas_d=[gas_d0,t,gas_d0(:,2),latlon,procode];
for i=1:lper % 24-50
    %for carbon tax
    pf(:,i)=op-pei.*cart(i,2);
    spf(:,i)=c+gas_d(:,5).*cart(i,2);%social cost considering all emissions
    tpf(:,i)=spf(:,i).*gas_d(:,2);% total profit B CNY
    gas_d=[gas_d2, tpf(:,i)];% 12 col
    %
    for j=1:l
        if gas_d(j,6)<=i% no resources
            %for carbon tax
            o_nc_e(j,i)=gas_d2(j,2);
            e_nc_e(j,i)=gas_d2(j,1);
            %
            gas_d(j,7)=0;
            gas_d(j,1)=0;
            gas_d(j,2)=0;
            gas_d(j,3)=0; 
            gas_d(j,12)=0;             
        end
        %for carbon tax
        if pf(j,i)<=0 && gas_d(j,6)>i% stranded by carbon tax
            o_pf_e(j,i)=gas_d2(j,2);
            e_pf_e(j,i)=gas_d2(j,1);
            gas_d(j,7)=0;
            gas_d(j,1)=0;
            gas_d(j,2)=0;
            gas_d(j,3)=0; 
            gas_d(j,12)=0;             
        end
        %
    end   
    gas_d_se=sortrows(gas_d,[7 -5 -4]);% sorted by status, emission intensity and cost
    gaseme(:,i)=gas_d_se(:,1);% emission matrix sorted by status, output and cost
    gasome(:,i)=gas_d_se(:,2);% output matrix sorted by status, output and cost
    gaslme(:,i)=gas_d_se(:,3);% output matrix sorted by status, output and cost 
    gaspfme(:,i)=gas_d_se(:,12);% output matrix sorted by status, output and cost     
    late(:,i)=gas_d_se(:,9);
    lone(:,i)=gas_d_se(:,10);
    procodee(:,i)=gas_d_se(:,11);
    gaseme0=gaseme;
    gasome0=gasome;
    gaslme0=gaslme;  
    gaspfme0=gaspfme;     
    gasomeb(:,i)=gas_d_se(:,8);
    gas_r(gas_r>sum(gasomeb(:,i)))=sum(gasomeb(:,i));% reduction cannot beyond inventory
    for j=1:l
        if sum(gasomeb(1:j,i))>=gas_r(i) && sum(gasomeb(1:(j-1),i))<gas_r(i) && gasome0(j,i)==1
            gaseme(1:(j-1),i)=0;
            gasome(1:(j-1),i)=0;
            gaslme(1:(j-1),i)=0;    
            gaspfme(1:(j-1),i)=0;    
            gasome(j,i)=sum(gasomeb(1:j,i))-gas_r(i);
            gaseme(j,i)=gaseme0(j,i)./gasomeb(j,i).*(sum(gasomeb(1:j,i))-gas_r(i));
            gaslme(j,i)=gaslme0(j,i)./gasomeb(j,i).*(sum(gasomeb(1:j,i))-gas_r(i));  
            gaspfme(j,i)=gaspfme0(j,i)./gasomeb(j,i).*(sum(gasomeb(1:j,i))-gas_r(i));              
        end
    end
    %for carbon tax
    for j=1:l
        if gaseme(j,i)==0 && gas_d_se(j,7)==1
            o_fo_e(j,i)=gasomeb(j,i);
            e_fo_e(j,i)=gaseme0(j,i);
        end
        if gaseme(j,i)~=0 && gaseme(j,i)<gaseme0(j,i)
            o_fo_e(j,i)=(gas_r(i)-sum(gasomeb(1:j-1,i)));
            e_fo_e(j,i)=gaseme0(j,i)./gasomeb(j,i).*o_fo_e(j,i);
        end
    end
    %
end
gas_se_ye=sum(gaseme,1)/1000;% yearly emission GtCO2e
gas_se_ye_cum=cumsum(gas_se_ye);% cumulative yearly emission
gas_se_yo=sum(gasome,1);% yearly output Mt

f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
x=2024:1:2050;
plot(x,gas_so_ye,'-','color','k','Markersize',6,'LineWidth',5)
hold on
plot(x,gas_se_ye,'-','color','k','Markersize',6,'LineWidth',5)
ylabel('Gas (GtCO_2e)','FontSize',50);
xlim([2024 2050]);
ylim([0 1]);
xticks([2025,2030,2035,2040,2045,2050]);
set(gca,'FontSize',100);
diffoe_cs=gas_so_ye_cum(1,size(gas_so_ye_cum,2))-gas_se_ye_cum(1,size(gas_se_ye_cum,2));
% output
xlswrite(pathtca,gasomo,'gasomo','a1');
xlswrite(pathtca,gasemo,'gasemo','a1');
xlswrite(pathtca,gaslmo,'gaslmo','a1');
xlswrite(pathtca,gaspfmo,'gaspfmo','a1');
xlswrite(pathtca,gasome,'gasome','a1');
xlswrite(pathtca,gaseme,'gaseme','a1');
xlswrite(pathtca,gaslme,'gaslme','a1');
xlswrite(pathtca,gaspfme,'gaspfme','a1');
xlswrite(pathtca,gas_so_ye,'gas_so_ye','a1');
xlswrite(pathtca,gas_se_ye,'gas_se_ye','a1');
xlswrite(pathtca,o_nc_o,'o_nc_o','a1');
xlswrite(pathtca,o_pf_o,'o_pf_o','a1');
xlswrite(pathtca,o_fo_o,'o_fo_o','a1');
xlswrite(pathtca,o_nc_e,'o_nc_e','a1');
xlswrite(pathtca,o_pf_e,'o_pf_e','a1');
xlswrite(pathtca,o_fo_e,'o_fo_e','a1');

xlswrite(pathtca,e_nc_o,'e_nc_o','a1');
xlswrite(pathtca,e_pf_o,'e_pf_o','a1');
xlswrite(pathtca,e_fo_o,'e_fo_o','a1');
xlswrite(pathtca,e_nc_e,'e_nc_e','a1');
xlswrite(pathtca,e_pf_e,'e_pf_e','a1');
xlswrite(pathtca,e_fo_e,'e_fo_e','a1');

xlswrite(pathtca,gas_r','overview','a11');
xlswrite(pathtca,gas_so_yo,'gas_so_yo','a1');
xlswrite(pathtca,gas_se_yo,'gas_se_yo','a1');

xlswrite(pathtca,lato,'lato','a1');
xlswrite(pathtca,lono,'lono','a1');
xlswrite(pathtca,procodeo,'procodeo','a1');
xlswrite(pathtca,late,'late','a1');
xlswrite(pathtca,lone,'lone','a1');
xlswrite(pathtca,procodee,'procodee','a1');

pceo=[zeros(31,27),(1:1:31)'];
for i=1:27
    for j=1:l
        for g=1:31
            if procodeo(j,i)==pceo(g,28)
                pceo(g,i)=pceo(g,i)+gasemo(j,i);
            end
        end
    end
end
tpceo=sum(pceo(:,1:27),2);
xlswrite(pathtca,tpceo,'tpceo','a1');

pcee=[zeros(31,27),(1:1:31)'];
for i=1:27
    for j=1:l
        for g=1:31
            if procodee(j,i)==pcee(g,28)
                pcee(g,i)=pcee(g,i)+gaseme(j,i);
            end
        end
    end
end
tpcee=sum(pcee(:,1:27),2);
xlswrite(pathtca,tpcee,'tpcee','a1');

%composition gas_so_yo
f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(o_nc_o);sum(o_fo_o)]';
area(x,y);
ylabel0=['APS' newline 'output reduction-gas (Bt)'];
ylabel(ylabel0)
xlim([2024 2050]);
legend({'Natural'},'Location','northwest')
set(gca,'FontSize',50);
print('F:\paper\science\figure\output reduction composition tax APS gas.tiff', '-dtiff', ['-r' num2str(dpi)]);

%composition gas_so_ye
f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(e_nc_o);sum(e_fo_o)]'./1000;
area(x,y);
legend({'Natural'},'Location','northwest')
ylabel0=['APS-PSO' newline 'emission reduction-gas (Mt)'];
ylabel(ylabel0)
xlim([2024 2050]);
set(gca,'FontSize',50);
print('F:\paper\science\figure\emission reduction composition tax APS Scale gas.tiff', '-dtiff', ['-r' num2str(dpi)]);

%composition gas_se_ye
f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
y=[sum(e_nc_e);sum(e_fo_e)]'./1000;
area(x,y);
legend({'Natural'},'Location','northwest')
ylabel0=['APS-PHEI' newline 'emission reduction-gas (Mt)'];
ylabel(ylabel0)
xlim([2024 2050]);
set(gca,'FontSize',50);
print('F:\paper\science\figure\emission reduction composition tax APS Emission gas.tiff', '-dtiff', ['-r' num2str(dpi)]);


%emission intensity tco2e/t
eiao=round(gas_so_ye,6).*1000./round(gas_so_yo,6);
eiae=round(gas_se_ye,6).*1000./round(gas_se_yo,6);

f=figure;
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
x=2024:1:2050;
plot(x,eiso,'-','color','k','Markersize',6,'LineWidth',5)
%hold on
%plot(x,eise,'-','color',[235, 23, 23]./255,'Markersize',6,'LineWidth',5)
%hold on
%plot(x,eiao,'-','color',[6, 107, 82]./255,'Markersize',6,'LineWidth',5)
%hold on
%plot(x,eiae,'-','color',[6, 107, 82]./255,'Markersize',6,'LineWidth',5)
ylabel('Emission intensity-gas-tax (tCO_2e/t)','FontSize',45);
xlim([2024 2050]);
ylim([3 4]);
xticks([2025,2030,2035,2040,2045,2050]);
%lgd = legend('STEPS-PSO','STEPS-PHEI','APS-PSO','APS-PHEI','Location','southwest','Orientation','vertical');
set(gca,'FontSize',45);
print('F:\paper\science\figure\Emission intensity-gas-tax.tiff', '-dtiff', ['-r' num2str(dpi)]);

%% output and emission intensity of gas
f=figure;
mycolors = [0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.8500 0.3250 0.0980];
set(f,'defaultAxesColorOrder',mycolors)
set(f, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
dpi = 300;
ee=gas_d0(:,5).*1390;
plot(gas_d0(:,2),ee,'ok','Markersize',15,'LineWidth',2)
ylabel('tCO_2e/t','FontSize',110);
xlabel('Output of gas field (Gm^3)','FontSize',110);
xlim([0 31]);
ylim([2.5 4.5]);
set(gca,'FontSize',110);
print('F:\paper\science\figure\gas output emission intensity.tiff', '-dtiff', ['-r' num2str(dpi)]);

close all
toc