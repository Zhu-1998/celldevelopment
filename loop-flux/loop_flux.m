clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mouse retina development %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=[1-(1/89+1/188+1/168) 1/89 1/188 1/168; 1/280 1-(1/280+1/408+1/275) 1/408 1/275; 1/508 1/300 1-(1/508+1/300+1/201) 1/201; 1/356 1/211 1/211 1-(1/356+1/211+1/211)];

n=size(M,1);
E=eye(n);
MM=[(M'-E);ones(1,n)];
N=[zeros(n,1);1];
P=MM\N;
% P=inv(MM)*N;

% P'*M

P=[0.1128;0.3541;0.2660;0.2671];

% for i=1:4
%     for j=1:4
%         C(i,j)=max(M(i,j)*P(i,1)-M(j,i)*P(j,1),0)/P(i,1);
%     end
% end
% 
% for i=1:4
%     C(i,i)=-(C(i,1)+C(i,2)+C(i,3)+C(i,4));
% end
% 
% 
% for i=1:4
%     for j=1:4
%         J(i,j)=C(i,j)*P(i,1);
%     end
% end
% 
% for i=1:4
%     J(i,i)=0;
% end


for i=1:4
    for j=1:4
        PJ(i,j)=M(i,j)*P(i,1);
    end
end


for i=1:4
    for j=1:4
        J(i,j)=M(i,j)*P(i,1)-min(M(i,j)*P(i,1), M(j,i)*P(j,1));
    end
end


J1241=min([J(1,2) J(2,4) J(4,1)]);

J_1=J;
J_1(1,2)=J_1(1,2)-J1241;
J_1(2,4)=J_1(2,4)-J1241;
J_1(4,1)=J_1(4,1)-J1241;


J13241=min([J_1(1,3) J_1(3,2) J_1(2,4) J_1(4,1)]);
J_2=J_1;
J_2(1,3)=J_2(1,3)-J13241;
J_2(3,2)=J_2(3,2)-J13241;
J_2(2,4)=J_2(2,4)-J13241;
J_2(4,1)=J_2(4,1)-J13241;

J1341=min([J_2(1,3) J_2(3,4) J_2(4,1)]);
J_3=J_2;
J_3(1,3)=J_3(1,3)-J1341;
J_3(3,4)=J_3(3,4)-J1341;
J_3(4,1)=J_3(4,1)-J1341;


xname = {'Progenitor','PR','AC/HC','RGC'};
yname = {'Progenitor','PR','AC/HC','RGC'};
h = heatmap(xname,yname,M);
h.CellLabelFormat = '%0.5f';
colormap(gca, 'parula')

b = heatmap(xname,yname,PJ);
colormap(othercolor('Greens3'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hsc %%%%%%%%%%%%%%%%%%%%%%%%
M_hsc=[1-(1/31+1/46+1/46+1/110+1/310) 1/31 1/46 1/46 1/110 1/310; 
   1/130 1-(1/130+1/120+1/130+1/150+1/180) 1/120 1/130 1/150 1/180; 
   1/150 1/55 1-(1/150+1/55+1/89+1/53+1/160) 1/89 1/53 1/160;
   1/270 1/100 1/210 1-(1/270+1/100+1/210+1/190+1/220) 1/190 1/220;
   1/120 1/66 1/55 1/72 1-(1/120+1/66+1/55+1/72+1/140) 1/140;
   1/320 1/110 1/210 1/160 1/180 1-(1/320+1/110+1/210+1/160+1/180)];

n_hsc=size(M_hsc,1);
E_hsc=eye(n_hsc);
MM_hsc=[(M_hsc'-E_hsc);ones(1,n_hsc)];
N_hsc=[zeros(n_hsc,1);1];
P_hsc=MM_hsc\N_hsc;
% P_hsc=inv(MM_hsc)*N_hsc;

% P_hsc'*M_hsc

P_hsc=[0.0612578948450237;
    0.279292098764672;
    0.125975439006672;
    0.263281217198114;
    0.112557516213203;
    0.157635833972315];

for i=1:6
    for j=1:6
        PJ_hsc(i,j)=M_hsc(i,j)*P_hsc(i,1);
    end
end

xname = {'HSC','Meg','Ery','Bas','Mon','Neu'};
yname = {'HSC','Meg','Ery','Bas','Mon','Neu'};
h = heatmap(xname,yname,M_hsc);
h.CellLabelFormat = '%0.5f';
hold on

figure(2)
b = heatmap(xname,yname,PJ_hsc);
b.CellLabelFormat = '%0.5f';
colormap(othercolor('Greens3'))



for i=1:6
    for j=1:6
        J_hsc(i,j)=M_hsc(i,j)*P_hsc(i,1)-min(M_hsc(i,j)*P_hsc(i,1), M_hsc(j,i)*P_hsc(j,1));
    end
end


J_hsc13421=min([J_hsc(1,3) J_hsc(3,4) J_hsc(4,2) J_hsc(2,1)]);
J_hsc_1=J_hsc;
J_hsc_1(1,3)=J_hsc_1(1,3)-J_hsc13421;
J_hsc_1(3,4)=J_hsc_1(3,4)-J_hsc13421;
J_hsc_1(4,2)=J_hsc_1(4,2)-J_hsc13421;
J_hsc_1(2,1)=J_hsc_1(2,1)-J_hsc13421;

J_hsc1351=min([J_hsc_1(1,3) J_hsc_1(3,5) J_hsc_1(5,1)]);
J_hsc_2=J_hsc_1;
J_hsc_2(1,3)=J_hsc_2(1,3)-J_hsc1351;
J_hsc_2(3,5)=J_hsc_2(3,5)-J_hsc1351;
J_hsc_2(5,1)=J_hsc_2(5,1)-J_hsc1351;

J_hsc1421=min([J_hsc_2(1,4) J_hsc_2(4,2) J_hsc_2(2,1)]);
J_hsc_3=J_hsc_2;
J_hsc_3(1,4)=J_hsc_3(1,4)-J_hsc1421;
J_hsc_3(4,2)=J_hsc_3(4,2)-J_hsc1421;
J_hsc_3(2,1)=J_hsc_3(2,1)-J_hsc1421;

J_hsc142351=min([J_hsc_3(1,4) J_hsc_3(4,2) J_hsc_3(2,3) J_hsc_3(3,5) J_hsc_3(5,1)]);
J_hsc_4=J_hsc_3;
J_hsc_4(1,4)=J_hsc_4(1,4)-J_hsc142351;
J_hsc_4(4,2)=J_hsc_4(4,2)-J_hsc142351;
J_hsc_4(2,3)=J_hsc_4(2,3)-J_hsc142351;
J_hsc_4(3,5)=J_hsc_4(3,5)-J_hsc142351;
J_hsc_4(5,1)=J_hsc_4(5,1)-J_hsc142351;

J_hsc142361=min([J_hsc_4(1,4) J_hsc_4(4,2) J_hsc_4(2,3) J_hsc_4(3,6) J_hsc_4(6,1)]);
J_hsc_5=J_hsc_4;
J_hsc_5(1,4)=J_hsc_5(1,4)-J_hsc142361;
J_hsc_5(4,2)=J_hsc_5(4,2)-J_hsc142361;
J_hsc_5(2,3)=J_hsc_5(2,3)-J_hsc142361;
J_hsc_5(3,6)=J_hsc_5(3,6)-J_hsc142361;
J_hsc_5(6,1)=J_hsc_5(6,1)-J_hsc142361;

J_hsc14251=min([J_hsc_5(1,4) J_hsc_5(4,2) J_hsc_5(2,5) J_hsc_5(5,1)]);
J_hsc_6=J_hsc_5;
J_hsc_6(1,4)=J_hsc_6(1,4)-J_hsc14251;
J_hsc_6(4,2)=J_hsc_6(4,2)-J_hsc14251;
J_hsc_6(2,5)=J_hsc_6(2,5)-J_hsc14251;
J_hsc_6(5,1)=J_hsc_6(5,1)-J_hsc14251;

J_hsc4254=min([J_hsc_6(4,2) J_hsc_6(2,5) J_hsc_6(5,4)]);
J_hsc_7=J_hsc_6;
J_hsc_7(4,2)=J_hsc_7(4,2)-J_hsc4254;
J_hsc_7(2,5)=J_hsc_7(2,5)-J_hsc4254;
J_hsc_7(5,4)=J_hsc_7(5,4)-J_hsc4254;

J_hsc14261=min([J_hsc_7(1,4) J_hsc_7(4,2) J_hsc_7(2,6) J_hsc_7(6,1)]);
J_hsc_8=J_hsc_7;
J_hsc_8(1,4)=J_hsc_8(1,4)-J_hsc14261;
J_hsc_8(4,2)=J_hsc_8(4,2)-J_hsc14261;
J_hsc_8(2,6)=J_hsc_8(2,6)-J_hsc14261;
J_hsc_8(6,1)=J_hsc_8(6,1)-J_hsc14261;

J_hsc1461=min([J_hsc_8(1,4) J_hsc_8(4,6) J_hsc_8(6,1)]);
J_hsc_9=J_hsc_8;
J_hsc_9(1,4)=J_hsc_9(1,4)-J_hsc1461;
J_hsc_9(4,6)=J_hsc_9(4,6)-J_hsc1461;
J_hsc_9(6,1)=J_hsc_9(6,1)-J_hsc1461;

J_hsc4654=min([J_hsc_9(4,6) J_hsc_9(6,5) J_hsc_9(5,4)]);
J_hsc_10=J_hsc_9;
J_hsc_10(4,6)=J_hsc_10(4,6)-J_hsc4654;
J_hsc_10(6,5)=J_hsc_10(6,5)-J_hsc4654;
J_hsc_10(5,4)=J_hsc_10(5,4)-J_hsc4654;

