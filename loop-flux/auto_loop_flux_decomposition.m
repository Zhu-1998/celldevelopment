%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Please cite: Proc Natl Acad Sci U S A. 2024;121(34):e2401540121. doi: 10.1073/pnas.2401540121.  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mouse retina development %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T=[0 1/89 1/188 1/168; 
    1/280 0 1/408 1/275; 
    1/508 1/300 0 1/201; 
    1/356 1/211 1/211 0];
M = T;
row_sums = sum(M, 2); 
for i=1:size(M,1)
    M(i,i) = -row_sums(i);
end

M=[-(1/89+1/188+1/168) 1/89 1/188 1/168; 
    1/280 -(1/280+1/408+1/275) 1/408 1/275; 
    1/508 1/300 -(1/508+1/300+1/201) 1/201; 
    1/356 1/211 1/211 -(1/356+1/211+1/211)];

row_sums_new = sum(M, 2);

n=size(M,1);
MM=[(M');ones(1,n)];
N=[zeros(n,1);1];
P=MM\N;
% P=inv(MM)*N;
sum(P)



for i=1:size(T,1)
    for j=1:size(T,1)
        C(i,j)=max(T(i,j)*P(i,1)-T(j,i)*P(j,1), 0)/P(i,1);
    end
end
row_sums_C = sum(C, 2); 
for i=1:size(C,1)
    C(i,i) = -(row_sums_C(i)-C(i,i));
end
row_sums_C_new = sum(C, 2); 

for i=1:size(T,1)
    for j=1:size(T,1)
        D(i,j)=min(T(i,j)*P(i,1), T(j,i)*P(j,1))/P(i,1);
    end
end
row_sums_D = sum(D, 2); 
for i=1:size(D,1)
    D(i,i) = -(row_sums_D(i)-D(i,i));
end
row_sums_D_new = sum(D, 2); 


for i=1:size(T,1)
    for j=1:size(T,1)
    F_ss(i,j)=T(j,i)*P(j,1)-T(i,j)*P(i,1);
    end
end
row_sums_F_ss = sum(F_ss, 2); 

for i=1:size(C,1)
    for j=1:size(C,1)
        F(i,j)=C(j,i)*P(j,1)-C(i,j)*P(i,1);
    end
end

for i=1:size(C,1)
    for j=1:size(C,1)
        J(i,j)=C(i,j)*P(i,1);
    end
end
for i=1:size(J,1)
    J(i,i)=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xname = {'Progenitor','PR','AC/HC','RGC'};
yname = {'Progenitor','PR','AC/HC','RGC'};
h = heatmap(xname,yname,M);
h.CellLabelFormat = '%0.5f';
colormap(gca, 'parula')

b = heatmap(xname,yname,F_ss);
colormap(othercolor('Greens3'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% loop-flux decomposition %%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%% auto loop-flux decomposition %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position = [];
for i = 1:size(J, 1)
    for j = 1:size(J, 2)
        if J(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position
positions_record = [];  % To record each position

% Loop through each element in the specified row
for j = 1:size(J, 2)  % Iterate through columns of the specified row
    if J(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)
    % Ensure found_position(2) is within matrix bounds
    if isempty(found_position) || found_position(2) > size(J, 1)
        break;  % Exit the loop if not found or out of bounds
    end

    % Continue looping through the found_position(2) row
    for j = 1:size(J, 2)
        if J(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end

J1241=min([J(positions_record(1,1),positions_record(1,2)) J(positions_record(2,1),positions_record(2,2)) J(positions_record(3,1),positions_record(3,1))]);
J_1=J;
J_1(positions_record(1,1),positions_record(1,2))=J_1(positions_record(1,1),positions_record(1,2))-J1241;
J_1(positions_record(2,1),positions_record(2,2))=J_1(positions_record(2,1),positions_record(2,2))-J1241;
J_1(positions_record(3,1),positions_record(3,1))=J_1(positions_record(3,1),positions_record(3,1))-J1241;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position = [];
for i = 1:size(J_1, 1)
    for j = 1:size(J_1, 2)
        if J_1(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_1, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position

% Loop through each element in the specified row
for j = 1:size(J_1, 2)  % Iterate through columns of the specified row
    if J_1(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)
    % Ensure found_position(2) is within matrix bounds
    if isempty(found_position) || found_position(2) > size(J_1, 1)
        break;  % Exit the loop if not found or out of bounds
    end

    % Continue looping through the found_position(2) row
    for j = 1:size(J_1, 2)
        if J_1(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end

J13241=min([J_1(positions_record(4,1),positions_record(4,2)) J_1(positions_record(5,1),positions_record(5,2)) J_1(positions_record(6,1),positions_record(6,2)) J_1(positions_record(7,1),positions_record(7,2))]);
J_2=J_1;
J_2(positions_record(4,1),positions_record(4,2))=J_2(positions_record(4,1),positions_record(4,2))-J13241;
J_2(positions_record(5,1),positions_record(5,2))=J_2(positions_record(5,1),positions_record(5,2))-J13241;
J_2(positions_record(6,1),positions_record(6,2))=J_2(positions_record(6,1),positions_record(6,2))-J13241;
J_2(positions_record(7,1),positions_record(7,2))=J_2(positions_record(7,1),positions_record(7,2))-J13241;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position = [];
for i = 1:size(J_2, 1)
    for j = 1:size(J_2, 2)
        if J_2(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_2, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position

% Loop through each element in the specified row
for j = 1:size(J_2, 2)  % Iterate through columns of the specified row
    if J_2(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)
    % Ensure found_position(2) is within matrix bounds
    if isempty(found_position) || found_position(2) > size(J_2, 1)
        break;  % Exit the loop if not found or out of bounds
    end

    % Continue looping through the found_position(2) row
    for j = 1:size(J_2, 2)
        if J_2(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end

J1341=min([J_2(positions_record(8,1),positions_record(8,2)) J_2(positions_record(9,1),positions_record(9,2)) J_2(positions_record(10,1),positions_record(10,2))]);
J_3=J_2;
J_3(positions_record(8,1),positions_record(8,2))=J_3(positions_record(8,1),positions_record(8,2))-J1341;
J_3(positions_record(9,1),positions_record(9,2))=J_3(positions_record(9,1),positions_record(9,2))-J1341;
J_3(positions_record(10,1),positions_record(10,2))=J_3(positions_record(10,1),positions_record(10,2))-J1341;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position = [];
for i = 1:size(J_3, 1)
    for j = 1:size(J_3, 2)
        if J_3(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_3, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position

% Loop through each element in the specified row
for j = 1:size(J_3, 2)  % Iterate through columns of the specified row
    if J_3(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)
    % Ensure found_position(2) is within matrix bounds
    if isempty(found_position) || found_position(2) > size(J_3, 1)
        break;  % Exit the loop if not found or out of bounds
    end

    % Continue looping through the found_position(2) row
    for j = 1:size(J_3, 2)
        if J_3(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hsc %%%%%%%%%%%%%%%%%%%%%%%%
clear all;

T_hsc=[0 1/31 1/46 1/46 1/110 1/310; 
   1/130 0 1/120 1/130 1/150 1/180; 
   1/150 1/55 0 1/89 1/53 1/160;
   1/270 1/100 1/210 0 1/190 1/220;
   1/120 1/66 1/55 1/72 0 1/140;
   1/320 1/110 1/210 1/160 1/180 0];
M_hsc = T_hsc;
row_sums = sum(M_hsc, 2); 
for i=1:size(M_hsc,1)
    M_hsc(i,i) = -row_sums(i);
end
row_sums_new = sum(M_hsc, 2); 

M_hsc=[-(1/31+1/46+1/46+1/110+1/310) 1/31 1/46 1/46 1/110 1/310; 
   1/130 -(1/130+1/120+1/130+1/150+1/180) 1/120 1/130 1/150 1/180; 
   1/150 1/55 -(1/150+1/55+1/89+1/53+1/160) 1/89 1/53 1/160;
   1/270 1/100 1/210 -(1/270+1/100+1/210+1/190+1/220) 1/190 1/220;
   1/120 1/66 1/55 1/72 -(1/120+1/66+1/55+1/72+1/140) 1/140;
   1/320 1/110 1/210 1/160 1/180 -(1/320+1/110+1/210+1/160+1/180)];

n_hsc=size(M_hsc,1);
MM_hsc=[(M_hsc');ones(1,n_hsc)];
N_hsc=[zeros(n_hsc,1);1];
P_hsc=MM_hsc\N_hsc;
% P_hsc=inv(MM_hsc)*N_hsc;
sum(P_hsc)

for i=1:size(T_hsc,1)
    for j=1:size(T_hsc,1)
    F_ss_hsc(i,j)=T_hsc(j,i)*P_hsc(j,1)-T_hsc(i,j)*P_hsc(i,1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xname = {'HSC','Meg','Ery','Bas','Mon','Neu'};
yname = {'HSC','Meg','Ery','Bas','Mon','Neu'};
h = heatmap(xname,yname,M_hsc);
h.CellLabelFormat = '%0.5f';
hold on

figure(2)
b = heatmap(xname,yname,F_ss_hsc);
b.CellLabelFormat = '%0.5f';
colormap(othercolor('Greens3'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:size(T_hsc,1)
    for j=1:size(T_hsc,1)
        J_hsc(i,j)=T_hsc(i,j)*P_hsc(i,1)-min(T_hsc(i,j)*P_hsc(i,1), T_hsc(j,i)*P_hsc(j,1));
    end
end
for i=1:size(J_hsc,1)
    J_hsc(i,i)=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% loop flux decomposition %%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%% auto loop-flux decomposition %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position = [];
for i = 1:size(J_hsc, 1)
    for j = 1:size(J_hsc, 2)
        if J_hsc(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_hsc, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position
positions_record = [];  % To record each position

% Loop through each element in the specified row
for j = 1:size(J_hsc, 2)  % Iterate through columns of the specified row
    if J_hsc(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)

    % Continue looping through the found_position(2) row
    for j = 1:size(J_hsc, 2)
        if J_hsc(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end

J_hsc13421=min([J_hsc(positions_record(1,1),positions_record(1,2)) J_hsc(positions_record(2,1),positions_record(2,2)) J_hsc(positions_record(3,1),positions_record(3,2)) J_hsc(positions_record(4,1),positions_record(4,2))]);
J_hsc_1=J_hsc;
J_hsc_1(positions_record(1,1),positions_record(1,2))=J_hsc_1(positions_record(1,1),positions_record(1,2))-J_hsc13421;
J_hsc_1(positions_record(2,1),positions_record(2,2))=J_hsc_1(positions_record(2,1),positions_record(2,2))-J_hsc13421;
J_hsc_1(positions_record(3,1),positions_record(3,2))=J_hsc_1(positions_record(3,1),positions_record(3,2))-J_hsc13421;
J_hsc_1(positions_record(4,1),positions_record(4,2))=J_hsc_1(positions_record(4,1),positions_record(4,2))-J_hsc13421;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position = [];
for i = 1:size(J_hsc_1, 1)
    for j = 1:size(J_hsc_1, 2)
        if J_hsc_1(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_hsc_1, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position

% Loop through each element in the specified row
for j = 1:size(J_hsc_1, 2)  % Iterate through columns of the specified row
    if J_hsc_1(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)

    % Continue looping through the found_position(2) row
    for j = 1:size(J_hsc_1, 2)
        if J_hsc_1(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end

J_hsc1351=min([J_hsc_1(positions_record(5,1),positions_record(5,2)) J_hsc_1(positions_record(6,1),positions_record(6,2)) J_hsc_1(positions_record(7,1),positions_record(7,2))]);
J_hsc_2=J_hsc_1;
J_hsc_2(positions_record(5,1),positions_record(5,2))=J_hsc_2(positions_record(5,1),positions_record(5,2))-J_hsc1351;
J_hsc_2(positions_record(6,1),positions_record(6,2))=J_hsc_2(positions_record(6,1),positions_record(6,2))-J_hsc1351;
J_hsc_2(positions_record(7,1),positions_record(7,2))=J_hsc_2(positions_record(7,1),positions_record(7,2))-J_hsc1351;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position = [];
for i = 1:size(J_hsc_2, 1)
    for j = 1:size(J_hsc_2, 2)
        if J_hsc_2(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_hsc_2, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position

% Loop through each element in the specified row
for j = 1:size(J_hsc_2, 2)  % Iterate through columns of the specified row
    if J_hsc_2(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)

    % Continue looping through the found_position(2) row
    for j = 1:size(J_hsc_2, 2)
        if J_hsc_2(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end

J_hsc1421=min([J_hsc_2(positions_record(8,1),positions_record(8,2)) J_hsc_2(positions_record(9,1),positions_record(9,2)) J_hsc_2(positions_record(10,1),positions_record(10,2))]);
J_hsc_3=J_hsc_2;
J_hsc_3(positions_record(8,1),positions_record(8,2))=J_hsc_3(positions_record(8,1),positions_record(8,2))-J_hsc1421;
J_hsc_3(positions_record(9,1),positions_record(9,2))=J_hsc_3(positions_record(9,1),positions_record(9,2))-J_hsc1421;
J_hsc_3(positions_record(10,1),positions_record(10,2))=J_hsc_3(positions_record(10,1),positions_record(10,2))-J_hsc1421;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position = [];
for i = 1:size(J_hsc_3, 1)
    for j = 1:size(J_hsc_3, 2)
        if J_hsc_3(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_hsc_3, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position

% Loop through each element in the specified row
for j = 1:size(J_hsc_3, 2)  % Iterate through columns of the specified row
    if J_hsc_3(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)

    % Continue looping through the found_position(2) row
    for j = 1:size(J_hsc_3, 2)
        if J_hsc_3(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end

J_hsc142351=min([J_hsc_3(positions_record(11,1),positions_record(11,2)) J_hsc_3(positions_record(12,1),positions_record(12,2)) J_hsc_3(positions_record(13,1),positions_record(13,2)) J_hsc_3(positions_record(14,1),positions_record(14,2)) J_hsc_3(positions_record(15,1),positions_record(15,2))]);
J_hsc_4=J_hsc_3;
J_hsc_4(positions_record(11,1),positions_record(11,2))=J_hsc_4(positions_record(11,1),positions_record(11,2))-J_hsc142351;
J_hsc_4(positions_record(12,1),positions_record(12,2))=J_hsc_4(positions_record(12,1),positions_record(12,2))-J_hsc142351;
J_hsc_4(positions_record(13,1),positions_record(13,2))=J_hsc_4(positions_record(13,1),positions_record(13,2))-J_hsc142351;
J_hsc_4(positions_record(14,1),positions_record(14,2))=J_hsc_4(positions_record(14,1),positions_record(14,2))-J_hsc142351;
J_hsc_4(positions_record(15,1),positions_record(15,2))=J_hsc_4(positions_record(15,1),positions_record(15,2))-J_hsc142351;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position = [];
for i = 1:size(J_hsc_4, 1)
    for j = 1:size(J_hsc_4, 2)
        if J_hsc_4(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_hsc_4, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position

% Loop through each element in the specified row
for j = 1:size(J_hsc_4, 2)  % Iterate through columns of the specified row
    if J_hsc_4(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)

    % Continue looping through the found_position(2) row
    for j = 1:size(J_hsc_4, 2)
        if J_hsc_4(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end

J_hsc142361=min([J_hsc_4(positions_record(16,1),positions_record(16,2)) J_hsc_4(positions_record(17,1),positions_record(17,2)) J_hsc_4(positions_record(18,1),positions_record(18,2)) J_hsc_4(positions_record(19,1),positions_record(19,2)) J_hsc_4(positions_record(20,1),positions_record(20,2))]);
J_hsc_5=J_hsc_4;
J_hsc_5(positions_record(16,1),positions_record(16,2))=J_hsc_5(positions_record(16,1),positions_record(16,2))-J_hsc142361;
J_hsc_5(positions_record(17,1),positions_record(17,2))=J_hsc_5(positions_record(17,1),positions_record(17,2))-J_hsc142361;
J_hsc_5(positions_record(18,1),positions_record(18,2))=J_hsc_5(positions_record(18,1),positions_record(18,2))-J_hsc142361;
J_hsc_5(positions_record(19,1),positions_record(19,2))=J_hsc_5(positions_record(19,1),positions_record(19,2))-J_hsc142361;
J_hsc_5(positions_record(20,1),positions_record(20,2))=J_hsc_5(positions_record(20,1),positions_record(20,2))-J_hsc142361;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J_hsc_5(J_hsc_5<10e-15)=0;
position = [];
for i = 1:size(J_hsc_5, 1)
    for j = 1:size(J_hsc_5, 2)
        if J_hsc_5(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_hsc_5, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position

% Loop through each element in the specified row
for j = 1:size(J_hsc_5, 2)  % Iterate through columns of the specified row
    if J_hsc_5(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)

    % Continue looping through the found_position(2) row
    for j = 1:size(J_hsc_5, 2)
        if J_hsc_5(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end

J_hsc14251=min([J_hsc_5(positions_record(21,1),positions_record(21,2)) J_hsc_5(positions_record(22,1),positions_record(22,2)) J_hsc_5(positions_record(23,1),positions_record(23,2)) J_hsc_5(positions_record(24,1),positions_record(24,2))]);
J_hsc_6=J_hsc_5;
J_hsc_6(positions_record(21,1),positions_record(21,2))=J_hsc_6(positions_record(21,1),positions_record(21,2))-J_hsc14251;
J_hsc_6(positions_record(22,1),positions_record(22,2))=J_hsc_6(positions_record(22,1),positions_record(22,2))-J_hsc14251;
J_hsc_6(positions_record(23,1),positions_record(23,2))=J_hsc_6(positions_record(23,1),positions_record(23,2))-J_hsc14251;
J_hsc_6(positions_record(24,1),positions_record(24,2))=J_hsc_6(positions_record(24,1),positions_record(24,2))-J_hsc14251;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position = [];
for i = 1:size(J_hsc_6, 1)
    for j = 1:size(J_hsc_6, 2)
        if J_hsc_6(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_hsc_6, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position

% Loop through each element in the specified row
for j = 1:size(J_hsc_6, 2)  % Iterate through columns of the specified row
    if J_hsc_6(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while ~isempty(found_position) && (found_position(2) ~= position(1) && found_position(2) ~= position(2))

    % Continue looping through the found_position(2) row
    for j = 1:size(J_hsc_6, 2)
        if J_hsc_6(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % 如果满足条件，则更新 position(2)
    if isempty(found_position) || found_position(2) == position(1) || found_position(2) == position(2)
        break;  
    end    
end

J_hsc4254=min([J_hsc_6(positions_record(26,1),positions_record(26,2)) J_hsc_6(positions_record(27,1),positions_record(27,2)) J_hsc_6(positions_record(28,1),positions_record(28,2))]);
J_hsc_7=J_hsc_6;
J_hsc_7(positions_record(26,1),positions_record(26,2))=J_hsc_7(positions_record(26,1),positions_record(26,2))-J_hsc4254;
J_hsc_7(positions_record(27,1),positions_record(27,2))=J_hsc_7(positions_record(27,1),positions_record(27,2))-J_hsc4254;
J_hsc_7(positions_record(28,1),positions_record(28,2))=J_hsc_7(positions_record(28,1),positions_record(28,2))-J_hsc4254;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
position = [];
for i = 1:size(J_hsc_7, 1)
    for j = 1:size(J_hsc_7, 2)
        if J_hsc_7(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_hsc_7, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position

% Loop through each element in the specified row
for j = 1:size(J_hsc_7, 2)  % Iterate through columns of the specified row
    if J_hsc_7(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)
    % Ensure found_position(2) is within matrix bounds
    if isempty(found_position) || found_position(2) > size(J_hsc_7, 1)
        break;  % Exit the loop if not found or out of bounds
    end

    % Continue looping through the found_position(2) row
    for j = 1:size(J_hsc_7, 2)
        if J_hsc_7(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end

J_hsc14261=min([J_hsc_7(positions_record(29,1),positions_record(29,2)) J_hsc_7(positions_record(30,1),positions_record(30,2)) J_hsc_7(positions_record(31,1),positions_record(31,2)) J_hsc_7(positions_record(32,1),positions_record(32,2))]);
J_hsc_8=J_hsc_7;
J_hsc_8(positions_record(29,1),positions_record(29,2))=J_hsc_8(positions_record(29,1),positions_record(29,2))-J_hsc14261;
J_hsc_8(positions_record(30,1),positions_record(30,2))=J_hsc_8(positions_record(30,1),positions_record(30,2))-J_hsc14261;
J_hsc_8(positions_record(31,1),positions_record(31,2))=J_hsc_8(positions_record(31,1),positions_record(31,2))-J_hsc14261;
J_hsc_8(positions_record(32,1),positions_record(32,2))=J_hsc_8(positions_record(32,1),positions_record(32,2))-J_hsc14261;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J_hsc_8(J_hsc_8<10e-15)=0;
position = [];
for i = 1:size(J_hsc_8, 1)
    for j = 1:size(J_hsc_8, 2)
        if J_hsc_8(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_hsc_8, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position

% Loop through each element in the specified row
for j = 1:size(J_hsc_8, 2)  % Iterate through columns of the specified row
    if J_hsc_8(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)
    % Ensure found_position(2) is within matrix bounds
    if isempty(found_position) || found_position(2) > size(J_hsc_8, 1)
        break;  % Exit the loop if not found or out of bounds
    end

    % Continue looping through the found_position(2) row
    for j = 1:size(J_hsc_8, 2)
        if J_hsc_8(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end

J_hsc1461=min([J_hsc_8(positions_record(33,1),positions_record(33,2)) J_hsc_8(positions_record(34,1),positions_record(34,2)) J_hsc_8(positions_record(35,1),positions_record(35,2))]);
J_hsc_9=J_hsc_8;
J_hsc_9(positions_record(33,1),positions_record(33,2))=J_hsc_9(positions_record(33,1),positions_record(33,2))-J_hsc1461;
J_hsc_9(positions_record(34,1),positions_record(34,2))=J_hsc_9(positions_record(34,1),positions_record(34,2))-J_hsc1461;
J_hsc_9(positions_record(35,1),positions_record(35,2))=J_hsc_9(positions_record(35,1),positions_record(35,2))-J_hsc1461;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J_hsc_9(J_hsc_9<10e-15)=0;
position = [];
for i = 1:size(J_hsc_9, 1)
    for j = 1:size(J_hsc_9, 2)
        if J_hsc_9(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if row index is within matrix bounds
if position(2) > size(J_hsc_9, 1) || position(2) < 1
    error('Specified row index is out of matrix bounds');
end

% Initialize variables
found_position = [];  % To store the found position

% Loop through each element in the specified row
for j = 1:size(J_hsc_9, 2)  % Iterate through columns of the specified row
    if J_hsc_9(position(2), j) > 0  % Check if the element is greater than 0
        found_position = [position(2), j];  % Record the position
        break;  % Exit the loop
    end
end

% Record the initial position
positions_record = [positions_record; position; found_position];

% Check if found_position(2) is equal to position(1)
while isempty(found_position) || found_position(2) ~= position(1)
    % Ensure found_position(2) is within matrix bounds
    if isempty(found_position) || found_position(2) > size(J_hsc_9, 1)
        break;  % Exit the loop if not found or out of bounds
    end

    % Continue looping through the found_position(2) row
    for j = 1:size(J_hsc_9, 2)
        if J_hsc_9(found_position(2), j) > 0  % Check if the element is greater than 0
            found_position = [found_position(2), j];  % Record the position
            break;  % Exit the loop
        end
    end
    
    % Record the new position
    positions_record = [positions_record; found_position];

    % If a new position is found, check found_position(2)
    if ~isempty(found_position) && found_position(2) ~= position(1)
        % Update position(2) to the value of found_position(2)
        position(2) = found_position(2);
    else
        break;  % Exit the loop if equal or not found
    end
end

J_hsc4654=min([J_hsc_9(positions_record(36,1),positions_record(36,2)) J_hsc_9(positions_record(37,1),positions_record(37,2)) J_hsc_9(positions_record(38,1),positions_record(38,2))]);
J_hsc_10=J_hsc_9;
J_hsc_10(positions_record(36,1),positions_record(36,2))=J_hsc_10(positions_record(36,1),positions_record(36,2))-J_hsc4654;
J_hsc_10(positions_record(37,1),positions_record(37,2))=J_hsc_10(positions_record(37,1),positions_record(37,2))-J_hsc4654;
J_hsc_10(positions_record(38,1),positions_record(38,2))=J_hsc_10(positions_record(38,1),positions_record(38,2))-J_hsc4654;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J_hsc_10(J_hsc_10<10e-15)=0;
position = [];
for i = 1:size(J_hsc_10, 1)
    for j = 1:size(J_hsc_10, 2)
        if J_hsc_10(i,j) > 0
            position = [i, j];
            break;
        end
    end
    if ~isempty(position)
        break;
    end
end

% Check if position is empty
if isempty(position)
    error('Decomposition is empty');
end
