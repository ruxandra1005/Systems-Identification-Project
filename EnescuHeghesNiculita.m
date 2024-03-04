clear all
close all
clc

% CITIREA DATELOR

load('proj_fit_40.mat');

X1_id = id.X{1};
X2_id = id.X{2};
Y_id = id.Y;

X1_val = val.X{1};
X2_val = val.X{2};
Y_val = val.Y;


% GENERAREA MATRICEI DE REGRESORI PT IDENTIFICARE
grad = 30;
mse_id = ones(1 , grad);
mse_val = ones(1, grad);

for m=1:grad

phi_id = zeros(length(X1_id)*length(X2_id), nchoosek(m+2 , m));
y_v = zeros(length(X1_id)^2 , 1); % vectorul matricei Y_id

j=1;
while (j<=length(X2_id)) 
    for i = 1:length(X1_id)^2 % for asta e pt parcurgere X1_id (am pus la patrat pentru ca nr de linii din phi este 1681)
        k = 1;
        for a = 0:m
            for b = 0:m-a
                if a+b<=m && a <= m && b <= m 
                    if  mod(i,41) == 0 
                        phi_id(i, k) = (X1_id(41)^a)*(X2_id(j)^b);
                        y_v(i , 1) = Y_id(41, j); 
                    else
                        phi_id(i, k) = (X1_id(mod(i,41))^a)*(X2_id(j)^b);  
                        y_v(i , 1) = Y_id(mod(i, 41) , j);
                    
                    end
                    
                end
                k = k+1;
            end
        end
        if mod(i, 41) == 0
            j=j+1;
        end
    end
end


% IDENTIFICARE COEFICIENTI THETA

theta = ones(nchoosek(m+2, m) , 1);
theta = phi_id\y_v ;

% GENERAREA MATRICEI DE REGRESORI PT VALIDARE

phi_val = zeros( length(X1_val)*length(X1_val), nchoosek(m+2 , m));

j=1;
while (j<=length(X2_val)) 
    for i = 1:length(X1_val)^2 % for asta e pt parcurgere X1_val (am pus la patrat pentru ca nr de linii din phi este 961)
        k = 1;
        for a = 0:m
            for b = 0:m-a
                if a+b<=m && a <= m && b <= m 
                    if  mod(i,31) == 0 
                        phi_val(i, k) = (X1_val(31)^a)*(X2_val(j)^b);
                    else
                        phi_val(i, k) = (X1_val(mod(i,31))^a)*(X2_val(j)^b); 
                    end
                end
                k = k+1;
            end
        end
        if mod(i, 31) == 0
            j=j+1;
        end
    end
end

% GENERARE APROXIMARE YHAT_VALIDARE

Yhat_val_v = ones(length(X1_val)^2 , 1);
Yhat_val = ones(length(X1_val) , length(X2_val));
Yhat_val_v = phi_val * theta;

j =1;
while(j<= length(X1_val))
for i = 1:length(X1_val)^2
    if  mod(i,31) == 0 
        Yhat_val(31 , j) = Yhat_val_v(i ,1);
        j= j+1;
    else
        Yhat_val(mod(i , 31) , j) = Yhat_val_v(i ,1);              
    end
         
end
end


% GENERARE APROXIMARE YHAT_IDENTIFICARE

Yhat_id_v = ones(length(X1_id)^2 , 1);
Yhat_id = ones(length(X1_id) , length(X2_id));
Yhat_id_v = phi_id * theta;

j =1;
while(j<= length(X1_id))
for i = 1:length(X1_id)^2
    if  mod(i,41) == 0 
        Yhat_id(41 , j) = Yhat_id_v(i ,1);
        j= j+1;
    else
        Yhat_id(mod(i , 41) , j) = Yhat_id_v(i ,1);              
    end
         
end
end

% MSE IDENTIFICARE

MSE_id=0;

for p=1:41
    for q=1:41
        MSE_id = MSE_id + (Y_id(p, q)-Yhat_id(p, q))*(Y_id(p, q)-Yhat_id(p, q));
    end
end

MSE_id = MSE_id/(41*41);

% MSE VALIDARE

MSE_val=0;

for p=1:31
    for q=1:31
        MSE_val = MSE_val + (Y_val(p, q)-Yhat_val(p, q))*(Y_val(p, q)-Yhat_val(p, q));
    end
end

MSE_val = MSE_val/(31*31);

mse_id(1, m) = MSE_id;
mse_val(1 , m) = MSE_val;

end

stem(mse_id, 'b');
hold on
stem(mse_val, 'r');
grid;
xlabel('m');
ylabel('MSE');
title ('EVOLUTIA MSE');
legend('MSE id' , 'MSE val');

% GRADELE MSE MINIME 

mse_min_id = min(mse_id);
for i = 1:m
    if mse_min_id == mse_id(1, i)
        grad_mse_min_id = i; 
    end
end

mse_min_val = min(mse_val);
for i = 1:m
    if mse_min_val == mse_val(1, i)
        grad_mse_min_val = i; 
    end
end

% PENTRU GRADUL MINIM AL VALIDARII => 

m = grad_mse_min_val;
phi_id = zeros( length(X1_id)*length(X1_id), nchoosek(m+2 , m));
y_v = zeros(length(X1_id)^2 , 1); % vectorul matricei Y_id

j=1;
while (j<=length(X2_id)) 
    for i = 1:length(X1_id)^2 % for asta e pt parcurgere X1 (am pus la patrat pentru ca nr de linii din phi este 1681)
        k = 1;
        for a = 0:m
            for b = 0:m-a
                if a+b<=m && a <= m && b <= m 
                    if  mod(i,41) == 0 
                        phi_id(i, k) = (X1_id(41)^a)*(X2_id(j)^b);
                        y_v(i , 1) = Y_id(41, j); 
                    else
                        phi_id(i, k) = (X1_id(mod(i,41))^a)*(X2_id(j)^b);  
                        y_v(i , 1) = Y_id(mod(i, 41) , j);
                    
                    end
                    
                end
                k = k+1;
            end
        end
        if mod(i, 41) == 0
            j=j+1;
        end
    end
end


% IDENTIFICARE COEFICIENTI THETA

theta = ones(nchoosek(m+2, m) , 1);
theta = phi_id\y_v ;

% GENERAREA MATRICEI DE REGRESORI PT VALIDARE

phi_val = zeros( length(X1_val)*length(X1_val) , nchoosek(m+2 , m));
j=1;
while (j<=length(X2_val)) 
    for i = 1:length(X1_val)^2 % for asta e pt parcurgere X1 (am pus la patrat pentru ca nr de linii din phi este 1681)
        k = 1;
        for a = 0:m
            for b = 0:m-a
                if a+b<=m && a <= m && b <= m 
                    if  mod(i,31) == 0 
                        phi_val(i, k) = (X1_val(31)^a)*(X2_val(j)^b);
                    else
                        phi_val(i, k) = (X1_val(mod(i,31))^a)*(X2_val(j)^b); 
                    end
                end
                k = k+1;
            end
        end
        if mod(i, 31) == 0
           j=j+1;
        end
    end
end

% GENERARE APROXIMARE YHAT_VALIDARE

Yhat_val_v = ones(length(X1_val)^2 , 1);
Yhat_val = ones(length(X1_val) , length(X2_val));
Yhat_val_v = phi_val * theta;

j =1;
while(j<= length(X1_val))
for i = 1:length(X1_val)^2
    if  mod(i,31) == 0 
        Yhat_val(31 , j) = Yhat_val_v(i ,1);
        j= j+1;
    else
        Yhat_val(mod(i , 31) , j) = Yhat_val_v(i ,1);              
    end
         
end
end


% GENERARE APROXIMARE YHAT_IDENTIFICARE

Yhat_id_v = ones(length(X1_id)^2 , 1);
Yhat_id = ones(length(X1_id) , length(X2_id));
Yhat_id_v = phi_id * theta;

j =1;
while (j<= length(X1_id))
for i = 1:length(X1_id)^2
    if  mod(i,41) == 0 
        Yhat_id(41 , j) = Yhat_id_v(i ,1);
        j= j+1;
    else
        Yhat_id(mod(i , 41) , j) = Yhat_id_v(i ,1);              
    end
         
end
end

% AFISARI APROXIMARI

figure;

subplot(121);
mesh(X1_val , X2_val, Yhat_val);
%hold on
%mesh(X1_val , X2_val, Y_val);

xlabel('X1 val');
ylabel('X2 val');
zlabel('Yhat val');
title('APROXIMARE VALIDARE');

subplot(122);
mesh(val.X{1}, val.X{2}, val.Y);

xlabel('X1 val');
ylabel('X2 val');
zlabel('Y val');
title('DATE VALIDARE');

figure;

subplot(121);
mesh(X1_id , X2_id, Yhat_id);
%hold on
%mesh(X1_id , X2_id, Y_id);

xlabel('X1 id');
ylabel('X2 id');
zlabel('Yhat id');
title('APROXIMARE IDENTIFICARE');

subplot(122);
mesh(id.X{1}, id.X{2}, id.Y);
xlabel('X1 id');
ylabel('X2 id');
zlabel('Y id');
title('DATE IDENTIFICARE');

