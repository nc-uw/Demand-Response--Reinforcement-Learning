%sample of 24hr electricty cost in cents
run('tariff.m');
pi=randi(size(p_all,1),1);
p = transpose(p_all(pi,:));
p=[33.55	28.47	24.04	23.98	25.42	28.61	29.84	32.19	29.06	34.28	35.82	38.72	38.06	36.88	36.45	33.60	42.05	63.66	48.44	46.87	46.93	42.50	35.08	31.06]';
p=[35.33	31.36	32.27	32.35	30.80	33.87	43.19	48.24	43.47	42.13	39.22	37.35	34.77	33.20	31.39	31.54	35.84	47.29	45.17	39.98	35.65	34.07	34.32	32.66]';
%function 1 -
    %create all possible on/off states over horizon period h for -
    %h running time app_rt = (1,2,..h) combinations for 1 appliance
h=24;
maxx=9999;
app_states_all=cell(h,1);
app_states_all{1}=eye(h,h);
for i = 2:h
    app_states_all{i}=eye(h-(i-1),h)+app_states_all{i-1}(2:end,:);
end;

%function 2 -
    %specify appliance data (in code / user prompt)
%input applaince data 
app_data = [0.72 3 24; 3.15 2 24; 3.18 3 24; 10.5 1 24; 5.5 3 24; 17 1 24];
%retrieve appliance details
N = size(app_data,1);
app_kwh = app_data(:,1);
app_rt = app_data(:,2);
app_reqt = app_data(:,3) - app_data(:,2) +1;

%function 3 -
    %use outputs of function 1 and function 2 to create -  
    %1. possible? on/off states for specified applainces, 
    %2. kwh consumption per hour matrice
    %3. tariff matrice
app_states = cell(N,1);
app_states_kwh = [];
for i = 1:N
    app_states{i}=app_kwh(i).*app_states_all{app_rt(i)};
    temp=app_states{i}*p;
    app_states_kwh(:,i)=[temp ; zeros(h-size(temp,1),1)];
end;
%app_states_kwh(app_states_kwh==0) = maxx;

s = [0 1];
S = permn(s,N);
ST=[];

for i = 1:2^N
    for j = 1:2^N;
        x=find(S(i,:)==1);
        y=find(S(j,:)==1);
        if and(size(x,2)<size(y,2), ismember(x,y) == 1);
            ST(i,j)=1;
        else
            ST(i,j)=0;
        end;
    end;
end;
ST(1,:)=ones(1,2^N);
ST(:,1)=zeros(2^N,1);
C=[];

c=app_states_kwh*S';
for i = 1:2^N-1
    min_on=min(app_reqt(((S(i,:)==zeros).*app_reqt')~=0));
    CS=(c-repmat(c(:,i),1,2^N)).*repmat(ST(i,:),24,1);
    CS=[CS(1:min_on,:); zeros(h-min_on,2^N)];
    CS(min_on,1:2^N-1)=0;
    C = [C ; CS];
end;
    
C=[C; zeros(h,2^N)];
C(C==0) = maxx;
%C(C>=maxx) = maxx;
R=-1.*C;

%EPISODE Q-LEARNING / SARSA - ALPHA BACKUP

%for Q-LEARNING initialize: 
%epsilon=0
%decay = 1 else 

%for SARSA initialize: 
%0<epsilon<=1
%or decay=10/100/1000
alpha = 1;
gamma = 1;
decay=1000;

%initialization
iter=1;
check=0;
pol_state = [];
pol_hour = [];
pol_iter = [];
pol_tariff = [];
pol_a=[];
pol_b=[];
epsilon=0.0;
q_flg=1;
%while (check==0)
iter_total=15000;
for iter=1:iter_total

    if iter==1
        %QS1=-maxx*(R==-maxx);
        %QS2=-maxx*(R==-maxx);
        %QS1=zeros(size(R));
        %QS2=zeros(size(R));
        QS1=R;
        QS2=R;
    else
        QS1=QS2;
    end;  

    %episodic learning - intialize
    Qdash2=[];
    Qrnd=[];
    Qmax=[];
    Si=[];
    Hi=[];
    Ri=[]; 
    deadlock_size_rnd = [];  
    deadlock_size_max = [];   
    n=0;
    
    %episodic learning - start episode    
    opt_state_iter = [];
    opt_hour_iter = [];
    opt_tariff_iter = [];
    opt_a_iter=[];
    opt_b_iter=[];
    while (n < N)                                  
        
        rndm=rand;        
        %epsilon = epsilon - (1/decay);
        if epsilon < 0 
            epsilon= 0;
        end;  
        
        if n==0
            a=1;
            b=24;
            q=QS1(1:24,:);
            r=R(a:b,:);
            
            if rndm < epsilon
                pol_row_1=1;
                pol_col_1=1;
                while R(pol_row_1,pol_col_1)==-maxx
                    pol_row_1=randi(size(q,1));
                    pol_col_1=randi(size(q,2));
                end;    
            else           
                %max selection
                pol_max_q=max(q(r~=-maxx));
                [pol_row_1, pol_col_1]=find(q==pol_max_q);

                %random row, column tie breaker
                if size(pol_col_1,1) > 1;
                    deadlock_max=randi(size(pol_col_1,1));
                    deadlock_size_max = [deadlock_size_max; size(pol_col_1,1) pol_max_q];
                    pol_col_1 = pol_col_1(deadlock_max,1);
                    pol_row_1 = pol_row_1(deadlock_max,1);     
                end;                 
            end;
            %colons
            opt_state_iter=[opt_state_iter; pol_col_1];
            opt_hour_iter=[opt_hour_iter; pol_row_1];
            opt_tariff_iter=[opt_tariff_iter; R(pol_row_1,pol_col_1)];
            opt_a_iter=[opt_a_iter;a];
            opt_b_iter=[opt_b_iter;b];    
            
            pol_col_a=pol_col_1;
            pol_row_a=pol_row_1;
        else
            
            a=(pol_col_a-1)*h + (pol_row_a+1);
            b=pol_col_a*h;
            q=QS1(a:b,:);
            r=R(a:b,:);    

            if isempty(q(r~=-maxx))
                pol_max_q=0
                pol_rndm_q=0  
            else            
                %random selection            
                pol_rnd_q=datasample(q(r~=-maxx),1);
                [pol_rnd_row, pol_rnd_col]=find(q==pol_rnd_q);

                %random row, column tie breaker
                if size(pol_rnd_col,1) > 1;
                    deadlock_rnd=randi(size(pol_rnd_col,1));
                    deadlock_size_rnd = [deadlock_size_rnd; size(pol_rnd_col,1) pol_rnd_q];
                    pol_rnd_col = pol_rnd_col(deadlock_rnd,1);
                    pol_rnd_row = pol_rnd_row(deadlock_rnd,1)+pol_row_a;
                else
                    pol_rnd_row=pol_rnd_row+pol_row_a;
                end;    

                %max selection
                pol_max_q=max(q(r~=-maxx));
                [pol_max_row, pol_max_col]=find(q==pol_max_q);

                %random row, column tie breaker
                if size(pol_max_col,1) > 1;
                    deadlock_max=randi(size(pol_max_col,1));
                    deadlock_size_max = [deadlock_size_max; size(pol_max_col,1) pol_max_q];
                    pol_max_col = pol_max_col(deadlock_max,1);
                    pol_max_row = pol_max_row(deadlock_max,1)+pol_row_a;
                else
                    pol_max_row=pol_max_row+pol_row_a;
                end; 
                
                if rndm < epsilon
                    if q_flg==0                
                        TD = alpha*((R(pol_row_1,pol_col_1) + gamma*(QS1(pol_rnd_row+(pol_col_1-1)*h,pol_rnd_col))) - QS1(pol_row_1,pol_col_1));                                             
                    else         
                        TD = alpha*((R(pol_row_1,pol_col_1) + gamma*(QS1(pol_max_row+(pol_col_1-1)*h,pol_max_col))) - QS1(pol_row_1,pol_col_1));                         
                    end;                   
                    QS2(pol_row_1,pol_col_1) = QS1(pol_row_1,pol_col_1) + TD;                                        
                    pol_row_2 = pol_rnd_row;
                    pol_col_2 = pol_rnd_col;
                else
%                         pol_col_1
%                         pol_row_1
%                         R(pol_row_1,pol_col_1)
%                         QS1(pol_row_1,pol_col_1)
%                         pol_max_col
%                         pol_max_row
%                         pol_max_row+(pol_col_1-1)*h
%                         pol_max_q
%                         QS1(pol_max_row+(pol_col_1-1)*h,pol_max_col)
                    TD = alpha*((R(pol_row_1,pol_col_1) + gamma*(QS1(pol_max_row+(pol_col_1-1)*h,pol_max_col))) - QS1(pol_row_1,pol_col_1));
                    QS2(pol_row_1,pol_col_1) = QS1(pol_row_1,pol_col_1) + TD;
%                         TD
%                         QS2(pol_row_1,pol_col_1)
                    pol_row_2 = pol_max_row;
                    pol_col_2 = pol_max_col;   
                end; 
            end;  
            %colons
            opt_state_iter=[opt_state_iter; pol_col_2];
            opt_hour_iter=[opt_hour_iter; pol_row_2];
            opt_tariff_iter=[opt_tariff_iter; R(pol_row_2+(pol_col_1-1)*24,pol_col_2)];
            opt_a_iter=[opt_a_iter;a];
            opt_b_iter=[opt_b_iter;b];
            
            pol_row_1 = pol_row_2+(pol_col_1-1)*h;
            pol_col_1 = pol_col_2;
            
            pol_row_a = pol_row_2;
            pol_col_a = pol_col_2;
        end;
        
        n=sum(S(pol_col_a,:));
        Si = [Si; S(pol_col_a,:)];
        Hi = [Hi; pol_row_a];
               
    end;     
    %episodic learning - end episode          
    
    diff = abs(QS2-QS1);
    check = all(diff(:) <=0.0001);
    
    iter
    opt_state_iter = [opt_state_iter; zeros(N-size(opt_state_iter,1),1)];
    opt_hour_iter = [opt_hour_iter; zeros(N-size(opt_hour_iter,1),1)];
    opt_tariff_iter = [opt_tariff_iter; zeros(N-size(opt_tariff_iter,1),1)];
    opt_a_iter=[opt_a_iter;zeros(N-size(opt_a_iter,1),1)];
    opt_b_iter=[opt_b_iter;zeros(N-size(opt_b_iter,1),1)];

    if iter==1
        pol_state = [opt_state_iter];
        pol_hour = [opt_hour_iter];
        pol_iter = [iter];
        pol_tariff = [opt_tariff_iter];
        pol_a=[opt_a_iter];
        pol_b=[opt_b_iter];

    else
        pol_state = [pol_state opt_state_iter];
        pol_hour = [pol_hour opt_hour_iter];
        pol_iter = [pol_iter iter];
        pol_tariff = [pol_tariff opt_tariff_iter];
        pol_a=[pol_a opt_a_iter];
        pol_b=[pol_b opt_b_iter];
    end;
   
    %iter=iter+1;
end;

opt_state_iter
opt_hour_iter
opt_tariff_iter
QS2(1:24,:)
R(1:24,:)

figure            
plot(1:iter-1,-1*sum(pol_tariff(:,1:iter-1)));

