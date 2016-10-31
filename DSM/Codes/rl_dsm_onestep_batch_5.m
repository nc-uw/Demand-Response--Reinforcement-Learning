%sample of 24hr electricty cost in cents
run('tariff.m');
p_max=ceil(max(p_all(:)));

%function 1 -
    %create all possible on/off states over horizon period h for -
    %h running time app_rt = (1,2,..h) combinations for 1 appliance
h=24;
maxx=9999;
%create all possible on/off states over horizon period for all possible running time combinations
app_states_all=cell(h,1);
app_states_all{1}=eye(h,h);

for i = 2:h
    app_states_all{i}=eye(h-(i-1),h)+app_states_all{i-1}(2:end,:);
end;

%specify app data
app_data = [0.72 3 24; 3.15 2 24; 3.18 3 24];
%10.5 1 24; 5.5 3 24; 17 1 24; 0.72 3 24; 3.15 2 24; 3.18 3 24; 10.5 1 24];

N = size(app_data,1);
app_kwh = app_data(:,1);
app_rt = app_data(:,2);
app_reqt = app_data(:,3) - app_data(:,2) +1;
kwh_max=p_max*max(app_kwh);

for pi=1:1 %size(p_all,1)
    
    %p = transpose(p_all(pi,:));
    p=[35.33	31.36	32.27	32.35	30.80	33.87	43.19	48.24	43.47	42.13	39.22	37.35	34.77	33.20	31.39	31.54	35.84	47.29	45.17	39.98	35.65	34.07	34.32	32.66]';
    %extract possible on/off states for specified apps, proceed to derive
    %kwh consumption and tarrif states
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



    %BATCH SARSA - ALPHA BACKUP
    %for q learning initialize decay = 1 else decay=10/100/1000;
    %tweak
    alpha = 0.5;
    gamma = 1;
    epsilon=0.0;
    decay=1000;

    %initialization
    iter=1;
    check=0;
    pol_state = [];
    pol_hour = [];
    pol_iter = [];
    pol_tariff = [];

    while (check==0)

        if iter==1
            QS1=zeros(size(R));
            %QS1=R;
        else
            QS1=QS2;
        end;

        %epsilon = epsilon - (1/decay);
        rndm=rand;

        if epsilon < 0 
            epsilon= 0;
        end;   

        Qdash2 = [];
        for i = 1:h

            Qdash=[];
            for j = 0:2^N-1

                counter=j*h+1;
                a=counter+i;
                b=(j+1)*h;

                q = QS1(a:b,:); 

                if isempty(q(q~=-maxx))
                    pol_max_q=0;
                    pol_rndm_q=0;
                else
                    pol_max_q=max(q(:));
                    pol_rndm_q=datasample(q(q~=-maxx),1);
                end;    

                if rndm < epsilon
                    Qdash = [Qdash pol_rndm_q];
                else
                    Qdash = [Qdash pol_max_q];
                end;
            end;

        Qdash2 = [Qdash2; Qdash];
        end;

    %Qdash2(Qdash2==-maxx) = 0;
    %Qdash2 = reshape(transpose(Qdash2), [h 2^N]);
    STQdash=[];
    for i = 1:2^N
        STQdash = [STQdash ; Qdash2.*repmat(ST(i,:),24,1)];
    end;

    TD=alpha*((R+gamma*STQdash)-QS1);
    TD(TD==(-maxx*alpha)) = -maxx;
    QS2 = QS1 + TD;
    QS2(QS2<(-maxx)) = -maxx;
    diff = abs(QS2-QS1);
    check = all(diff(:) <=0.1);

    %greedy policy
    %intialize
    pol_row=0;
    pol_col=1;
    deadlock_size = [];
    opt_state_iter = [];
    opt_hour_iter = [];
    opt_tariff_iter = [];
    n=0;
    n2=1;

        while (n < N)

            a=(pol_col-1)*h + (pol_row+1);
            b=pol_col*h; 
            q=QS2(a:b,:);
            [pol_row, pol_col] = ind2sub(size(q),find(q==max(q(:))));
            deadlock=randi(size(pol_row,1));
                    if size(pol_row,1) > 1;
                    deadlock_size = [deadlock_size; size(pol_row,1) max(q(:))];
                    end;
                 pol_row=pol_row(deadlock,1);
                 pol_col=pol_col(deadlock,1);
            pol_row = pol_row+h-1-(b-a);
            opt_state_iter = [opt_state_iter; pol_col];
            opt_hour_iter = [opt_hour_iter; pol_row];
            if n2==1
                opt_tariff_iter(n2,:)=app_states_kwh(opt_hour_iter(n2),:)*S(opt_state_iter(n2),:)';  
            else
                opt_tariff_iter(n2,:)=app_states_kwh(opt_hour_iter(n2),:)*(S(opt_state_iter(n2),:)-S(opt_state_iter(n2-1),:))';
            end;

            n=sum(S(pol_col,:));
            n2=n2+1;
        end;

        opt_state_iter = [opt_state_iter; zeros(N-size(opt_state_iter,1),1)];
        opt_hour_iter = [opt_hour_iter; zeros(N-size(opt_hour_iter,1),1)];
        opt_tariff_iter = [opt_tariff_iter; zeros(N-size(opt_tariff_iter,1),1)];

        if iter==1
            pol_state = [opt_state_iter];
            pol_hour = [opt_hour_iter];
            pol_iter = [iter];
            pol_tariff = [opt_tariff_iter];

        else
            pol_state = [pol_state opt_state_iter];
            pol_hour = [pol_hour opt_hour_iter];
            pol_iter = [pol_iter iter];
            pol_tariff = [pol_tariff opt_tariff_iter];
        end;

    iter=iter+1;

    end;

    pol_total_tariff = sum(pol_tariff);
    policy_evol= [pol_iter; pol_state; pol_hour; pol_total_tariff];

    pi;
    opt_state_iter;
    opt_hour_iter;
    iter;

    for i = 1:find(opt_state_iter==2^N)
        if i==1              
            total_tariff=[];
            app_number=[];
            for j = 1:N
                if S(opt_state_iter(i,:),j)==1
                    total_tariff=[total_tariff; app_states{j}(opt_hour_iter(i,:),:).*S(opt_state_iter(i,:),j)];
                    app_number = [app_number j];
                end;
            end;
        else
            for j = 1:N
                if (S(opt_state_iter(i,:),j)-S(opt_state_iter(i-1,:),j))==1
                    total_tariff=[total_tariff; app_states{j}(opt_hour_iter(i,:),:).*(S(opt_state_iter(i,:),j)-S(opt_state_iter(i-1,:),j))];
                    app_number = [app_number j];
                end;
            end;
        end;
    end;
    total_tariff=[app_number' total_tariff];
    total_tariff=sortrows(total_tariff,1);
    total_tariff=total_tariff(:,2:end);
    x=1:N;
    y=x;
   
    %schedule- plot double y
    figure
    stackedbar = @(x, A) bar(x, A, 0.75, 'stack');
    prettyline = @(x, y) plot(x, y, 'k', 'LineWidth', 1);
    [ax, h1, h2] = plotyy(1:h, total_tariff', 1:h, p, stackedbar, prettyline);
    set(ax, 'XLim', [0 h+1]);
    set(ax, 'XTick', 1:h+1);
    set(ax(1), 'YLim', [0 60]);
    set(ax(1), 'YTick', 0:10:60);
    set(ax(2), 'YLim', [0 p_max]);
    set(ax(2), 'YTick', 0:10:p_max);
    colormap gray;
    title('Q-learning - Appliance Schedule at 30 iterations');
    xlabel('Hour of Day');
    ylabel(ax(1),'Estimated Tariff (in Cents)');
    ylabel(ax(2),'Price per Hour (in Cents)');
    legend('Appliance 1', 'Appliance 2', 'Appliance 3', 'Appliance 4', 'Appliance 5', 'Appliance 6', 'Tariff per KWh: 1st Jan 2016', 'Location', 'northeast');
pi    
end;

figure
plot(1:iter-1,sum(pol_tariff(:,1:iter-1)));


% for ii=1:3
%     saveas(ii,[filename '-' num2str(ii)],'pdf)
% end

        

        





