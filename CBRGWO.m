

% Grey Wolf Optimizer
function [Alpha_pos,Convergence_curve]=CBRGWO(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,dim);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,dim);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,dim);
Delta_score=inf; %change this to -inf for maximization problems

ML = 10;
gBestCollPos=zeros(ML,dim);
gBestCollCost=1.0e+100*ones(ML,1);
gBestCollCount=0;
gBestCollProb=zeros(ML,1);
MemoryLimit = ML;

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);
N=SearchAgents_no;
AllFitness = inf*ones(N,1);%record the fitness of all positions
it=1;
Convergence_curve=[];
FEs=0;

for i=1:size(Positions, 1)
    fitness=fobj(Positions(i,:));
    AllFitness(i)=fitness;
    % Update Alpha, Beta, and Delta
    if fitness<Alpha_score 
        Alpha_score=fitness; % Update alpha
        Alpha_pos=Positions(i,:);
        if gBestCollCount < ML && gBestCollCount>0
            gBestCollCount = gBestCollCount+1;
            gBestCollPos(gBestCollCount,:) = Positions(i,:);
            gBestCollCost(gBestCollCount) = AllFitness(i);
             %display(gBestCollCont);
        else
           % display(gBestCollCont);
            gBestCollCount = 1;
            gBestCollPos(gBestCollCount,:) = Positions(i,:);
            gBestCollCost(gBestCollCount) = AllFitness(i);
        end
    end

    if fitness>Alpha_score && fitness<Beta_score 
        Beta_score=fitness; % Update beta
        Beta_pos=Positions(i,:);
    end

    if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score 
        Delta_score=fitness; % Update delta
        Delta_pos=Positions(i,:);
    end
end

% Main loop
while  FEs < MaxFEs
% % Main loop
% while l<Max_iter
    for i=1:size(Positions,1)  
        
       % Return back the search agents that go beyond the boundaries of the search space
        %Positions=BB(N,Positions,Alpha_pos,dim,fobj);
        if FEs<MaxFEs/2
            % ------------------------------
            % 骨架机制
            %Elite_antlion_position = Alpha_pos.*0.5+Beta_pos.*0.3+Delta_pos.*0.2;
            Positions = BB(SearchAgents_no,Positions,Alpha_pos,dim,fobj);
        else
            % ------------------------------
            % 构建轮盘赌机制的轮盘
            F=zeros(ML,1);
            MeanCost = mean(gBestCollCost);
            for j=1:MemoryLimit
                F(j) = exp(-gBestCollCost(j)/MeanCost);
            end
            gBestCollProb(:)=F/sum(F);
            % 使用轮盘赌机制随机取一个值进行更新，可以跳出局部最优
            % 这个gbestX用作后面的Alpha_pos代替
            if range(gBestCollCost)~=0  % 不等于
                % 轮盘赌机制
                pr=RouletteWheelSelection(gBestCollProb);
                % Choose k randomly, not equal to i
                K=[1:pr-1 pr+1:MemoryLimit];
                k=K(randi([1 numel(K)]));
                gbestX=gBestCollPos(k,:);
            else
                % gBestCollCost==0说明没有更新过，Alpha_pos最好的
                gbestX = Alpha_pos;
            end
        end
        
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        if FEs<MaxFEs
            FEs=FEs+1;  
            % Calculate objective function for each search agent
            fitness=fobj(Positions(i,:));
            AllFitness(i)=fitness;
            % Update Alpha, Beta, and Delta
            if fitness<Alpha_score 
                Alpha_score=fitness; % Update alpha
                Alpha_pos=Positions(i,:);
                if gBestCollCount < ML && gBestCollCount>0
                    gBestCollCount = gBestCollCount+1;
                    gBestCollPos(gBestCollCount,:) = Positions(i,:);
                    gBestCollCost(gBestCollCount) = AllFitness(i);
                     %display(gBestCollCont);
                else
                   % display(gBestCollCont);
                    gBestCollCount = 1;
                    gBestCollPos(gBestCollCount,:) = Positions(i,:);
                    gBestCollCost(gBestCollCount) = AllFitness(i);
                end
            end

            if fitness>Alpha_score && fitness<Beta_score 
                Beta_score=fitness; % Update beta
                Beta_pos=Positions(i,:);
            end

            if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score 
                Delta_score=fitness; % Update delta
                Delta_pos=Positions(i,:);
            end
        else
            break;
        end
    end
     
    % ------------------------------
    
    %a=2-FEs*((2)/MaxFEs); % a decreases linearly fron 2 to 0
    %a = 1/(cos(FEs/MaxFEs)+0.001);
    %x = FEs*(2/MaxFEs);
    %a = 1/(x^3+x^2+1);
    n=1-FEs*((1)/MaxFEs);
    if FEs <= MaxFEs/2
        % 对a加了正余弦机制
        a = 0 + (2 - 0)*(1 + (cos(((FEs-1)*pi)/(MaxFEs-1)))^n)/2;
    else
        a = 0 + (2 - 0)*(1 - abs(cos(((FEs-1)*pi)/(MaxFEs-1)))^n)/2;
    end
    
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        V=zeros(1,dim);
        for j=1:size(Positions,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            if FEs<MaxFEs/2
                D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
                X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
            else
                D_alpha=abs(C1*gbestX(j)-Positions(i,j)); % Equation (3.5)-part 1
                X1=gbestX(j)-A1*D_alpha; % Equation (3.6)-part 1
            end
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2
            
            r1=rand();
            r2=rand();
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3
            
            %Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            % 加了动态权重
            w1 = abs(X1)/abs(X1+X2+X3);
            w2 = abs(X2)/abs(X1+X2+X3);
            w3 = abs(X3)/abs(X1+X2+X3);
            V(j)=(X1*w1+X2*w2+X3*w3)/3;% Equation (3.7)
            
        end
        %% CGO
        I = randi([1,2],1,12);
        Ir = randi([0,1],1,5);
        RandGroup = randperm(N,randperm(N,1));
        MG = mean(Positions(RandGroup,:)).*(length(RandGroup)~=1)+Positions(RandGroup(1,1),:)*(length(RandGroup)==1);

        Alpha(1,:) = rand(1,dim);
        Alpha(2,:) = 2*rand(1,dim)-1;
        Alpha(3,:) = Ir(1)*rand(1,dim)+1;
        Alpha(4,:) = Ir(2)*rand(1,dim)+(~Ir(2));
        ii=randi([1,4],1,3);
        SelectedAlpha = Alpha(ii,:);

        Seed1 = Positions(i,:) + SelectedAlpha(1,:).*(I(1).*Alpha_pos-I(2).*MG);
        Seed2 = Alpha_pos + SelectedAlpha(2,:).*(I(3).*MG-I(4).*Positions(i,:));
        Seed3 = MG + SelectedAlpha(3,:).*(I(5).*Positions(i,:)-I(6).*Alpha_pos);

        F_Seed1 = fobj(Seed1);
        F_Seed2 = fobj(Seed2);
        F_Seed3 = fobj(Seed3);
        F_V = fobj(V);
        FEs=FEs+4;

        Fitness_comb=[F_Seed1,F_Seed2,F_Seed3,F_V,AllFitness(i)];
        [~,m]=min(Fitness_comb);
        
        curFitness = Fitness_comb(m);
        if curFitness<Alpha_score 
            Alpha_score=curFitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end

        if curFitness>Alpha_score && curFitness<Beta_score 
            Beta_score=curFitness; % Update beta
            Beta_pos=Positions(i,:);
        end

        if curFitness>Alpha_score && curFitness>Beta_score && curFitness<Delta_score 
            Delta_score=curFitness; % Update delta
            Delta_pos=Positions(i,:);
        end
        
        
        if m==1   
            Positions(i,:)=Seed1;
        elseif m==2
            Positions(i,:)=Seed2; 
        elseif m==3
            Positions(i,:)=Seed3; 
        elseif m==4
            Positions(i,:)=V; 
        end
        %%
    end
    Convergence_curve(it)=Alpha_score;
    it=it+1;
end

end

