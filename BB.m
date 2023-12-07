
%function [P,fes] = BB(N,ant_position,Elite_antlion_position,dim,fobj,fes)
function P = BB(N,ant_position,Elite_antlion_position,dim,fobj)
%%骨架变异
%N = n;
    for i=1:N
        Objective_values(1,i)=fobj(ant_position(i,:));
%           fes = fes + 1;
%           CR = 1-exp(-abs(Objective_values(1,i)-Destination_fitness));%%变异控制变量
        CR = 0.8;
        [ k1,k2,k3 ] = GetRan3( i,N );
        k=rand();
        for j=1:dim
            if rand()<CR
	      % 正态函数的方式生成特征值。
                mu=(Elite_antlion_position(j)+ant_position(i,j))/2;
                sigma=abs(Elite_antlion_position(j)-ant_position(i,j));
                V(i,j)=normrnd(mu,sigma);
            else
                %%DE/best/1
                % 选择三个个体的特征值进当前个体的特征值更新。
                V(i,j)=ant_position(k1,j)+k.*(ant_position(k2,j)-ant_position(k3,j));
            end
        end
        %fes = fes + 1;
        Vfitness(i)=fobj(V(i,:));
        if Vfitness(i)<Objective_values(1,i)
            ant_position(i,:)=V(i,:);
        end
    end
    %%骨架变异结束
    
    P = ant_position;
end