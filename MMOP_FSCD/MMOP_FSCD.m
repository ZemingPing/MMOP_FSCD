function [ ps,pf] = MMOP_FSCD( func_name,VRmin,VRmax,n_obj,NP,Max_Gen )

n_var=size(VRmin,2);    %Obtain the dimensions of decision space
Max_FES=Max_Gen*NP;     %Maximum fitness evaluations
%% Initialize population
VRmin=repmat(VRmin,NP,1);
VRmax=repmat(VRmax,NP,1);
pos=VRmin+(VRmax-VRmin).*rand(NP,n_var); %initialize the positions of the individuals
%% Evaluate the population
fitness=zeros(NP,n_obj);
for i=1:NP
    fitness(i,:)=feval(func_name,pos(i,:));
end
fitcount=NP;            % count the number of fitness evaluations
G=0;

%% main loop
while fitcount <= Max_FES
    G=G+1;
    %设置参数控制差分算子
    P_temp = 1-((G-1)/Max_Gen);

    Parentstar=[pos,fitness];
    K=10;
    Parentstar = non_domination_scd_kmeans_sort(Parentstar,n_obj,n_var,K);% n_var+n_obj+1: the level of front   +2：the decision space distance +3：CSCD value  +4：the objective space distance
    
    % calculate the number of individals in each front
    Total_pf = max(Parentstar(:,n_var+n_obj+1));
    pop = cell(Total_pf,1);
    num= zeros(Total_pf,1);
    
    % divide the population into multiple cells based on their fornt
    for i=1:NP
        pop{Parentstar(i,n_obj+n_var+1)} = [ pop{Parentstar(i,n_obj+n_var+1)};Parentstar(i,:)];
        num(Parentstar(i,n_obj+n_var+1)) = num(Parentstar(i,n_obj+n_var+1)) + 1;
    end
    
    % sort the individuals in each front based on CSCD value
    for i=1:Total_pf
        [~,index_sorted_on_scd]= sort(pop{i,1}(:,n_var+n_obj+3),'descend');  
        pop{i,1}= pop{i,1}(index_sorted_on_scd,:);
    end
    
    offspring = [];
    if rand < P_temp       
       for i=1: Total_pf
           u = zeros(num(i),n_var);
           for j= 1 : size(pop{i,1},1)
               % select exemplar based on the proposed DBESM mechanism
               if i==1
                   r_pf(1:3) = 1;
               else
                   r_pf = ceil(rand(1,3)*(i-1));
               end

               if  num(r_pf(1)) < 3
                   r1 = ceil(rand* num(r_pf(1))); 
               else
                   r1 = min(ceil(0.1 * num(r_pf(1))),3);  % p = 0.1;
               end
              
               distance = sum((repmat(pop{i,1}(j,1:n_var), r1 ,1) - pop{r_pf(1),1}(1:r1,1:n_var)).^2,2);
               if length(distance) == 1
                   select_index = r1;
               else
                   if length(find(max(distance)==distance))== r1
                       select_index = ceil(rand*r1);
                   else
                       distance = max(distance) - distance;
                       select_pro =  distance./(sum(distance));
                       true_pro = cumsum(select_pro);
                       select_index = find(rand<true_pro);
                   end
               end
                
               r2 = ceil(rand*num(r_pf(2)));
               r3 = ceil(rand*num(r_pf(3)));
               % Mutation
               v = pop{i,1}(j,1:n_var) + 0.8*(pop{r_pf(1),1}(select_index(1),1:n_var) - pop{i,1}(j,1:n_var) + pop{r_pf(2),1}(r2,1:n_var) -pop{r_pf(3),1}(r3,1:n_var));
               v = boundConstraint(v, pop{i,1}(j,:), [VRmin(1,:);VRmax(1,:)]);
               % Crossover
               jrand = ceil(rand*n_var);
               for h=1:n_var
                   if rand < 1 || h==jrand
                       u(j,h) = v(1,h);
                   else
                       u(j,h) = pop{i,1}(j,h);
                   end
               end
           end
           offspring = [offspring;u];
       end
    else
       for i=1: Total_pf
           if i < 3
               u = zeros(num(i),n_var);
               for j= 1 : size(pop{i,1},1)  
                   r_pf(1:3) = i;
                                                 
                   %pop{i,1}(j,1:n_var)
                   if (pop{i,1}(j,n_var + n_obj + 2)) > (pop{i,1}(j,n_var + n_obj + 4))
                       distance = sum((repmat(pop{i,1}(j,1:n_var), num(r_pf(1)) ,1) - pop{r_pf(1),1}(1:num(r_pf(1)),1:n_var)).^2,2);
                   else
                       distance = sum((repmat(pop{i,1}(j,n_var:n_var + n_obj), num(r_pf(1)) ,1) - pop{r_pf(1),1}(1:num(r_pf(1)),n_var:n_var + n_obj)).^2,2);
                   end
                   [sorted_distance_value,sorted_distance] = sort(distance);
                  
                   if  num(r_pf(1)) < 7   %1 2 3 4 5 6
                       candidate = sorted_distance;
                   else
                       r1 = 5; 
                       candidate = sorted_distance(2:r1+1);
                   end
                   r4_temp = find(pop{i,1}(candidate,n_var + n_var + 3) == max(pop{i,1}(candidate,n_var + n_var + 3)));
                   r4 = r4_temp(1);
                   r4 = candidate(r4);
                   %选择邻居中CSC
                   r2 = ceil(rand*num(r_pf(2)));
                   r3 = ceil(rand*num(r_pf(3)));
                   % Mutation
                   v = pop{i,1}(j,1:n_var) + 0.8*(pop{r_pf(1),1}(r4,1:n_var) - pop{i,1}(j,1:n_var) + pop{r_pf(2),1}(r2,1:n_var) -pop{r_pf(3),1}(r3,1:n_var));
                   v = boundConstraint(v, pop{i,1}(j,:), [VRmin(1,:);VRmax(1,:)]);
                   % Crossover
                   jrand = ceil(rand*n_var);
                   for h=1:n_var
                       if rand<1 || h==jrand
                           u(j,h) = v(1,h);
                       else
                           u(j,h) = pop{i,1}(j,h);
                       end
                    end
                end
                offspring = [offspring;u];
           else

                u = zeros(num(i),n_var);
                for j= 1 : size(pop{i,1},1) 
                    r_pf = ceil(rand(1,3)*(i-1));
                    if  num(r_pf(1)) < 6
                        r1 = num(r_pf(1)); 
                    else
                        r1 = 5;  % p = 0.1;
                    end
                  
                    conv = sum(pop{r_pf(1),1}(1:num(r_pf(1)),n_var+1:n_var+n_obj),2);
                    [sorted_conv_value,sorted_conv] = sort(conv);
                    candidate = sorted_conv(1:r1);
                    r4 = ceil(rand*length(candidate));
                    r4 = candidate(r4);
   
                    r2 = ceil(rand*num(r_pf(2)));
                    r3 = ceil(rand*num(r_pf(3)));
                    % Mutation
                    v = pop{i,1}(j,1:n_var) + 0.8*(pop{r_pf(1),1}(r4,1:n_var) - pop{i,1}(j,1:n_var) + pop{r_pf(2),1}(r2,1:n_var) -pop{r_pf(3),1}(r3,1:n_var));
                    v = boundConstraint(v, pop{i,1}(j,:), [VRmin(1,:);VRmax(1,:)]);
                    % Crossover
                    jrand = ceil(rand*n_var);
                    for h=1:n_var
                        if rand<1 || h==jrand
                            u(j,h) = v(1,h);
                        else
                            u(j,h) = pop{i,1}(j,h);
                        end
                    end
                end
                offspring = [offspring;u];
            end
       end
    end
     
    % Evaluate offspring
    for i=1:size(offspring,1);
        offspring_fitness(i,:)=feval(func_name,offspring(i,:));
        fitcount = fitcount +1 ;
    end
    temp_pop = [pos,fitness;offspring,offspring_fitness];
    temp_pop = non_domination_scd_kmeans_sort(temp_pop ,n_obj,n_var,K);
    
    % Environmental seceltion
    Parent = DSS_of_replace_decision_chromosome_kmeans(temp_pop,n_obj,n_var, NP);
    %Parent = replace_decision_chromosome_kmeans(temp_pop,n_obj,n_var, NP);
    pos = Parent(:,1:n_var);
    fitness = Parent(:,n_var+1:n_var+n_obj);
    
end
% Output
ps = Parent(Parent(:,n_var+n_obj+1)==1,1:n_var);
pf = Parent(Parent(:,n_var+n_obj+1)==1,n_var+1:n_var+n_obj);
end

