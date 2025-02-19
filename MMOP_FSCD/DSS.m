                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           function f = DSS(current_pop,temp_NP,M,V)
   
    cscd_back = current_pop(:,V + M + 3 );
    cscd_min = min(cscd_back);
    cscd_max = max(cscd_back);
    Archives_pop = current_pop(find(cscd_back == cscd_max),:);
    current_pop(find(cscd_back == cscd_max),:) = []; 
    while size(Archives_pop,1) < temp_NP
        sum_distance_dec = 0;
        sum_distance_obj = 0;
        for k1 = 1 : size(current_pop,1)
            distance_dec = min(sqrt(sum((current_pop(k1,1 : V) - Archives_pop(:,1 : V)).^2,2)));
            sum_distance_dec = sum_distance_dec + distance_dec;     
            current_pop(k1,V + M + 5) = distance_dec;
            distance_obj = min(sqrt(sum((current_pop(k1,V + 1:V + M) - Archives_pop(:,V + 1:V + M)).^2,2)));
            sum_distance_obj = sum_distance_obj + distance_obj;
            current_pop(k1,V + M + 7) = distance_obj;

        end
        avg_distance_dec = sum_distance_dec / size(current_pop,1);
        avg_distance_obj = sum_distance_obj / size(current_pop,1);
        for k1 = 1 : size(current_pop,1)
           if  current_pop(k1,V + M + 5) > avg_distance_dec || current_pop(k1,V + M + 7) > avg_distance_obj
               current_pop(k1,V + M + 6)=max(current_pop(k1,V + M + 5),current_pop(k1,V + M + 7)); % Eq. (6) in the paper
           else
               current_pop(k1,V + M + 6)=min(current_pop(k1,V + M + 5),current_pop(k1,V + M + 7)); % Eq. (7) in the paper
           end
        end
        
        current_pop(:,V + M + 6) = current_pop(:,V + M + 6) + current_pop(:,V + M + 3);
        [temp,index] = sort(current_pop(:,M + V + 6),"descend");
        Archives_pop = [Archives_pop;current_pop(index(1),1:V + M + 4)];    
        current_pop(index(1),:) = [];
    end
    f = Archives_pop;  
end

