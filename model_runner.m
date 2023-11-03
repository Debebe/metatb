
%Author: Debebe Shaweno
%The University of Melbourne
%This code calls the Ode function

close all
clear variables
clc
rate_birth_withgrowth = 0.033 * [1, 1, 1];                                         
rate_fastprogression = 0.4*ones(1, 3);                                   
rate_reactivation = 0.002*ones(1, 3);                                 
rate_stabilization = 3.6*ones(1, 3);                                           
rate_naturalmortality = 0.0154*ones(1, 3);                    
rate_untreatedmortality = 0.125 * ones(1, 3);                      
rate_naturalrecovery = 0.205 *ones(1, 3);                          
proportion_infectious =0.4*ones(1, 3);                        
proportion_detected= [0.65, 0.6, 0.6];
rate_casedetection = (proportion_detected)/(1-proportion_detected) * ...
                                (rate_untreatedmortality+ rate_naturalrecovery) ;     
rate_relapse = 0.003*(0.1*rate_naturalrecovery/(rate_naturalrecovery + ...
    rate_casedetection) + 0.9*rate_casedetection/(rate_casedetection + rate_naturalrecovery))* ones(1, 3);
proportion_partialimmune = [0.79, 0.79,  0.79];                     
% ---------------------------------------
pop= 230000;

%    Initial conditions for compartments
time_start = 0;
time_final = 1000;
compartment_susceptible_initial = [0.22*pop, 0.46*pop,  0.32*pop];
compartment_earlylatent_initial = zeros(1, 3);
compartment_latelatent_initial = zeros(1, 3);
compartment_active_initial = [1, 1, 1 ];
compartment_recovered_initial = zeros(1,3);
                              
% time span
time_span = [time_start, time_final];
%initial state values 
compartment_initial_matrix = ... 
    [compartment_susceptible_initial; ...
    compartment_earlylatent_initial; ...
    compartment_latelatent_initial; ...
    compartment_active_initial; ...
    compartment_recovered_initial; ...
    zeros(6,3)];

compartment_initial_matrix=reshape(compartment_initial_matrix,11,3);
% Parameters
  parameter = [proportion_infectious;
                proportion_partialimmune; 
                rate_naturalmortality; ...
                rate_fastprogression; ...
                rate_stabilization; ...
                rate_naturalrecovery; ...
                rate_reactivation;  ...
                rate_untreatedmortality; 
                rate_casedetection; ...
                rate_birth_withgrowth; ... 
                 rate_relapse];
        options = odeset('RelTol',1e-5);
        
        %
        beta1=54;
        beta2=15;
        beta3=15;
        p=0.1;
        
        mixingmatrix=[beta1, beta2*p, beta3*p*p; beta1*p, beta2, beta3*p; beta1*p*p, beta2*p, beta3];
        
        [time_output, population_output] = spatial_mathematical_model(time_span, ...
         compartment_initial_matrix, parameter, mixingmatrix);
         population_output_reshaped ...
            = reshape(population_output, size(population_output, 1), 11, 3); 

     % Prevalent cases
    total_prevalentcases_dstb=sum([population_output(:,4), ...
    population_output(:,16),population_output(:,28)], 2); 
    % Incident cases
    change_newcases = diff(population_output_reshaped);
    change_time =   diff(time_output);
    % Output
    population_hotspots = sum(population_output_reshaped(:,1:5,1),2); 
    population_adjacent = sum(population_output_reshaped(:,1:5,2),2); 
    population_remote = sum(population_output_reshaped(:,1:5,3),2); 
    population_total = population_hotspots + population_adjacent + population_remote;
                    
    rate_notification_hotspots = 1e5 *sum(change_newcases(:, 7, 1)) ./ sum(change_time)./sum(population_hotspots(end-1, :));
    rate_notification_adjacent = 1e5 *sum(change_newcases(:, 7, 2))./ sum(change_time )./sum( population_adjacent(end-1, :));
    rate_notification_remote =    1e5 *sum(change_newcases(:, 7, 3))./ sum(change_time) ./ sum(population_remote(end-1, :));
    
    mean_hot_notification = mean(rate_notification_hotspots );
    mean_adj_notification = mean(rate_notification_adjacent );
    mean_rem_notification = mean(rate_notification_remote );
    mean_notification = [ mean_hot_notification; mean_adj_notification; mean_rem_notification];
  
 %   Develp an algorithm that finds the beta value that leads the model  best fitting model to notifications data 
    notification_data=[137, 181,168, 220, 184; 42, 66,	93,	108, 110;  61,	91,	57,	75, 87]';
    Pop_data = compartment_susceptible_initial;
    per_pop = diag(1./Pop_data);
    notification_per_hundredk=100000*notification_data*per_pop;
  
  loglikelihood_current = -Inf; 
 
  % iterate to find best mixing matrix
    parameter_current = [beta1; beta2;  beta3; p];
  
  parameter_candidate = parameter_current; itmat = []; burnin=110000; endit= 120000;
   for it=1:endit
       if rem(it, 10)==0, disp(it), end
       acceptbeta1=0; acceptbeta2=0; acceptbeta3= 0; acceptp=0;
        for pchange = 1: 4
            parameter_candidate = parameter_current;
            if pchange==1
            parameter_candidate(1) = parameter_candidate(1) + 0.1 *randn;
            elseif pchange==2
             parameter_candidate(2) = parameter_candidate(2) + 0.1*randn; 
            elseif pchange==3
             parameter_candidate(3) = parameter_candidate(3) + 0.1*randn;       
            else
              parameter_candidate(4) = parameter_candidate(4) + 0.001* randn; 
             end
              if min(parameter_candidate) >= 0, 
        beta1 = parameter_candidate(1);
        beta2 = parameter_candidate(2);
        beta3 = parameter_candidate(3);
        p       =   parameter_candidate(4);
         mixingmatrix = [beta1, beta2*p, beta3*p*p; beta1*p, beta2, beta3*p; beta1*p*p, beta2*p, beta3];
        
 % set up for the 5 years
         compartment_initial_matrix=population_output_reshaped(end, :, :);
         compartment_initial_matrix=reshape(compartment_initial_matrix,11, 3);
         compartment_initial_matrix(6:11, :)=zeros(6,3);
         
         time_span = [0, 500];
            
        [time_output, population_output] = spatial_mathematical_model(time_span, ...
         compartment_initial_matrix, parameter, mixingmatrix);
     
        time_start = 0;
        time_final = 1;
        time_span = [time_start, time_final];
        for year=1:5
         compartment_initial_matrix=population_output_reshaped(end, :, :);
         compartment_initial_matrix=reshape(compartment_initial_matrix,11, 3);
         compartment_initial_matrix(6:11, :)=zeros(6,3);
        
        [time_output, population_output] = spatial_mathematical_model(time_span, ...
         compartment_initial_matrix, parameter, mixingmatrix);
     
        population_output_reshaped ...
            = reshape(population_output, size(population_output, 1), 11, 3); 
        % ................................................................
        total_prevalentcases_dstb=sum([population_output(:,4), ...
        population_output(:,16),population_output(:,28)], 2); 
 
        change_newcases = diff(population_output_reshaped);
        change_time =   diff(time_output);
        population_hotspots = sum(population_output_reshaped(:,1:5,1),2); 
        population_adjacent = sum(population_output_reshaped(:,1:5,2),2); 
        population_remote = sum(population_output_reshaped(:,1:5,3),2); 
        population_total = population_hotspots + population_adjacent + ...
                          population_remote;

        rate_notification_hotspots = 1e5 *sum(change_newcases(:, 7, 1)) ./ sum(change_time)./sum(population_hotspots(end-1, :));
        rate_notification_adjacent = 1e5 *sum(change_newcases(:, 7, 2))./ sum(change_time )./sum( population_adjacent(end-1, :));
        rate_notification_remote =    1e5 *sum(change_newcases(:, 7, 3))./ sum(change_time) ./ sum(population_remote(end-1, :));

        mean_hot_notification = mean(rate_notification_hotspots );
        mean_adj_notification = mean(rate_notification_adjacent );
        mean_rem_notification = mean(rate_notification_remote );
        mean_notification = [ mean_hot_notification; mean_adj_notification; mean_rem_notification];
        mean_notification_year(:, year) = mean_notification;
        end
%
       loglikelihood_candidate= 0;
       nrows=size(notification_per_hundredk, 1);
       ncols=size(notification_per_hundredk, 2);
       mean_notification_year_place=mean_notification_year';
  for a=1:nrows                       %years
      for b=1:ncols                     %regions
          data_output_year_place= notification_per_hundredk(a, b);
          model_number_year_place=mean_notification_year_place(a, b);
          compare_data_with_model=poisspdf(round(data_output_year_place), model_number_year_place);
          loglikelihood_candidate = loglikelihood_candidate + log(compare_data_with_model);
      end
  end
       
  if exp(loglikelihood_candidate - loglikelihood_current) > rand, loglikelihood_current= loglikelihood_candidate; parameter_current = parameter_candidate;  
   if pchange==1, acceptbeta1=acceptbeta1+1; elseif pchange==2,acceptbeta2=acceptbeta2+1;elseif pchange==3,acceptbeta3=acceptbeta3+1; else  acceptp=acceptp+1; end 
  end
               end
        end
              
  itmat = [itmat; acceptbeta1, acceptbeta2, acceptbeta3, acceptp,  loglikelihood_current,  parameter_current', mean_notification',  loglikelihood_candidate];
    
     end
   
   
   
