theta_E_vals=[1,4,16,64];
theta_I_vals=[1,4,16,64];
% input_vals=5:0.5:12;
% input_vals=[7,9,10,11,12];
input_vals=13:17;
for i=1:length(theta_E_vals)
    for j=1:length(theta_I_vals)
        for k=1:length(input_vals)
            theta_E=theta_E_vals(i);
            theta_I=theta_I_vals(j);
            input=input_vals(k);
            input_E=input;
            input_I=input;
            %%%% FG %%%%%%%%
            folder_name="Model_sim_FG_theta_"+string(theta_E)+"_"+string(theta_I)+"_input_"+string(input);
            % folder_name="Model_sim_FG_input_"+string(input);
            indicator=1;
            generate_WC_OU_model_gamma(indicator,folder_name,theta_E,theta_I,input_E,input_I);
            %%%% SG %%%%%%%%
            folder_name="Model_sim_SG_theta_"+string(theta_E)+"_"+string(theta_I)+"_input_"+string(input);
            % folder_name="Model_sim_SG_input_"+string(input);
            indicator=0;
            generate_WC_OU_model_gamma(indicator,folder_name,theta_E,theta_I,input_E,input_I);
        end
    end 
end 