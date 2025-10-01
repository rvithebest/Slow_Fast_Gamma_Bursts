theta_E_vals=[1,4,16,64];
theta_I_vals=[1,4,16,64];
for i=1:length(theta_E_vals)
    for j=1:length(theta_I_vals)
        theta_E=theta_E_vals(i);
        theta_I=theta_I_vals(j);
        %%%% FG %%%%%%%%
        folder_name="Model_sim_FG_theta_"+string(theta_E)+"_"+string(theta_I);
        indicator=1;
        generate_WC_OU_model_gamma(indicator,folder_name,theta_E,theta_I);
        %%%% SG %%%%%%%%
         folder_name="Model_sim_SG_theta_"+string(theta_E)+"_"+string(theta_I);
        indicator=0;
        generate_WC_OU_model_gamma(indicator,folder_name,theta_E,theta_I);
    end 
end 