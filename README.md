# HeterogenousAgentModel

The final project for Heterogenous Agent Model course. Please let me know if I should keep the channel be private given the majority of the code is the sample code we have in class. 

To get the IRF of an MIT shock
  in TFP: run ge_ayagari_irf.m, set param.tfpshock_size = -0.05, keep other shocksize be 0, and change the plot_path to local directory
  in Discount Rate: run ge_ayagari_irf.m, set param.discountshock_size = -0.001, keep other shocksize be 0, and change the plot_path to local directory
  in Preference: run ge_ayagari_irf.m, set param.preferenceshock_size = 0.001, keep other shocksize be 0, and change the plot_path to local directory
  in Borrowing Constraint (approximate with asset quality): run ge_ayagari_irf_borrowlimit2.m, and change the plot_path to local directory
