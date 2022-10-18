function []=strain_mvc_plots_v2Siemens(struct,roi_name,trigger)

dt=trigger(1);
mkdir(roi_name)
cd(roi_name)


% displacement
limits=[-10, 10];
v=struct.dx;
sd=struct.dx_sd;
strain_plot_v2Siemens(10*v,10*sd,'$\Delta_{x} \quad [\mathrm{mm}]$','dx',limits,dt);

v=struct.dy;
sd=struct.dy_sd;
strain_plot_v2Siemens(10*v,10*sd,'$\Delta_{y} \quad [\mathrm{mm}]$','dy',limits,dt);

v=struct.dz;
sd=struct.dz_sd;
strain_plot_v2Siemens(10*v,10*sd,'$\Delta_{z} \quad [\mathrm{mm}]$','dz',limits,dt);


%velocity
limits=[-3, 3];
v=struct.vx;
sd=struct.vx_sd;
strain_plot_v2Siemens(v,sd,'$v_{x} \quad [\mathrm{cm/s}]$','vx',limits,dt);

v=struct.vy;
sd=struct.vy_sd;
strain_plot_v2Siemens(v,sd,'$v_{y} \quad [\mathrm{cm/s}]$','vy',limits,dt);

v=struct.vz;
sd=struct.vz_sd;
strain_plot_v2Siemens(v,sd,'$v_{z} \quad [\mathrm{cm/s}]$','vz',limits,dt);


% euler strain
limits=[-0.6, 0.6];
v=struct.E_lambda(:,1);
sd=struct.E_lambda_sd(:,1);
strain_plot_v2Siemens(v,sd,'$E_{\lambda_{1}} \quad [\mathrm{mm/mm}]$','E1',limits,dt);

v=struct.E_lambda(:,2);
sd=struct.E_lambda_sd(:,2);
strain_plot_v2Siemens(v,sd,'$E_{\lambda_{2}} \quad [\mathrm{mm/mm}]$','E2',limits,dt);

v=struct.E_lambda(:,3);
sd=struct.E_lambda_sd(:,3);
strain_plot_v2Siemens(v,sd,'$E_{\lambda_{3}} \quad [\mathrm{mm/mm}]$','E3',limits,dt);

v=struct.ShearE_max;
sd=struct.ShearE_max_sd(:,1);
strain_plot_v2Siemens(v,sd,'$E_{max} \quad [\mathrm{mm/mm}]$','Emax',limits,dt);

v=struct.E_Volumetric;
sd=struct.E_Volumetric_sd;
strain_plot_v2Siemens(v,sd,'$E_{vol} \quad [\mathrm{mm^3/mm^3}]$','EV',limits,dt);

% lagrangian strain
limits=[-0.6, 0.6];
v=struct.L_lambda(:,1);
sd=struct.L_lambda_sd(:,1);
strain_plot_v2Siemens(v,sd,'$L_{\lambda_{1}} \quad [\mathrm{mm/mm}]$','L1',limits,dt);

v=struct.L_lambda(:,2);
sd=struct.L_lambda_sd(:,2);
strain_plot_v2Siemens(v,sd,'$L_{\lambda_{2}} \quad [\mathrm{mm/mm}]$','L2',limits,dt);

v=struct.L_lambda(:,3);
sd=struct.L_lambda_sd(:,3);
strain_plot_v2Siemens(v,sd,'$L_{\lambda_{3}} \quad [\mathrm{mm/mm}]$','L3',limits,dt);

v=struct.ShearL_max;
sd=struct.ShearL_max_sd(:,1);
strain_plot_v2Siemens(v,sd,'$L_{max} \quad [\mathrm{mm/mm}]$','Lmax',limits,dt);

v=struct.L_Volumetric;
sd=struct.L_Volumetric_sd;
strain_plot_v2Siemens(v,sd,'$L_{vol} \quad [\mathrm{mm^3/mm^3}]$','LV',limits,dt);

% euler strain rate
limits=[-1000,1000];
v=struct.SR_E_lambda(:,1);
sd=struct.SR_E_lambda_sd(:,1);
strain_plot_v2Siemens(v,sd,'$SR_{\lambda_{1}} \quad [\mathrm{m s^{-1}}]$','SRE1',limits,dt);

v=struct.SR_E_lambda(:,2);
sd=struct.SR_E_lambda_sd(:,2);
strain_plot_v2Siemens(v,sd,'$SR_{\lambda_{2}} \quad [\mathrm{m s^{-1}}]$','SRE2',limits,dt);

v=struct.SR_E_lambda(:,3);
sd=struct.SR_E_lambda_sd(:,3);
strain_plot_v2Siemens(v,sd,'$SR_{\lambda_{3}} \quad [\mathrm{m s^{-1}}]$','SRE3',limits,dt);

limits=[-1500, 1500];
v=struct.ShearSR_E_max;
sd=struct.ShearSR_E_max_sd(:,1);
strain_plot_v2Siemens(v,sd,'$SR_{max} \quad [\mathrm{m s^{-1}}]$','SREmax',limits,dt);


% lagrangian strain rate Lagrangian
limits=[-1000,1000];

v=struct.SR_lambda(:,1);
sd=struct.SR_lambda_sd(:,1);
strain_plot_v2Siemens(v,sd,'$SR_{\lambda_{1}} \quad [\mathrm{m s^{-1}}]$','SR1',limits,dt);

v=struct.SR_lambda(:,2);
sd=struct.SR_lambda_sd(:,2);
strain_plot_v2Siemens(v,sd,'$SR_{\lambda_{2}} \quad [\mathrm{m s^{-1}}]$','SR2',limits,dt);

v=struct.SR_lambda(:,3);
sd=struct.SR_lambda_sd(:,3);
strain_plot_v2Siemens(v,sd,'$SR_{\lambda_{3}} \quad [\mathrm{m s^{-1}}]$','SR3',limits,dt);


limits=[-1500, 1500];
v=struct.ShearSR_max;
sd=struct.ShearSR_max_sd(:,1);
strain_plot_v2Siemens(v,sd,'$SR_{max} \quad [\mathrm{m s^{-1}}]$','SRmax',limits,dt);




cd ..

end

