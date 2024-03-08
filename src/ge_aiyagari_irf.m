% Aiyagari model
% Finite difference implicit updating with Poisson income
% IRF
% Greg Kaplan 2024

clear;
close all;
warning('off','AutoDiff:maxmin');

%plot_path = ['/Users/judyyue/Dropbox/course_material/2024_winter/HeterogeneousHouseholds/FinalProject/graphs/discount_rate_'];
plot_path = ['/Users/judyyue/Dropbox/course_material/2024_winter/HeterogeneousHouseholds/FinalProject/graphs/preference_'];
%plot_path = ['/Users/judyyue/Dropbox/course_material/2024_winter/HeterogeneousHouseholds/FinalProject/graphs/TFP_'];
addpath(plot_path);

%% OPTIONS
options.Display                 = 1;
options.MakePlotsSteadyState    = 0;
options.MakePlotsIRF            = 1;
options.UseAutoDiffSteadyState  = 0;
options.UseAutoDiffIRF          = 1;
options.AlgorithmIRF            = 4; % 1 = fsolve:levenberg-marquardt, 2 = fsolve:trust-region, 3 = fsolve:trust-region-dogleg, 4= simple broyden


%% PARAMETERS

% preferences
param.risk_aver   = 1;
param.preference = 1;
% param.rho         = 0.01; %0.005; %quarterly
param.rhoguess    = 0.006; %0.005; %quarterly; if not calibrate this is actual

%production
param.deprec      = 0.015; %0.025;
param.capshare    = 0.4;

% income risk: discretized N(mu,sigma^2)
param.arrivalrate_y   = 0.25; % on average once per year
param.mu_y            = 1;
param.sd_y            = 0.2;
param.ny              = 5;

% asset grids
param.na          = 50;
param.amax        = 100;   %multiple of aggregate output
param.borrow_lim  = -1;    %multiple of aggregate output
param.agrid_par   = 0.4;    %1 for linear, 0 for L-shaped

%time grid
param.Ttrans    = 400; %75; % quarters for transition 400 Ttrans, 100 N trans, dt_init =0.05 works well for helictopter drop; for permanent shock use Ttrans = 50 and dt larger, like 0.2
param.nT        = 50; %750; %50; % number time steps for transtion, i.e size of Tgrid;
param.dt_initial   = 4/12; %4/52; % length of first time step (in quarters) - set to Ttrans/nTfor equally spaced steps, else dt steps will increase at a constant exponential rate;
param.Tplotmax = 100; 

% computation
param.maxiter_hjb = 100;
param.tol_hjb     = 1.0e-8;
param.delta_hjb   = 1.0e5;
param.mindV      = 1.0e-10; % for numerical stability
param.tol_irf       = 1e-7;

%tfp shock
param.tfpshock_size        = 0; %-0.05; 
param.tfpshock_decay       = 0.4; % quarterly decay rate

%discount rate shock
param.discountshock_size    = 0; %-0.001; %shock to discount rate
param.discountshock_decay   = 0.5; % quarterly decay rate


%preference shock
param.preferenceshock_size    = -0.001; %-0.001; %shock to preference
param.preferenceshock_decay   = 0.5; % quarterly decay rate



%guesses
param.rguess    = 1-exp(-param.rhoguess) - 0.0005; % a bit lower than inverse of discount rate


%% UTILITY FUNCTION

if param.risk_aver==1
    param.u = @(c)log(c);
else    
    param.u = @(c)(c.^(1-param.risk_aver)-1)./(1-param.risk_aver);
end    
param.u1 = @(c) c.^(-param.risk_aver);
param.u1inv = @(u) u.^(-1./param.risk_aver);

%% SET UP GRIDS

% assets
grids.agrid = linspace(0,1, param.na)';
grids.agrid = grids.agrid.^(1./param.agrid_par);
grids.agrid = param.borrow_lim + (param.amax - param.borrow_lim).* grids.agrid;

% asset grid spacing: for partial derivatives
grids.dagrid  =  diff(grids.agrid);
grids.dagridf = [grids.dagrid; grids.dagrid(param.na-1)];
grids.dagridb = [grids.dagrid(1); grids.dagrid];

% trapezoidal rule: for KFE and moments
grids.adelta          = zeros(param.na,1);
grids.adelta(1)       = 0.5*grids.dagrid(1);
grids.adelta(2:param.na-1)  = 0.5*grids.dagrid(1:param.na-2) + 0.5*grids.dagrid(2:param.na-1);
grids.adelta(param.na)      = 0.5*grids.dagrid(param.na-1);

% income: disretize normal distribution
if param.ny>1
    width = fzero(@(x)discrete_normal(param.ny, param.mu_y, param.sd_y,x),2);
    [temp,grids.ygrid,grids.ydist] = discrete_normal(param.ny, param.mu_y, param.sd_y,width);
else
    grids.ygrid = 1;
    grids.ydist = 1;
end

%grids on (a,y) space
grids.aagrid = repmat(grids.agrid,1, param.ny);
grids.yygrid = repmat(grids.ygrid', param.na,1);
grids.aydelta = repmat(grids.adelta,1, param.ny);
grids.daagrid = repmat(grids.dagrid,1, param.ny);
grids.daagridf = repmat(grids.dagridf,1,param.ny);
grids.daagridb = repmat(grids.dagridb,1,param.ny);
    

% income continuous time transition matrix
grids.ytrans = -param.arrivalrate_y.*eye(param.ny) + param.arrivalrate_y .* repmat(grids.ydist',param.ny,1);
grids.yytrans = kron(grids.ytrans, speye(param.na));

 % time grid
 % Ttrans is in quarters; 
if param.dt_initial == param.Ttrans/param.nT
    dt_growth = 1;
else
    f = @(x) param.Ttrans.*(x-1)/(x.^param.nT -1) - param.dt_initial;
    dt_growth = fzero(f,1.5);
end    

grids.dt_trans = param.dt_initial.*dt_growth.^(0:param.nT-1)';
grids.Tgrid = cumsum(grids.dt_trans);

%% SOLVE FOR STEADY STATE EQUILIBRIUM

%define market clearing function
fSS = @(x)market_clearing_steadystate(x,param,grids);

% rguess_AD = myAD(rguess);
% resid_AD = fSS(rguess_AD);


if options.UseAutoDiffSteadyState ==0
    options_fsolve = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-7,'SpecifyObjectiveGradient',false);
    [x,fval,exitflag,output] = fsolve(fSS,param.rguess,options_fsolve);   

elseif options.UseAutoDiffSteadyState ==1    
    options_fsolve = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-7,'SpecifyObjectiveGradient',true);
    [x,fval,exitflag,output] = fsolve(@(x)fAD(x,fSS),param.rguess,options_fsolve);    
end

%evaluate at equm prices
[~,steadystate] = market_clearing_steadystate(x,param,grids);



%% SOLVE FOR IMPULSE RESPONSE

%set shocks
shocks.tfpshock = steadystate.tfp * param.tfpshock_size .* exp(-param.tfpshock_decay.*grids.Tgrid) ;
shocks.discountshock = param.discountshock_size.*exp(-param.discountshock_decay.*(grids.Tgrid-param.dt_initial));

shocks.prefshock = param.preferenceshock_size.*exp(-param.preferenceshock_decay.*(grids.Tgrid-param.dt_initial));

%define market clearing function
fIRF = @(x)market_clearing_irf(x,param,grids,steadystate,shocks);

%guess steady state prices; (instead of getting the steady state r, first guess the steadystate k, then make your guess as the first order derivate of the productivity)
xguess = steadystate.r .* ones(param.nT,1);
ltfp = steadystate.tfp + shocks.tfpshock;
xguess = param.capshare.*ltfp.*steadystate.KLratio.^(param.capshare-1) - param.deprec;


%set optimization options

if options.AlgorithmIRF < 4 % use fsolve

    options_fsolve = optimoptions('fsolve','Display','iter','FunctionTolerance',param.tol_irf);
    if options.AlgorithmIRF==1
        options_fsolve = optimoptions(options_fsolve,'Algorithm','levenberg-marquardt');
    elseif options.AlgorithmIRF==2
        options_fsolve = optimoptions(options_fsolve,'Algorithm','trust-region-dogleg');
    elseif options.AlgorithmIRF==3
        options_fsolve = optimoptions(options_fsolve,'Algorithm','trust-region');
    end

    if options.UseAutoDiffIRF ==0
        options_fsolve = optimoptions(options_fsolve,'SpecifyObjectiveGradient',false);
        [xirf,fval,exitflag,output] = fsolve(fIRF,xguess,options_fsolve);   
    
    elseif options.UseAutoDiffIRF  ==1
        options_fsolve = optimoptions(options_fsolve,'SpecifyObjectiveGradient',true);
        [xirf,fval,exitflag,output] = fsolve(@(x)fAD(x,fIRF),xguess,options_fsolve);    
    end

elseif options.AlgorithmIRF ==4  % use simple Broyden
    %evaluate jacobian at initial guess using autodiff
    xguess_AD = myAD(xguess);
    resid_AD = fIRF(xguess_AD);
    initjac = getderivs(resid_AD);

    xirf = fn_simple_broyden(fIRF,xguess,inv(initjac),20,param.tol_irf,2);
end

%evaluate at equm prices
[res,irf] = market_clearing_irf(xirf,param,grids,steadystate,shocks);


%% MAKE PLOTS
if options.MakePlotsSteadyState ==1 
    PlotResultsSteadyState(1,param,grids,steadystate)
end

if options.MakePlotsIRF==1 
    PlotResultsIRF(2,param,grids,steadystate,irf, plot_path)
end

%% FUNCTIONS


%%
function [res,steadystate] = market_clearing_steadystate(r,param,grids)
     
    %extract some parameters and grids
    adelta = grids.adelta;
    ygrid = grids.ygrid;
    ydist = grids.ydist;
    na = param.na;
    ny = param.ny;
    agrid = grids.agrid;
    aydelta = grids.aydelta;
    pref = param.preference;
    %implied wages
    tfp = 1;
    KLratio = ((r+param.deprec)/ (tfp*param.capshare))^(1/(param.capshare-1));
    wage = tfp*(1-param.capshare).* KLratio^param.capshare;
    
    % rescale efficiency units of labor so that output = 1
    yscale = (KLratio^-param.capshare) ./ (grids.ygrid' * grids.ydist);
    
    % initialize value function
    Vguess = zeros(na,ny) + 0.*r;
    for iy = 1:ny
        Vguess(:,iy) = pref.*param.u(r.*agrid + wage.*yscale.*ygrid(iy))./ param.rhoguess;
    end


    % solve HJB
    V = Vguess;
    
    Vdiff = 1;
    iter = 0;

    while iter <= param.maxiter_hjb && Vdiff>param.tol_hjb
        iter = iter + 1;
        
        %update value function
        passedvars.r = r;
        passedvars.wage = wage;
        passedvars.yscale = yscale;
        passedvars.delta = param.delta_hjb;
        passedvars.rho = param.rhoguess;
        passedvars.preference = param.preference;
        [Vnew,A,con,sav] = UpdateHJB(V,param,grids,passedvars);

        Vdiff = max(abs(Vnew(:)-V(:)));
        
        V = Vnew;
    end 


    % solve KFE
    lgvecadj = [A' ones(na*ny,1); ones(1,na*ny) 0] \ [zeros(na*ny,1); 1];
    gvecadj = lgvecadj(1:na*ny);
    gmatadj = reshape(gvecadj,na,ny);
    gmat    = gmatadj./aydelta;

 
    %asset supply
    gamarg = sum(gmat,2);
    Ea = sum(gamarg.*agrid.*adelta);
    L = yscale.*ygrid'*ydist;
    
    %asset demand
    K = KLratio*L;
    
    res = Ea - K;

    %package output
    steadystate.r = r;
    steadystate.wage = wage;   
    steadystate.con = con;
    steadystate.sav = sav;
    steadystate.gamarg = gamarg;
    steadystate.gmat = gmat;
    steadystate.A = A;
    steadystate.V = V;
    steadystate.rho = param.rhoguess;
    steadystate.yscale = yscale;
    steadystate.tfp = tfp;
    steadystate.Ea = Ea;
    steadystate.K = K;
    steadystate.output = tfp * (K^param.capshare) *(L^(1-param.capshare));
    steadystate.Ec = sum(con(:) .* gmat(:) .* aydelta(:));
    steadystate.KLratio = KLratio;
    steadystate.preference = pref;
  
end

%%
function [res,irf] = market_clearing_irf(x,param,grids,steadystate,shocks)
     
    %extract some parameters and grids
    na = param.na;
    ny = param.ny;
    nT = param.nT;
    
    
    % shock sequence
    irf.tfp = steadystate.tfp  + shocks.tfpshock;
    
    irf.rho = steadystate.rho + shocks.discountshock; 

    irf.preference = steadystate.preference + shocks.prefshock; 



    %interest rate and wages
    irf.r = x;
    irf.KLratio = ((irf.r + param.deprec) ./ (irf.tfp.*param.capshare)).^(1/(param.capshare-1));
    irf.wage = irf.tfp.* (1-param.capshare).* irf.KLratio .^param.capshare;
        
    
    %intialize arrays for HJB
    irf.V = cell(nT,1);
    irf.A = cell(nT,1);
    irf.con = cell(nT,1);
    irf.sav = cell(nT,1);
    
    % solve HJB backward from steady state

    for it = nT:-1:1
        
        passedvars.r = irf.r(it);
        passedvars.wage = irf.wage(it);
        passedvars.yscale = steadystate.yscale;
        passedvars.delta = grids.dt_trans(it);
        passedvars.rho = irf.rho(it);
        passedvars.preference = irf.preference(it);

        if it==nT
            Vnext = steadystate.V;
        elseif it<nT
            Vnext = irf.V{it+1};
        end

        [irf.V{it},irf.A{it},irf.con{it},irf.sav{it}] = UpdateHJB(Vnext,param,grids,passedvars);
    
    end
    
    %initialize arrays for KFE
    irf.gmat = cell(nT,1);
    irf.gvecadj = cell(nT,1);

    irf.Ea = zeros(nT,1) + 0.*x(1);
    irf.L = zeros(nT,1) + 0.*x(1);
    irf.Ec = zeros(nT,1) + 0.*x(1);

    %initialize at steady state;
    irf.gmat{1} = steadystate.gmat;
    irf.gvecadj{1} = reshape(grids.aydelta .* irf.gmat{1},na*ny,1);

    % solve KFE forward
    for it = 1:nT - 1
        irf.gvecadj{it+1} = (speye(na*ny) - irf.A{it}' * grids.dt_trans(it)) \ irf.gvecadj{it};
        irf.gmat{it+1}    = reshape(irf.gvecadj{it+1},na,ny) ./ grids.aydelta;
    end
 
    %asset supply
    for it = 1:nT
        irf.gamarg{it} = sum(irf.gmat{it},2);
        irf.Ea(it) = sum(irf.gamarg{it} .* grids.agrid .* grids.adelta);
        irf.L(it) = steadystate.yscale .* grids.ygrid'*grids.ydist;
        irf.Ec(it) = sum(irf.con{it}(:) .* irf.gmat{it}(:) .* grids.aydelta(:));
    end

    %asset demand
    irf.K = irf.KLratio .* irf.L;

    % output
    irf.output = irf.tfp .* (irf.K.^param.capshare) .*(irf.L.^(1-param.capshare));
    
    % market clearing residual
    res = irf.Ea - irf.K;
%     res = irf.Ea(2:nT) - irf.K(2:nT);


end

%%
function [Vnew,A,con,sav] = UpdateHJB(V,param,grids,passedvars)

    %extract some variables
    na = param.na;
    ny = param.ny;
    agrid = grids.agrid;
    ygrid = grids.ygrid;
    daagridf = grids.daagridf;
    daagridb = grids.daagridb;
    aagrid = grids.aagrid;
    yygrid = grids.yygrid;
    yytrans = grids.yytrans;
    r = passedvars.r;
    wage = passedvars.wage;
    yscale = passedvars.yscale;
    delta = passedvars.delta;
    rho = passedvars.rho;

    pref = passedvars.preference;
    pref_inv = passedvars.preference^(-1/param.risk_aver);

    % initialize arrays as AD objects
    dVf = zeros(na,ny) + 0.*passedvars.r;
    dVb = zeros(na,ny) + 0.*passedvars.r;
    
    % forward difference         
    dVf(1:na-1,:) = (V(2:na,:)-V(1:na-1,:))./daagridf(1:na-1,:);
    
    dVf(na,:) = pref*param.u1(r.*agrid(na) + wage.*yscale.*ygrid); %state constraint
    
    % backward difference
    dVb(2:na,:) = (V(2:na,:)-V(1:na-1,:))./daagridb(2:na,:);
    dVb(1,:) = pref*param.u1(r.*agrid(1) + wage.*yscale.*ygrid); %state constraint
    
    %consumption and savings with forward difference
    conf = pref_inv*param.u1inv(max(dVf,param.mindV));
    savf = r.*aagrid + wage.*yscale.*yygrid - conf;
    Hf = pref*param.u(conf) + dVf.*savf;
    
    %consumption and savings with backward difference
    conb = pref_inv*param.u1inv(max(dVb,param.mindV));
    savb = r.*aagrid + wage.*yscale.*yygrid - conb;
    Hb = pref*param.u(conb) + dVb.*savb;
    
    %consumption and derivative with adot = 0
    con0 = r.*aagrid + wage.*yscale.*yygrid;
    dV0 = pref*param.u1(con0);
    H0 = pref*param.u(con0);
    
    % choice of forward or backward differences based on sign of drift    
    Ineither = (1-(savf>0)) .* (1-(savb<0));
    Iunique = (savb<0).*(1-(savf>0)) + (1-(savb<0)).*(savf>0);
    Iboth = (savb<0).*(savf>0);
    Ib = Iunique.*(savb<0).*(Hb>H0) + Iboth.*(Hb>Hf).*(Hb>H0);
    If = Iunique.*(savf>0).*(Hf>H0) + Iboth.*(Hf>Hb).*(Hf>H0);
    I0 = 1-Ib-If;
    
    
    %consumption, savings and utility
    con  = conf.*If + conb.*Ib + con0.*I0;
    sav  = savf.*If + savb.*Ib;    
    util = pref*param.u(con);
    utilvec = reshape(util,na*ny,1);
    
    %construct A matrix: tri-diagonal elements
    Alowdiag = -Ib.*savb./daagridb;
    Adiag = -If.*savf./daagridf + Ib.*savb./daagridb;
    Aupdiag = If.*savf./daagridf;
    
    
    %use spdiags to create A matrix for each income value
    Amat = cell(ny,1);
    for iy = 1:ny
        Amat{iy} = spdiags(Adiag(:,iy),0,na,na) + ...
                    spdiags(Alowdiag(2:na,iy),-1,na,na) + ...
                    spdiags([0;Aupdiag(1:na-1,iy)],1,na,na);
    end
    
    %combine to create large sparse A matrix
    if ny > 1
        A  = [Amat{1} sparse(na,na*(ny-1))];
        if ny>2    
            for iy = 2:ny-1
                A = [A; sparse(na,na*(iy-1)) Amat{iy} sparse(na,na*(ny-iy))];
            end
        end    
        A = [A; sparse(na,na*(ny-1)) Amat{ny}];
    elseif ny==1
        A = Amat{1};
    end
    
    % add Poisson income switches;
    A = A + yytrans;
    
    B = (rho + 1./delta)*speye(na*ny) - A;
    
    % solve linear system
    Vvec = reshape(V,na*ny,1);
    Vvecnew = B \ (utilvec + Vvec./delta);
    Vnew    = reshape(Vvecnew,na,ny);


end

%%
function PlotResultsSteadyState(fignum,param,grids,equm)
    % extract some objects
    agrid = grids.agrid;    
    ny = param.ny;
    na = param.na;
    amax = param.amax;
    borrow_lim = param.borrow_lim;
    con = equm.con;
    sav = equm.sav;
    gamarg = equm.gamarg;
    
    figure(fignum);
    
    % consumption policy function
    subplot(2,4,1);
    plot(agrid,con(:,1),'b-',agrid,con(:,ny),'r-','LineWidth',1);
    grid;
    xlim([borrow_lim amax]);
    title('Consumption Policy Function');
    legend('Lowest income state','Highest income state');

    % savings policy function
    subplot(2,4,2);
    plot(agrid,sav(:,1),'b-',agrid,sav(:,ny),'r-','LineWidth',1);
    hold on;
    plot(agrid,zeros(na,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim([borrow_lim amax]);
    title('Savings Policy Function');
    
    % consumption policy function: zoomed in
    subplot(2,4,3);
    plot(agrid,con(:,1),'b-o',agrid,con(:,ny),'r-o','LineWidth',2);
    grid;
    xlim(borrow_lim + [0 1]);
    title('Consumption: Zoomed');
    
     % savings policy function: zoomed in
    subplot(2,4,4);
    plot(agrid,sav(:,1),'b-o',agrid,sav(:,ny),'r-o','LineWidth',2);
    hold on;
    plot(agrid,zeros(na,1),'k','LineWidth',0.5);
    hold off;
    grid;
    xlim(borrow_lim + [0 1]);
    title('Savings: Zoomed');
    
    %income distribution
    subplot(2,4,5);
    bar(grids.ygrid,grids.ydist);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[0 0.5 0.5],'EdgeColor','blue','LineStyle','-');
    ylabel('')
    title('Income distribution');
    
    
    %asset distribution: CDF
    abin        = (borrow_lim:0.01:amax)'; %multiples of 1% of annual income
    nabin       = length(abin);
    
    interpacdf  = griddedInterpolant(agrid,cumsum(grids.adelta.*gamarg),'pchip');
    amass       = zeros(size(abin));
    amass(1)    = grids.adelta(1)*gamarg(1);
    amass(2:nabin) = interpacdf(abin(2:nabin)) - interpacdf(abin(1:nabin-1));
    
    % asset distribution: inverse CDF
    acdf = cumsum(grids.adelta.*gamarg);
    iacdf1 = find(acdf>0.9999999,1);
    if isempty(iacdf1)
        iacdf1 = na;
    end
    interpainvcdf  = griddedInterpolant(acdf(1:iacdf1),agrid(1:iacdf1),'linear');
    
    %plot asset distribution
    subplot(2,4,6:8);
    % bar(abin(2:nabin),amass(2:nabin),'histc');
    bar(abin(2:nabin),amass(2:nabin));
    sh  = findall(gcf,'marker','*'); delete(sh);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',[.7 .7 .7],'EdgeColor',[0.7 0.7 0.7],'LineStyle','-','LineWidth',0.01);
    
    xlim([param.borrow_lim agrid(min(iacdf1+1,na)) ]);
    
    hold on;
    bar(abin(1),amass(1),0.1,'FaceColor','black','EdgeColor','black','LineWidth',0.001);
    hold off;
    ylabel('')
    grid on;
    title('Asset distribution');
 
    
    disp(['Equilibrium interest rate: ' num2str(equm.r*100) '% (quarterly)']);
    %disp(['Discount rate: ' num2str(param.rho*100) '%  (quarterly)']);
    disp(['Mean assets: ' num2str(sum(gamarg.*grids.agrid.*grids.adelta))]);
    disp(['Fraction borrowing constrained: ' num2str(gamarg(1).*grids.adelta(1) * 100) '%']);
    disp(['10th Percentile: ' num2str(interpainvcdf(0.1))]);
    disp(['50th Percentile: ' num2str(interpainvcdf(0.5))]);
    disp(['90th Percentile: ' num2str(interpainvcdf(0.9))]);
    disp(['99th Percentile: ' num2str(interpainvcdf(0.99))]);

end

%%
function PlotResultsIRF(fignum,param,grids,steadystate,irf, plot_path)
    %plots raw IRFs, i.e. no time aggregation

    % extract some objects
    Tplotmax = param.Tplotmax; %max plot in quarters
    Tgrid = grids.Tgrid;    
    nT = param.nT;
    
    % transition
    figure(fignum);    
    
    % tfp 
    subplot(2,4,1);
    hold on;
    plot(Tgrid,irf.tfp,'b-','LineWidth',1);
    plot(Tgrid,steadystate.tfp.*ones(nT,1),'k--','LineWidth',1);
    hold off;
    grid on;
    xlim([0 Tplotmax]);
    title('TFP');
    
    
    % capital 
    subplot(2,4,2);
    hold on;
    plot(Tgrid,irf.K,'b-','LineWidth',1);
    plot(Tgrid,irf.Ea,'r--','LineWidth',1);
    plot(Tgrid,steadystate.K.*ones(nT,1),'k--','LineWidth',1);
    hold off;
    grid on;
    xlim([0 Tplotmax]);
    legend('Capital','Household Assets','Location','best')
    title('Assets');


    % output
    subplot(2,4,3);
    hold on;
    plot(Tgrid,irf.output,'b-','LineWidth',1);
    plot(Tgrid,steadystate.output.*ones(nT,1),'k--','LineWidth',1);
    hold off;
    grid on;
    xlim([0 Tplotmax]);
    title('Output');
    
    
    % consumption
    subplot(2,4,4);
    hold on;
    plot(Tgrid,irf.Ec,'b-','LineWidth',1);
    plot(Tgrid,steadystate.Ec.*ones(nT,1),'k--','LineWidth',1);
    hold off;
    grid on;
    xlim([0 Tplotmax]);
    title('Consumption');

    % interest rate
    subplot(2,4,5);
    hold on;
    plot(Tgrid,irf.r.*400,'b-','LineWidth',1);
    plot(Tgrid,steadystate.r.*ones(nT,1).*400,'k--','LineWidth',1);
    hold off;
    grid on;
    xlim([0 Tplotmax]);
    title('Interest Rate (% p.a.)');

    % wage rate
    subplot(2,4,6);
    hold on;
    plot(Tgrid,irf.wage,'b-','LineWidth',1);
    plot(Tgrid,steadystate.wage.*ones(nT,1),'k--','LineWidth',1);
    hold off;
    grid;
    xlim([0 Tplotmax]);
    title('Wage');

    set(gcf, 'Position', [100, 100, 1000, 500]);
    saveas(gcf, [plot_path, 'plot_irf.png']);


end
%%
function [f,x,p] = discrete_normal(n,mu,sigma,width)
% creates equally spaced approximation to normal distribution
% n is number of points
% mu is mean
% sigma is standard deviation
% width is the multiple of stand deviation for the width of the grid
% f is the error in the approximation
% x gives the location of the points
% p is probabilities


x = linspace(mu-width*sigma,mu+width*sigma,n)';

if n==2
    p = 0.5.*ones(n,1);
elseif n>2    
    p  = zeros(n,1);
    p(1) = normcdf(x(1) + 0.5*(x(2)-x(1)),mu,sigma);
    for i = 2:n-1
        p(i) = normcdf(x(i) + 0.5*(x(i+1)-x(i)),mu,sigma) - normcdf(x(i) - 0.5*(x(i)-x(i-1)),mu,sigma);
    end
    p(n) = 1 - sum(p(1:n-1));
end

Ex = x'*p;
SDx = sqrt((x'.^2)*p - Ex.^2);

f = SDx-sigma;
end

%% functions
function [val,jac] = fAD(x,f)
    xAD = myAD(x);
    resAD = f(xAD);
    val = getvalues(resAD);
    jac = getderivs(resAD);
end