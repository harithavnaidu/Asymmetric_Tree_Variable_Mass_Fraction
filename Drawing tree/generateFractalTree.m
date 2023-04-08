function [Ngen,Radius,Lenght,ID] = generateFractalTree
% Asymetric Fractal Tree Generation - by Vasilina in August 20,2018
% ------------------- Fractal tree parameters description------------------
% Given:
% eta - area ratio, eta = (r_d1^2 + r_d2^2)/r_parent^2 
% xi -  Power exponent for Murray's law
% Ngen -   Total number of all generations
% Rroot, Lroot - Root vessel size
% lrr - length-to-radius ratio, lrr = 50 in [Olufsen]
%
% Computed:
% gamma - asymmetry ratio, gamma = (r_d2/r_d1)^2, 0 < gamma <= 1
% alpha, beta - Scaling factors for asymmetric branching
% 
%--------------------------------------------------------------------------
% clc;
% clear all;
% close all;
% format long;
% workspace

    %----------------- Parameters------------------------
    Rroot = 0.0055;  % root vessel
    Ngen = 15; %6; % number of generations
    symflag = 1; %symmetry flag

    % values from [Olufsen, 2000]
    eta = 1.16; % area ratio >1
    xi = 2.76;  % Murray's law exponent

    %---------------- Compute asymmetry ratio-------------
    % get gamma, solve algebraic nonlinear Eq.(21) in [Olufsen, 2000]
    syms y 
    f = eta*(1.0 + y^(xi/2.0))^(2.0/xi) - 1.0 - y;
    sol = vpasolve(f,y,[0,2]);
    gamma = double(sol);

    %for symmetric tree
    if symflag==1
    gamma=1.0;
    end

    %---------------- Get bifurcation parameters (radius rule) ---------
    alpha = (1.0 + gamma^(xi/2.0))^(-1.0/xi);
    beta = alpha * sqrt(gamma);

    %- -------------- ID relates unique branch to radius indexing------
    maxs = 2^(Ngen-1);  %need dynamic allocation for large Ngen
    ID = zeros(Ngen,maxs);
    ID(1,1)=1;
    for k=2:Ngen
        % # branches at previous gen, current gen has 2*nprevgen=2^(k-1)
        nprevgen = 2^(k-2);   
        for s = 1: nprevgen  %half of all branches in k-th gen.
            i = ID(k-1,s);
            ID(k,s) = i;                %-> R(i,k+1-i)
            ID(k,nprevgen+s) = i + 1;   %-> R(i+1,k-i)
        end
    end

    %------ Get Radius array----------
    Radius=zeros(Ngen,Ngen);
    Lenght = zeros(Ngen,Ngen);
    for k = 1:Ngen
        for i=1:k
            j=k+1-i;
            Radius(i,j) = alpha^(i-1)*beta^(j-1)*Rroot;
            Lenght(i,j) = LengthSegmentk(Radius(i,j));
        end
    end
end