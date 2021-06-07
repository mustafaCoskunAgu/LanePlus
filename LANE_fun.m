function H = LANE_fun(methods,Net,Attri,LabelY,d,alpha1,alpha2,varargin)
%Jointly embed labels and attriuted network into embedding representation H
%     H = LANE_fun(Net,Attri,LabelY,d,alpha1,alpha2,numiter);
%     H = AANE_fun(Net,Attri,d,alpha1,alpha2,numiter);
% 
%          Net   is the weighted adjacency matrix
%         Attri  is the attribute information matrix with row denotes nodes
%        LabelY  is the label information matrix
%          d     is the dimension of the embedding representation
%         alpha1 is the weight for node attribute information
%         alpha2 is the weight for label information
%        numiter is the max number of iteration

%   Copyright 2017, Xiao Huang and Jundong Li.
%   $Revision: 1.0.0 $  $Date: 2017/10/18 00:00:00 $



n = size(Net,1);

option = methods;
switch(option)
    case 'Lap'
         %fprintf('Laplacian Based \n');
         LG = norLap(Net); % Normalized Network Laplacian
         LA = norLap(Attri); % Normalized Node Attributes Laplacian
    case 'CN'
        %fprintf('Common Neighbor Based\n');
        LG = CN(Net);
        LA = CN(Attri);
    case 'AA'
        %fprintf('Adamic-Adar Based\n');
        LG = AA(Net);
        LA = AA(Attri);
    case 'RA'
        %fprintf('Resource Allocation BAsed\n');
        LG = RA(Net);
        LA = RA(Attri);
    case 'HDI'
        %fprintf('HDI Based\n');
        LG = HDI(Net);
        LA = HDI(Attri);
    case 'HPI'
        %fprintf('HPI Based\n');
        LG = HPI(Net);
        LA = HPI(Attri);
    case 'Salton'
        %fprintf('Salton Based\n');
        LG = Salton(Net);
        LA = Salton(Attri);
    case 'RWR'
        %fprintf('Smoothed RWR Based\n');
        LG = RWR(Net);
        LA = RWR(Attri);
    otherwise
         fprintf('Give a proper method\n');
      
end
        
% LG = Salton(Net);
% LA = Salton(Attri);
UAUAT = zeros(n,n); % UA*UA^T
opts.disp = 0;

if isempty(varargin)
    %% Unsupervised attriuted network embedding
    % Input of Parameters
    numiter = alpha2; % the max number of iteration
    beta1 = d; % the weight for node attribute information
    beta2 = alpha1; % the weight for the correlations
    d = LabelY; % the dimension of the embedding representation
    H = zeros(n,d);
    for i = 1:numiter
        HHT = H*H';
        TotalLG1 = LG+beta2*UAUAT+HHT;
        [UG,~] = eigs(.5*(TotalLG1+TotalLG1'),d,'LA',opts);
        UGUGT = UG*UG';
        
        TotalLA = beta1*LA+beta2*UGUGT+HHT;
        [UA,~] = eigs(.5*(TotalLA+TotalLA'),d,'LA',opts);
        UAUAT = UA*UA';
        
        TotalLH = UAUAT+UGUGT;
        [H,~] = eigs(.5*(TotalLH+TotalLH'),d,'LA',opts);
    end
else
    %% Supervised attriuted network embedding
    numiter = varargin{1}; % the max number of iteration
    H = zeros(n,d);
    switch(option)
        case 'Lap'
             %fprintf('Laplacian Based \n');
             %LG = norLap(Net); % Normalized Network Laplacian
             %LA = norLap(Attri); % Normalized Node Attributes Laplacian
             LY = norLap(LabelY*LabelY');
        case 'CN'
            %fprintf('Common Neighbor Based\n');
            %LG = CN(Net);
            %LA = CN(Attri);
            %LY = CN(LabelY*LabelY');
            LY = norLap(LabelY*LabelY');
        case 'AA'
            %fprintf('Adamic-Adar Based\n');
            %LG = AA(Net);
            %LA = AA(Attri);
            %LY = AA(LabelY*LabelY');
            LY = norLap(LabelY*LabelY');
        case 'RA'
            %fprintf('Resource Allocation Based\n');
            %LG = RA(Net);
            %LA = RA(Attri);
            %LY = RA(LabelY*LabelY');
            LY = norLap(LabelY*LabelY');
        case 'HDI'
            %fprintf('HDI Based\n');
            %LG = HDI(Net);
            %LA = HDI(Attri);
            %LY = HDI(LabelY*LabelY');
            LY = norLap(LabelY*LabelY');
        case 'HPI'
            %fprintf('HPI Based\n');
            %LG = HPI(Net);
            %LA = HPI(Attri);
            %LY = HPI(LabelY*LabelY');
            LY = norLap(LabelY*LabelY');
        case 'Salton'
            %fprintf('Salton Based\n');
            %LG = Salton(Net);
            %LA = Salton(Attri);
            %LY = Salton(LabelY*LabelY');
            LY = norLap(LabelY*LabelY');
        case 'RWR'
            %fprintf('Smoothed RWR Based\n');
            %LG = RWR(Net);
            %LA = RWR(Attri);
            %LY = RWR(LabelY*LabelY');
            LY = norLap(LabelY*LabelY');
        otherwise
             fprintf('Give a proper method\n');

    end
     % Normalized Label Laplacian
    %LY = Salton(LabelY*LabelY');
    UYUYT = zeros(n,n); % UY*UY^T
    % Iterations
    for i = 1:numiter
        HHT = H*H';
        TotalLG1 = LG+alpha1*UAUAT+alpha2*UYUYT+HHT;
        [UG,~] = eigs(.5*(TotalLG1+TotalLG1'),d,'LA',opts);
        UGUGT = UG*UG';
        
        TotalLA = alpha1*(LA+UGUGT)+HHT;
        [UA,~] = eigs(.5*(TotalLA+TotalLA'),d,'LA',opts);
        UAUAT = UA*UA';
        
        TotalLY = alpha2*(LY+UGUGT)+HHT;
        [UY,~] = eigs(.5*(TotalLY+TotalLY'),d,'LA',opts);
        UYUYT = UY*UY';
        
        TotalLH = UAUAT+UGUGT+UYUYT;
        [H,~] = eigs(.5*(TotalLH+TotalLH'),d,'LA',opts);
    end
end

end


    function LapX = norLap(InpX)
    % Compute the normalized graph Laplacian of InpX
        InpX = InpX'; % Transpose for speedup
        InpX = bsxfun(@rdivide,InpX,sum(InpX.^2).^.5); % Normalize
        InpX(isnan(InpX)) = 0;
        SX = InpX'*InpX;
        nX = length(SX);
        SX(1:nX+1:nX^2) = 1+10^-6;
        DXInv = spdiags(full(sum(SX,2)).^(-.5),0,nX,nX);
        LapX = DXInv*SX*DXInv;
        LapX = .5*(LapX+LapX');
    end
    
    function [ sim ] = Salton( InpX )
            % Compute the normalized graph Laplacian of InpX
        InpX = InpX'; % Transpose for speedup
        InpX = bsxfun(@rdivide,InpX,sum(InpX.^2).^.5); % Normalize
        InpX(isnan(InpX)) = 0;
        SX = InpX'*InpX;
        nX = length(SX);
        SX(1:nX+1:nX^2) = 1+10^-6;
        tempdeg = repmat((sum(SX,2)).^0.5,[1,size(SX,1)]);       
        tempdeg = tempdeg .* tempdeg';            
        sim = SX * SX;              

        sim = sim./tempdeg;                 
        sim(isnan(sim)) = 0; sim(isinf(sim)) = 0;
    end

    
    function [ sim ] = HDI( InpX )
        InpX = InpX'; % Transpose for speedup
        InpX = bsxfun(@rdivide,InpX,sum(InpX.^2).^.5); % Normalize
        InpX(isnan(InpX)) = 0;
        SX = InpX'*InpX;
        nX = length(SX);
        SX(1:nX+1:nX^2) = 1+10^-6;
        sim = SX * SX;      
        deg_row = repmat(sum(SX,1), [size(SX,1),1]);
        deg_row = deg_row .* spones(sim);
        deg_row = max(deg_row, deg_row');  
        sim = sim ./ deg_row; clear deg_row;     
        sim(isnan(sim)) = 0; sim(isinf(sim)) = 0;   
    end
    function [ sim ] = CN( InpX )
    %sim = train * train;  
        InpX = InpX'; % Transpose for speedup
        InpX = bsxfun(@rdivide,InpX,sum(InpX.^2).^.5); % Normalize
        InpX(isnan(InpX)) = 0;
        SX = InpX'*InpX;
        nX = length(SX);
        SX(1:nX+1:nX^2) = 1+10^-6;
        sim = SX * SX; 
    end
    
    function [ sim ] = AA( InpX )
    
        InpX = InpX'; % Transpose for speedup
        InpX = bsxfun(@rdivide,InpX,sum(InpX.^2).^.5); % Normalize
        InpX(isnan(InpX)) = 0;
        SX = InpX'*InpX;
        nX = length(SX);
        SX(1:nX+1:nX^2) = 1+10^-6;

        train1 = SX ./ repmat(log(sum(SX,2)),[1,size(SX,1)]); 
        train1(isnan(train1)) = 0; 

        sim = SX * train1;   clear train1;  
    end
    
    function [ sim ] = RA( InpX )   
        InpX = InpX'; % Transpose for speedup
        InpX = bsxfun(@rdivide,InpX,sum(InpX.^2).^.5); % Normalize
        InpX(isnan(InpX)) = 0;
        SX = InpX'*InpX;
        nX = length(SX);
        SX(1:nX+1:nX^2) = 1+10^-6;
        train1 = SX ./ repmat(sum(SX,2),[1,size(SX,1)]); 
        train1(isnan(train1)) = 0; 
        train1(isinf(train1)) = 0;
        sim = SX * train1;  clear train1; 
    end
    
    function [ sim ] = HPI( InpX )
        InpX = InpX'; % Transpose for speedup
        InpX = bsxfun(@rdivide,InpX,sum(InpX.^2).^.5); % Normalize
        InpX(isnan(InpX)) = 0;
        SX = InpX'*InpX;
        nX = length(SX);
        SX(1:nX+1:nX^2) = 1+10^-6;
        sim = SX * SX;      
        deg_row = repmat(sum(SX,1), [size(SX,1),1]);
        deg_row = deg_row .* spones(sim);
        deg_row = min(deg_row, deg_row');  
        sim = sim ./ deg_row; clear deg_row;     
        sim(isnan(sim)) = 0; sim(isinf(sim)) = 0;   
    end
    
    
    function [sim] = RWR( InpX)
        InpX = InpX'; % Transpose for speedup
        InpX = bsxfun(@rdivide,InpX,sum(InpX.^2).^.5); % Normalize
        InpX(isnan(InpX)) = 0;
        SX = InpX'*InpX;
        nX = length(SX);
        SX(1:nX+1:nX^2) = 1+10^-6;
        DXInv = spdiags(full(sum(SX,2)).^(-1),0,nX,nX);
        LapX = DXInv*SX;
        %LapX = .5*(LapX+LapX');
        
        n = size(LapX,1);
        restart_prob = 0.5;
        sim = (eye(n) - (1 - restart_prob) * LapX) \ (restart_prob * eye(n));
        %R = log(Q + 1/n);
        %sim = R*R';
    end