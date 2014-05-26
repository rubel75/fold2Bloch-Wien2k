function ubs
% Plot undolded band structure

%% Init. parameters
KPATH = [0 1/2 0; ...
         0 0 0; ...
         0 0 1/2]; % k-point path
nk = 100*2; % number of divisions along k-path
sk = 0.007; % smearing in k-space
ERANGE = [-1 0]; % energy range for plot (Ry)
ne = 80*2; % number of divisions in energy space
se = 0.002*3; % smearing in the energy
finpt = '6-atom2D.unfolded'; % input file name
Ef = 0.0; % Fermi energy (Ry)
ry2ev = 13.605698066; % Ry -> eV conversion factor
pwr = 1; % power for result plotting
         % 1 - linear scale, 1/2 - sqrt, etc.

%% INITIALIZATION
KPT = kgen(KPATH,nk); % mesh the k-poin path
ENE = linspace( ERANGE(1) , ERANGE(2) , ne );
[KEIG, EIG, W] = readinput(finpt); % read input data from file
% EIG - energy eigenvalues
% KEIG - k-list for eigenvalues
% W - list of characters

%% MAIN
SD = zeros(ne,nk); % spectral density
for ik = 1 : nk
    progrbar( ik/nk , 'Progress:' );
    for ie = 1 : ne
        for j = 1 : length(EIG)
            if EIG(j) > ERANGE(1) && EIG(j) < ERANGE(2)
                DKV = KEIG(j,:) - KPT(ik,:);
                dk = sqrt(sum( DKV.^2 ));
                de = abs( EIG(j)-ENE(ie) );
                wk = gauss(dk,0,sk);
                we = gauss(de,0,se);
                SD(ie,ik) = SD(ie,ik) + wk*we*W(j);
            end
        end
    end
end

%% Plot results
imagesc((1:nk),(ENE-Ef)*ry2ev,SD.^pwr)
set(gca,'YDir','normal')
xlabel('Wave vector')
ylabel('Energy (eV)')

% -------------------------------------------------------------------------
function KLIST = kgen(KPATH,nk)
% produce a mesh along the k-point path
KLIST = zeros(nk,3);
npath = size(KPATH,1) - 1;
lpath = zeros(npath,1);
for i = 1 : npath
    dk = KPATH(i,:) - KPATH(i+1,:);
    dl = sqrt(sum(dk.*dk));
    if i ~= 1
        lpath(i) = lpath(i-1) + dl;
    else
        lpath(i) = dl;
    end
end
ndivsum = round(nk*lpath/lpath(end));
ndiv(1) = ndivsum(1);
for i = 2 : npath
    ndiv(i) = ndivsum(i) - sum(ndivsum(1:i-1));
end
for i = 1 : npath
    if i == 1
        ni = 1;
    else
        ni = ndivsum(i-1)+1;
    end
    nj = ndivsum(i);
    nk = ndiv(i);
    if i > 1
        ni = ni - 1;
        nk = nk+1;
    end
    for dim = 1 : 3
        KLIST(ni:nj,dim) = ...
            linspace( KPATH(i,dim) , KPATH(i+1,dim) , nk );
    end
end
% -------------------------------------------------------------------------
function [KEIG, EIG, W] = readinput(filename)
% read input data
DATA = importdata(filename);
KEIG = DATA(:,1:3);
EIG = DATA(:,4);
W = DATA(:,5);
% -------------------------------------------------------------------------
function RES = gauss(X,x0,s) % bound. conditions
RES = 1/(s*sqrt(2*pi)) * exp(-((X-x0).^2)/(2*s^2));
% -------------------------------------------------------------------------
function  progrbar( f , fname )
% produces progress bar
persistent flogicold % keep variables in memory between calls
persistent prvcallprogr
persistent fdiap

if isempty(prvcallprogr) % variable is not assigned
    fprintf( 1 , '%s [' , fname );
    prvcallprogr = true;
    fdiap = 0 : 0.1 : 1;
    flogicold = f >= fdiap;
elseif prvcallprogr == false % start a new progress bar
    fprintf( 1 , '%s [' , fname );
    prvcallprogr = true;
    flogicold = f >= fdiap;
elseif prvcallprogr == true
    flogic = f >= fdiap;
    if any( flogicold ~= flogic ) % any change in progress?
        fprintf( 1 , '=' ); % add to the progress bar
        flogicold = flogic; % update the logic array
        if flogic(end) % last done
            fprintf( 1 , ']\n' ); % new line
            prvcallprogr = false; % reset
        end
    end
end