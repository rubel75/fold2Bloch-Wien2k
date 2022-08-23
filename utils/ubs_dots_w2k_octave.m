% Plot undolded band structure. The plot generated is not as nice as 
% the Matlab plot. The script was developed for tutorial purposes and 
% is not actively supported. For publication-quality plots use 
% Matlab ubs_dots.m
%
% Usage:
%
%    first luch Octave GUI
%    $ octave
%
%    next within Octave execute the script
%    >> ubs_dots_w2k_octave
%
% Update history:
% - ported from Matlab (Sept 12, 2017)
% - updated for compatibility with Octave 4.0.x (Aug 3, 2019)
% - added periodic images of the BZ (May 26, 2020)
%
% (c) Oleg Rubel, McMaster University

function ubs_dots_w2k_octave

%% Init. parameters
KPATH = [1/2 0 0;...
         0 0 0; ...
         1/2 1/2 0]; % k-point path
Dp2s = [1 0 -1
        0 2 0
        1 0 1]; % transformation matrix used to transform a primitive cell to a supercell
KLABEL = {'L'; 'G'; 'X'};
finpt = 'SiGe.f2b'; % input file name
Ef = 0.385799; % Fermi energy (Ry)
ERANGE = [Ef-0.2 Ef+0.4]; % energy range for plot (Ry)
ry2ev = 13.605698066; % Ry -> eV conversion factor
pwr = 1/1; % power for result plotting
         % 1 - linear scale, 1/2 - sqrt, etc.
         % 0 - folded bands (needs wth = 0)
msz = 5; % marker size for plot
lwdth = 0.5; % plot line width
fontSize = 9; % points
PLTSZ = [1 1 600 400]; % plot size
wth = 0.05; % threshold weight
clrmp = jet;    % flipud(gray)
                % flipud(pink)
                % flipud(hot)
                % flipud(autumn)
                % cool
                % flipud(bone)
                % flipud(jet)
                % jet
G = [   0.083726 -0.027909 -0.027909;
   0.000000  0.078938 -0.039469;
   0.000000  0.000000  0.068362];  % Reciprocal latt. vect. from *.outputkgen
G = G'; % transpose G matrix (Wien2k only!)

%% INITIALIZATION
%[KEIG, EIG, W] = readinput(finpt); % read input data from file
S = load("-ascii",finpt);
KEIG = S(:,1:3); % KEIG - k-list for eigenvalues
EIG = S(:,4); % EIG - energy eigenvalues
W = S(:,5); % W - list of characters
clear S % deallocate S

%% MAIN
L = [];
ENE = [];
WGHT = [];
G=Dp2s*G; % rescale reciprocal lattice vectors 
          % from supercell to primitive cell
dl = 0; % cumulative length of the path
KPATH = coordTransform(KPATH,G);
KEIG = coordTransform(KEIG,G);
XTICKS = [0];
for ikp = 1 : size(KPATH,1)-1
    B = KPATH(ikp,:) - KPATH(ikp+1,:);
    dk = sqrt(dot(B,B));
    XTICKS = [XTICKS; XTICKS(ikp)+dk];
    for j = 1 : length(EIG)
        if EIG(j) > ERANGE(1) && EIG(j) < ERANGE(2) && W(j) >= wth
            dist = Inf; % initialize distance to a path
            quit = false;
            for ikx = -1:1 % include periodic images of the BZ
                for iky = -1:1
                    for ikz = -1:1
                        KPERIOD = [ikx iky ikz]; % periodic shift
                        % transform to Cartezian coords
                        KPERIOD = coordTransform(KPERIOD,G);
                        % evaluate distance to the path
                        dist2 = dp2l( KEIG(j,:) + KPERIOD , ...
                            KPATH(ikp,:) , KPATH(ikp+1,:) , epsk );
                        % select smallest distance
                        dist = min(dist,dist2);
                        if dist < eps % k-point is on the path
                            A = KPATH(ikp,:) - KEIG(j,:) - KPERIOD;
                            x = dot(A,B)/dk;
                            if x >= 0  &&  x <= dk+eps % k-point is within the path range
                                L = [L; x+dl]; % append k-point coordinate along the path
                                ENE = [ENE; EIG(j)]; % append energy list
                                WGHT = [WGHT; W(j)];
                            end
                            quit = true;
                            break; % exit ikz loop
                        end
                    end % ikz loop
                    if quit
                        break; % exit iky loop
                    end
                end % iky loop
                if quit
                    break; % exit ikx loop
                end
            end % ikx loop
        end
    end
    dl = dl + dk;
end
if isempty(L)
    msg = ['No eigenvalues are selected for the plot. ', ...
        'The likely reason is that the energy range is ', ...
        'too restrictive (check ERANGE), or no k-points are located ', ...
        'on the path selected (check KPATH)'];
    error(msg);
end

%% Plot results
hFig = figure(1);

% Fig 1(a)
subplot(1,2,1);
set(gca,'FontSize',fontSize);
set(hFig, 'Position', PLTSZ, 'PaperPositionMode','auto')
map = colormap(clrmp);
WGHTRS = rescale(WGHT,pwr);
%scatter(L,(ENE-Ef)*ry2ev, WGHTRS*msz, WGHTRS,'LineWidth',lwdth); % compatible with octave 4.2.1
scatter(L,(ENE-Ef)*ry2ev, WGHTRS*msz, WGHTRS); % comaptible with octave 4.0.3
hold on;
axis([XTICKS(1) XTICKS(end) min((ENE-Ef)*ry2ev) max((ENE-Ef)*ry2ev)])
yticks = get(gca,'ytick');
set(gca,'YTick',yticks);
for i = 1 : length(yticks)
    newYTick{i} = sprintf('%1.1f',yticks(i));
end
set(gca,'YTickLabel',newYTick);
hline = plot([0 XTICKS(end)],[0 0]); % Fermi level
set(hline,'Color','k','LineStyle','--');
set(gca,'XTick',XTICKS);
set(gca,'XTickLabel',KLABEL);
set(gca,'XGrid','on', 'GridLineStyle','-');
caxis([0 1]); % normalize intensities to 1
xlabel('Wave vector')
ylabel('Energy (eV)')
box on
hold off

% Fig 1(b)
subplot(1,2,2);
set(gca,'FontSize',fontSize);
DAT = linspace(0,1,10);
DATX = ones(size(DAT));
DATRS = rescale(DAT,pwr);
%scatter(DATX,DAT, DATRS*msz, DATRS,'LineWidth',lwdth); % compatible with octave 4.2.1
scatter(DATX,DAT, DATRS*msz, DATRS); % comaptible with octave 4.0.3
caxis([0 1])
ylabel('Spectral weight')

%pause % hold on before exiting

endfunction
% -------------------------------------------------------------------------
function W = coordTransform(V,G)
% transform vector V(:,3) in G(3,3) coord. system -> W(:,3) in Cartesian coordinates
% G vector elements are in columns!
W = zeros(size(V));
for i = 1:size(V,1)
    W(i,:) = G(1,:)*V(i,1) + G(2,:)*V(i,2) + G(3,:)*V(i,3);
end;
endfunction
% -------------------------------------------------------------------------
function WRESCL = rescale(W,pwr)
% rescale weights using a power functio W^pwr
WRESCL=W.^(pwr); % rescale if needed to enhance
WRESCL = WRESCL + eps; % need eps to make plot "heapy"
endfunction
% -------------------------------------------------------------------------
function [KEIG, EIG, W] = readinput(filename)
% read input data
DATA = importdata(filename);
KEIG = DATA(:,1:3);
EIG = DATA(:,4);
W = DATA(:,5);
endfunction
% -------------------------------------------------------------------------
function RES = dp2l(X0,X1,X2,epsk) % distance from point {X0} to line {X1}-{X2}
% see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
denom = X2 - X1;
denomabs = sqrt(dot(denom,denom));
if denomabs < eps
    display(X1); display(X2);
    error('X1 = X2');
end;
numer = cross( X0-X1 , X0-X2 );
numerabs = sqrt(dot(numer,numer));
RES = numerabs/denomabs;
if RES <= epsk
    if norm(X0-X1) + norm(X0-X2) - norm(X1-X2) >  2*epsk
        RES = Inf; % exclude collinear points that are not on the segment
    end;
end;
endfunction
% -------------------------------------------------------------------------
