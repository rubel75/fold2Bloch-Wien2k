function ubs_dots
% Plot undolded band structure
%
% Update history:
% - added periodic images of the BZ (May 26, 2020)
%
% (c) Oleg Rubel, McMaster University

%% Init. parameters
KPATH = [0 1/2 0; ...
        0 0 0; ...
        0 0 1/2]; % k-point path
Dp2s = [1 0 -1
        0 2 0
        1 0 1]; % transformation matrix used to transform a primitive cell to a supercell
KLABEL = {'Y'; 'G'; 'Z'};
finpt = '6-atom2D.f2b'; % input file name
Ef = 0.0460363511; % Fermi energy (Ry)
ERANGE = [Ef-1 Ef+0.4]; % energy range for plot (Ry)
ry2ev = 13.605698066; % Ry -> eV conversion factor
pwr = 1/1; % power for result plotting
         % 1 - linear scale, 1/2 - sqrt, etc.
         % 0 - folded bands (needs wth = 0)
msz = 10; % marker size for plot
lwdth = 0.5; % plot line width
fontSize = 9; % points
PLTSZ = [1 1 600/1.5 300/1.5]; % plot size
wth = 0.05; % threshold weight
clrmp = jet;    % flipud(gray)
                % flipud(pink)
                % flipud(hot)
                % flipud(autumn)
                % cool
                % flipud(bone)
                % flipud(jet)
                % jet
G = [ 0.333333  0.000000  0.000000;
      0.000000  0.166667  0.000000;
      0.000000  0.000000  0.111111]; % Reciprocal latt. vect. from *.outputkgen
roundOffErrK = 0.000001; % this is the round off error 1/3 = 0.333333 + err


%% INITIALIZATION
[KEIG, EIG, W] = readinput(finpt); % read input data from file
% EIG - energy eigenvalues
% KEIG - k-list for eigenvalues
% W - list of characters

%% Convert energy units [Ry] -> [eV]
EIG = EIG*ry2ev;
Ef = Ef*ry2ev;
ERANGE = ERANGE*ry2ev;

%% MAIN
L = [];
ENE = [];
WGHT = [];
G = G'; % transpose G matrix (need for Wien2k)
for i=1 : 3
    G(i,:)=Dp2s*transpose(G(i,:)); % rescale reciprocal lattice vectors 
end                                % from supercell to primitive cell
dl = 0; % cumulative length of the path
KPATH = coordTransform(KPATH,G);
KEIG = coordTransform(KEIG,G);
epsk = [roundOffErrK roundOffErrK roundOffErrK]; % k rounding error
epsk = coordTransform(epsk,G); % transform to Cart. coords
epsk = sqrt(dot(epsk,epsk)); % get magnitude of the vector
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
                        if dist < epsk % k-point is on the path
                            A = KPATH(ikp,:) - KEIG(j,:);
                            x = dot(A,B)/dk;
                            if x >= 0  &&  x <= dk+epsk % k-point is within the path range
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
scatter(L,(ENE-Ef), WGHTRS*msz, WGHTRS,'LineWidth',lwdth);
hold on;
axis([XTICKS(1) XTICKS(end) ERANGE(1)-Ef ERANGE(2)-Ef])
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
scatter(DATX,DAT, DATRS*msz, DATRS,'LineWidth',lwdth);
caxis([0 1])
ylabel('Spectral weight')

% SAVE plot as *.eps
print( [finpt '.eps'], '-depsc')

% -------------------------------------------------------------------------
function W = coordTransform(V,G)
% transform vector V(:,3) in G(3,3) coord. system -> W(:,3) in Cartesian coordinates
% G vector elements are in columns!
W = zeros(size(V));
for i = 1:size(V,1)
    W(i,:) = G(1,:)*V(i,1) + G(2,:)*V(i,2) + G(3,:)*V(i,3);
end
% -------------------------------------------------------------------------
function WRESCL = rescale(W,pwr)
% rescale weights using a power functio W^pwr
WRESCL=W.^(pwr); % rescale if needed to enhance
WRESCL = WRESCL + eps; % need eps to make plot "heapy"
% -------------------------------------------------------------------------
function [KEIG, EIG, W] = readinput(filename)
% read input data
DATA = importdata(filename);
KEIG = DATA(:,1:3);
EIG = DATA(:,4);
W = DATA(:,5);
% -------------------------------------------------------------------------
function RES = dp2l(X0,X1,X2,epsk) % distance from point {X0} to line {X1}-{X2}
% see http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
denom = X2 - X1;
denomabs = sqrt(dot(denom,denom));
if denomabs < eps
    display(X1); display(X2);
    error('X1 = X2');
end
numer = cross( X0-X1 , X0-X2 );
numerabs = sqrt(dot(numer,numer));
RES = numerabs/denomabs;
if RES <= epsk
    if norm(X0-X1) + norm(X0-X2) - norm(X1-X2) >  2*epsk
        RES = Inf; % exclude collinear points that are not on the segment
    end
end
% -------------------------------------------------------------------------
