%% Traveling network simulation
% Nate Cira 3/12/2024
%% Begin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

%% Establish parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%branching angle (input in degrees, converted to radians)
theta=deg2rad(60);

%equilibrium size (sum of edge lengths)
eqS=100;

%growth rate (normalized by retraction time scale)
kG=.2;

%branching rate (normalized by retraction time scale)
kB=.1;

%retraction rate (normalized to 1)
kR=1;

%switching rate (normalized by retraction time scale, and equal to kB for steady state)
kS=kB;

%hill coefficient
n=10;

%minimum number of free leaves
fleafmin=2;

% discretization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%distance step of branching, growth, retraction
dd=1;

% #time steps per retraction time (higher value indicates more time steps per retraction time scale)
%below sets an optimal dt
dt=max(3*kB+kG,1);

%total time
runtime=400;

%probabilities
PB = kB/dt;
PG = kG/dt;
PS = kS/dt;
PR = kR/dt;

% Display related %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display the simulation? if no, 0, if yes, 1, if every x frames, >1
display=1;

%Move display frame with the network? if no, 0, if yes, 1
moveframe=1;

%display size
xdisp = 80;
ydisp = 80;

% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize two vertices
% vertices each occupy a row in the variable "vertices"
% their data structure is identity, partner1, partner2, partner3, xcoord,
% ycoord, retraction state
vertices=[1 2 0 0 0 0+dd/2 0; 2 1 0 0 0 -dd/2 0];

% namecounter keeps track of which vertex identities have been used
namecounter=3;

%current network size
S=dd;

%% Run the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:ceil(runtime*dt)

    elims=[];
    newret=0;

    %find the leaves
    leaves=find((vertices(:,3)==0 & vertices(:,4)==0) | (vertices(:,2)==0 & vertices(:,3)==0) | (vertices(:,2)==0 & vertices(:,4)==0));

    %count how many retracting leaves
    retleaves=find(vertices(:,7)==1);
    numret=length(retleaves);

    %manipulate each leaf
    for x=1:length(leaves)

        %find vertices which are not retracting
        if vertices(leaves(x),7)==0

            %modify the switching probability by size
            mPS=2*PS*1/(1+(eqS/(S))^n);

            %roll dice
            dice=rand;

            %% Branching %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if dice<PB

                %find the coordinates of the partner to establish an angle
                coords = vertices(find(vertices(:,1)==vertices(leaves(x),2)),5:6);

                %calculate an angle between the two points
                angle = atan((vertices(leaves(x),6)-coords(2))/(vertices(leaves(x),5)-coords(1)));

                if (vertices(leaves(x),6)-coords(2)>=0 && vertices(leaves(x),5)-coords(1)<0) ||(vertices(leaves(x),6)-coords(2)<0 && vertices(leaves(x),5)-coords(1)<0)
                    angle=angle+pi;
                end

                %calculate the new vertices' positions
                v1c = [dd*cos(theta/2)*cos(angle)-dd*sin(theta/2)*sin(angle) dd*cos(theta/2)*sin(angle)+dd*sin(theta/2)*cos(angle)]+vertices(leaves(x),5:6);
                v2c = [dd*cos(theta/2)*cos(angle)+dd*sin(theta/2)*sin(angle) dd*cos(theta/2)*sin(angle)-dd*sin(theta/2)*cos(angle)]+vertices(leaves(x),5:6);

                %establish the new vertices
                vertices=[vertices; namecounter, vertices(leaves(x),1),0,0,v1c(1),v1c(2),0;namecounter+1, vertices(leaves(x),1),0,0,v2c(1),v2c(2),0];

                %establish partners in the branching vertex
                vertices(leaves(x),3:4)=[namecounter,namecounter+1];

                namecounter=namecounter+2;
                S=S+2*dd;

            %% Growth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif dice<PB+PG

                %find the coordinates of the partner to establish an angle
                coords = vertices(find(vertices(:,1)==vertices(leaves(x),2)),5:6);

                %calculate an angle between the two points
                angle = atan((vertices(leaves(x),6)-coords(2))/(vertices(leaves(x),5)-coords(1)));

                if (vertices(leaves(x),6)-coords(2)>=0 && vertices(leaves(x),5)-coords(1)<0) ||(vertices(leaves(x),6)-coords(2)<0 && vertices(leaves(x),5)-coords(1)<0)
                    angle=angle+pi;
                end

                %Move the vertex
                vertices(leaves(x),5:6)=[dd*cos(angle) dd*sin(angle)]+vertices(leaves(x),5:6);
                S=S+dd;

            %% Switching %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif dice<PB+PG+mPS && length(leaves)-numret-newret>fleafmin

                newret=newret+1; % do not allow switching if it would result in less than the minimum number of free leaves
                vertices(leaves(x),7)=1;

            %% Do nothing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif dice<PB+PG+mPS && length(leaves)-numret-newret<=fleafmin
            else
            end

        %% retract if already retracting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            dice=rand;
            if dice<PR
                %find the coordinates of the partner to establish an angle
                partnerposition = find(vertices(leaves(x),2:4)~=0)+1;
                partnerid=find(vertices(:,1) == vertices(leaves(x),partnerposition));

                coords = vertices(partnerid,5:6);

                %calculate an angle between the two points
                angle = atan((vertices(leaves(x),6)-coords(2))/(vertices(leaves(x),5)-coords(1)));

                if (vertices(leaves(x),6)-coords(2)>=0 && vertices(leaves(x),5)-coords(1)<0) ||(vertices(leaves(x),6)-coords(2)<0 && vertices(leaves(x),5)-coords(1)<0)
                    angle=angle+pi;
                end
                %Move the vertex
                vertices(leaves(x),5:6)=[-dd*cos(angle) -dd*sin(angle)]+vertices(leaves(x),5:6);
                S=S-dd;

                %remove the vertex and switch the retraction counter (if a new leaf is created) once fully retracted
                if round(5*(coords(1)-vertices(leaves(x),5)))==0 && round(5*(coords(2)-vertices(leaves(x),6)))==0
                    elims = [elims; leaves(x)];

                    for p=2:4
                        if vertices(partnerid,p)==vertices(leaves(x),1)
                            vertices(partnerid,p)=0;
                        end
                    end

                    if sum(vertices(partnerid,2:4)==0)==2
                        vertices(partnerid,7)=1;
                    else
                    end
                end
            end
        end

    end

    % eliminate the removed vertices
    for e=1:length(elims)
        vertices(elims(e),:)=[];
        elims=elims-1;
    end

    %% display network
    if rem(t-1,display)==0
        %    render image
        clf
        hold on
        for v=1:size(vertices,1)
            for p=2:4
                if vertices(v,p)~=0
                    for v2=1:size(vertices,1)
                        if vertices(v2,1)==vertices(v,p)
                            plot ([vertices(v,5) vertices(v2,5)], [vertices(v,6) vertices(v2,6)],'k')
                        end
                    end
                end
            end
        end

        if moveframe==1
            %calculate centers
            avgx=mean(vertices(:,5));
            avgy=mean(vertices(:,6));
            axis([avgx-xdisp/2 avgx+xdisp/2 avgy-ydisp/2 avgy+ydisp/2])
            axis square
        else
            axis([-xdisp/2 xdisp/2 -ydisp/2 ydisp/2])
            axis square
        end
        hold off
        pause(.05)

    end

end
