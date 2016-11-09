%% Volcanic Eruptions - No Drag, Iterating with time
clc,clear

%Changing Parameters:
v0=40; %initial velocity in meters per second
w=1000; %width in meters
N=1000; %number of particles
nphases=500; %amount of phases
% uniform_a = 180*rand(1,N); %Uniform random numbers between 0 and 180
% centereduniform_a = 90*rand(1,N) + 45; %Uniform random numbers between 45 and 135
% normal_a = 5*randn(1,N) + 90; %Normal random numbers with mean=90 and sd=5
% Normal random numbers with two means, mean=70° & 110° and Std. Dev.=5:
%  a1=5*randn(1,N/2)+70;
%  a2=5*randn(1,N/2)+110;
%  a=[a1 a2];

%Constants
gftime=1000000; %grain flow time
y0=0; %initial y position is zero
g=9.8; %modified gravity
gs=0.1; %grain size in meters
r=w/2; %horizontal distance (radius from the cone)
x0=r; %crater is in the middle
y=0; %startxing y position
n=0; %n=number of iterations
x(w)=0; %1 meter bin values with length = width (initial topogrpahy surface)
binsize=1; %1 meter bin sizes
hc=0.6; %minimum criteria for grain flow
kmu=0.05; %coefficient of kinetic friction
smu=0.62; %coefficient of static frition


figure;

%Part 1 - Ballistic (Particle landing):
%Preliminary Layers:
while x~=hc
    a1=5*randn(1,N/2)+80;
    a2=5*randn(1,N/2)+100;
    a=[a1 a2]; %angle distribution
    v0x=v0.*cosd(a); %initial x-velocity
    v0y=v0.*sind(a); %initial y-velocity

    y=@(t) y0+v0y.*t-0.5*g.*t.^2; %vertical particle path
    ypath=zeros(5/0.5,N); %makes empty matrix as placeholder

    for t=0:0.5:50 %looping with time with 0.5 second increments
        n=n+1; %counts number of loops
        ypath(n,:)=y(t);
    end

    %Finds the time each particle hits the ground
    yindex = zeros(size(ypath(1,:)));
    for p=1:N
        yindex(p)=min(find(ypath(:,p)<0));
        groundtime(p)=ypath(yindex(p),p);
    end
    tg=yindex/2; %time particle hits the ground. Time before the y-distance becomes negative
    xland=x0+v0x.*tg; %finds the horizontal landing distance from crater
    xlandbin=round(xland); %+1 %integer number for horizontal landing distance. Index in bin

%         figure
%         topo=histogram(a,N); %histogram of angles (for last layer...)


    [occurrences,unique_val]=hist(xlandbin,unique(xlandbin)); %finds any repitition of landing location
    x(unique_val)=x(unique_val)+(occurrences.*gs); %counts particles that fall into same bin
    n=0; %resets loops for every layer
end


%Part 2: Grain Flow & Successive Layers in time steps:

%Morphing (downhill only... - temp)
%disp(x) %original topography
xflow=0;
dt=0.5; %delta t (time step increments)
for t2=1:dt:gftime
    gfabs=find(x>0); %indeces of bins that only contain particles
    gxindex=min(gfabs):max(gfabs); %array of topography with containing only particles
    gpoints=find(abs(x(gxindex)-x(gxindex+1)>=hc) | abs(x(gxindex)-x(gxindex-1)>=hc)); %particles that meet the criteria for grain flow and that will move
    gpoints=gpoints+min(gfabs)-1;
    if ~isempty(gpoints)

        for gfcrit=gpoints

            %Heights are equal on both sides of particle
            if x(gfcrit)-x(gfcrit+1) == x(gfcrit)-x(gfcrit-1)
                gfval=x(gfcrit)-x(gfcrit+1) == x(gfcrit)-x(gfcrit-1);
                
                %Particle lands on right side of crater and moves right
                if gfcrit > r
                    gfcond=gfcrit > r;
                    aflow=atand((x(gfcrit)-x(gfcrit+1))/binsize).*gfval.*gfcond; %angle of moving particle at surface with the width=binsize
                    xflow=round(dt^2*g*(sind(aflow)-kmu*cosd(aflow))); %distance particle moves within the time step
                    %disp('Right side of crater and moves right')
                    if gfcrit+xflow>0
                        x(gfcrit)=x(gfcrit)-gs; %removes particle from original landing
                        x(gfcrit+xflow)=x(gfcrit+xflow)+gs; %adds new particle location to topography array
                    end
                    
                %Particle lands on the left side of the crater and moves left   
                elseif gfcrit < r
                    gfcond=gfcrit < r;
                    aflow=atand((x(gfcrit)-x(gfcrit-1))/binsize).*gfval.*gfcond;
                    xflow=-1*round(dt^2*g*(sind(aflow)-kmu*cosd(aflow))); %particle moves left
                    %disp('Left side of crater and moves left')
                    if gfcrit+xflow>0
                        x(gfcrit)=x(gfcrit)-gs; %removes particle from original landing
                        x(gfcrit+xflow)=x(gfcrit+xflow)+gs; %adds new particle location to topography array
                    end
                
                %Particle lands in crater and moves left or right randomly
                else x(gfcrit)-x(gfcrit+1) == x(gfcrit)-x(gfcrit-1);
                    aflow=atand((x(gfcrit)-x(gfcrit+1))/binsize).*gfval;
                    xflow=randi([-1 1])*round(dt^2*g*(sind(aflow)-kmu*cosd(aflow)));
                    %disp('Inside of crater and moves left or right')
                    if gfcrit+xflow>0
                        x(gfcrit)=x(gfcrit)-gs; %removes particle from original landing
                        x(gfcrit+xflow)=x(gfcrit+xflow)+gs; %adds new particle location to topography array
                    end
                end
                
                
            %Compares heights to the right of each bin
            elseif x(gfcrit)-x(gfcrit+1)>=hc
                gfval=x(gfcrit)-x(gfcrit+1)>=hc;
                aflow=atand((x(gfcrit)-x(gfcrit+1))/binsize).*gfval; %angle of moving particle at surface with the width=binsize
                xflow=round(dt^2*g*(sind(aflow)-kmu*cosd(aflow))); %distance particle moves within the time step
                %disp('grain flow right')
                if gfcrit+xflow>0
                    x(gfcrit)=x(gfcrit)-gs; %removes particle from original landing
                    x(gfcrit+xflow)=x(gfcrit+xflow)+gs; %adds new particle location to topography array
                end

            %Compares heights to the left of each bin
            elseif x(gfcrit)-x(gfcrit-1)>=hc
                gfval=x(gfcrit)-x(gfcrit-1)>=hc;
                aflow=atand((x(gfcrit)-x(gfcrit-1))/binsize).*gfval;
                xflow=-1*round(dt^2*g*(sind(aflow)-kmu*cosd(aflow)));
                %disp('grain flow left')
                if gfcrit+xflow>0
                    x(gfcrit)=x(gfcrit)-gs; %removes particle from original landing
                    x(gfcrit+xflow)=x(gfcrit+xflow)+gs; %adds new particle location to topography array
                end
                
            end
        end
    end
end


%Phases for Successive Layers:
for phases=1:nphases
        
    %Back to Part 1:

    a1=5*randn(1,N/2)+80;
    a2=5*randn(1,N/2)+100;
    a=[a1 a2]; %angle distribution
    v0x=v0*cosd(a); %initial x-velocity
    v0y=v0*sind(a); %initial y-velocity

    y=@(t) y0+v0y.*t-0.5*g.*t.^2; %vertical particle path
    ypath=zeros(5/0.5,N); %makes empty matrix as placeholder

    for t=0:0.5:50 %looping with time with 0.5 second increments
        n=n+1; %counts number of loops
        ypath(n,:)=y(t);
    end

    %Finds the time each particle hits the ground
    yindex = zeros(size(ypath(1,:)));
    for p=1:N
        yindex(p)=min(find(ypath(:,p)<0));
        groundtime(p)=ypath(yindex(p),p);
    end
    tg=yindex/2; %time particle hits the ground. Time before the y-distance becomes negative
    xland=x0+v0x.*tg; %finds the horizontal landing distance from crater
    xlandbin=round(xland); %+1 %integer number for horizontal landing distance. Index in bin

%     figure
%     topo=histogram(a,w); %histogram of angles (for last layer...)


    [occurrences,unique_val]=hist(xlandbin,unique(xlandbin)); %finds any repitition of landing location
    x(unique_val)=x(unique_val)+(occurrences.*gs); %counts particles that fall into same bin
    n=0;

        
    %Part 2 - Grain Flow    
    for t2=1:dt:gftime
        gfabs=find(x>0); %indeces of bins that only contain particles
        gxindex=min(gfabs):max(gfabs); %array of topography with containing only particles
        gpoints=find( abs(x(gxindex)-x(gxindex+1)>=hc) | abs(x(gxindex)-x(gxindex-1)>=hc)); %particles that meet the criteria for grain flow and that will move
        gpoints=gpoints+min(gfabs)-1;
        if ~isempty(gpoints)

            for gfcrit=gpoints

                %Heights are equal on both sides of particle
                if x(gfcrit)-x(gfcrit+1) == x(gfcrit)-x(gfcrit-1)
                    gfval=x(gfcrit)-x(gfcrit+1) == x(gfcrit)-x(gfcrit-1);

                    %Particle lands on right side of crater and moves right
                    if gfcrit > r
                        gfcond=gfcrit > r;
                        aflow=atand((x(gfcrit)-x(gfcrit+1))/binsize).*gfval.*gfcond; %angle of moving particle at surface with the width=binsize
                        xflow=round(dt^2*g*(sind(aflow)-kmu*cosd(aflow))); %distance particle moves within the time step
                        %disp('Right side of crater and moves right')
                        if gfcrit+xflow>0
                            x(gfcrit)=x(gfcrit)-gs; %removes particle from original landing
                            x(gfcrit+xflow)=x(gfcrit+xflow)+gs; %adds new particle location to topography array
                        end

                    %Particle lands on the left side of the crater and moves left   
                    elseif gfcrit < r
                        gfcond=gfcrit < r;
                        aflow=atand((x(gfcrit)-x(gfcrit-1))/binsize).*gfval.*gfcond;
                        xflow=-1*round(dt^2*g*(sind(aflow)-kmu*cosd(aflow))); %particle moves left
                        %disp('Left side of crater and moves left')
                        if gfcrit+xflow>0
                            x(gfcrit)=x(gfcrit)-gs; %removes particle from original landing
                            x(gfcrit+xflow)=x(gfcrit+xflow)+gs; %adds new particle location to topography array
                        end

                    %Particle lands in crater and moves left or right randomly
                    else x(gfcrit)-x(gfcrit+1) == x(gfcrit)-x(gfcrit-1);
                        aflow=atand((x(gfcrit)-x(gfcrit+1))/binsize).*gfval;
                        xflow=randi([-1 1])*round(dt^2*g*(sind(aflow)-kmu*cosd(aflow)));
                        %disp('Inside of crater and moves left or right')
                        if gfcrit+xflow>0
                            x(gfcrit)=x(gfcrit)-gs; %removes particle from original landing
                            x(gfcrit+xflow)=x(gfcrit+xflow)+gs; %adds new particle location to topography array
                        end
                    end

                    
                %Compares heights to the right of each bin
                elseif x(gfcrit)-x(gfcrit+1)>=hc
                    gfval=x(gfcrit)-x(gfcrit+1)>=hc;
                    aflow=atand((x(gfcrit)-x(gfcrit+1))/binsize).*gfval; %angle of moving particle at surface with the width=binsize
                    xflow=round(dt^2*g*(sind(aflow)-kmu*cosd(aflow))); %distance particle moves within the time step
                    %disp('grain flow right')
                    if gfcrit+xflow>0
                        x(gfcrit)=x(gfcrit)-gs; %removes particle from original landing
                        x(gfcrit+xflow)=x(gfcrit+xflow)+gs; %adds new particle location to topography array
                    end

                    
                %Compares heights to the left of each bin
                elseif x(gfcrit)-x(gfcrit-1)>=hc
                    gfval=x(gfcrit)-x(gfcrit-1)>=hc;
                    aflow=atand((x(gfcrit)-x(gfcrit-1))/binsize).*gfval;
                    xflow=-1*round(dt^2*g*(sind(aflow)-kmu*cosd(aflow)));
                    %disp('grain flow left')
                    if gfcrit+xflow>0
                        x(gfcrit)=x(gfcrit)-gs; %removes particle from original landing
                        x(gfcrit+xflow)=x(gfcrit+xflow)+gs; %adds new particle location to topography array
                    end

                end
            end
        end
    end

    %Graph of volcano
    if rem(phases,50)==0
        topo=plot(x)
    end
    hold on
    xlabel('Topography Width (meters)')
    ylabel('Height (meters)')
    axis equal
    %axis([0 w 0 w])
    %axis auto
    %axis image
end
%% Volcanic Eruptions - No Drag, looping with height
clc,clear

x0=0; y0=0; %initial x and y positions are zero
v0=10; %initial velocity
g=9.8; %modified gravity
th=10; %thickness in meters
w=40; %width in meters
gs=1; %grain size in meters
h=0; %accumulating height
r=w/2; %horizontal distance (radius from the cone)
N=0; %N=number of particles
x(r)=0; %1 meter bin values with the width of radius

while h<th
    
%     for a = 90.*rand(1,1); %Random number generator for angles between 0 and 90
%         nbins = 1;
%         hist(a,nbins)
%     end

    %a = 90.*rand(1,1); %Random number generator for angles between 0 and 90
    a=45;

    v0x=v0*sind(a); %initial x-velocity
    v0y=v0*cosd(a); %initial y-velocity
    
    syms tg
    f = @(tg) v0y*tg-4.9*tg^2; %particle's vertical path
    tg=solve(f,tg); %tg=time particle hits the ground
    tg=double(tg(2));

    xloc = x0+v0x.*tg; %horizontal location of particle after landing
    xbin=round(xloc);
    x(xbin)=gs; %modify to save last iteration here!!!

    N=N+1;
    h=h+1
end
