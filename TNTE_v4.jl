#region Documentation
# Hosts start with two effectors and fitness is determined in a balanced, immune dependent, or dev dependent way
# v4 added a batch functionality to the code to allow for greater than 100 run sims
using Core: GeneratedFunctionStub
using Base: RangeStepIrregular, end_base_include, Float64, StatusActive

#Import packages
using Random
using Trapz
using Dates
using JLD2
using FileIO
using StatsBase
using Statistics

#Define Agent Structures, which contain the specific information for each agent necessary for simulation
mutable struct Organism
    Network :: Matrix{Float64}
    Concentrations :: Matrix{Float64}
    Net1_Size :: Int64
    Net2_Size :: Int64
    Lineage :: Int64
end
mutable struct Parasite
    Network :: Array{Float64}
    InfTime :: Int64
    Concentrations :: Vector{Float64}
end

#region Declare constant variables, some variables are defined as constants in this region for performance needs
const NumHosts = 500 #declares number of hosts to be used in each generation
const Runs = 100 #number of runs the sim goes through
const CredLimits = [.1,1]#[.1,1,Inf] limits how much protein can be activated in each time step of a simulation
const workDir = pwd() #gets current directory for saving data
#endregion

#region Funtions
function culpoint(Host,a,NumHosts)# function for finding if a host should be killed
    return count(x->x < a[Host],a)/NumHosts
end
function Simulation(NumHosts,NumPars,Generations,Credits,OrgFitBalance,Timing,Start_Time,Run)# function that executes a single simulation consisting of "Generations" number of generations 
    DeathCoef = .3 #amount of parasites and hosts that die each generation
    DeathThreshold = .9 #if the parasite area of infection exceeds this value, hosts die 
    Adol = 50 #number of Timesteps taken for adolescence
    LifeSpan = 150 #number of timesteps taken for full life
    gen = 1

    #define host and parasite arrays, populate with randomly generated individuals
    Hosts = []
    Parasites = []
    BlankParasite = Parasite(zeros(5),0,[.5])
    for i in 1:NumHosts
        tmp = rand(0:1,(5,5))
        tmp[[1,5,21,25]].= 0
        tmp1 = (rand(Float64,(5,5))*2).-1
        tmpImmNet = tmp.*tmp1
        
        tmp2 = rand(0:1,(5,5))
        tmp2[[1,5,21,25]].= 0
        tmp3 = (rand(Float64,(5,5))*2).-1
        tmpNSNet = tmp2.*tmp3
        
        zer = zeros(5,5)
        Nets = hcat(vcat(tmpImmNet,zer),vcat(zer,tmpNSNet))

        push!(Hosts,Organism(Nets,[0 0 0 0 0 0 0 0 0 0; .5 .5 .5 .5 .5 .5 .5 .5 .5 .5],5,5,i))
        if i <= NumPars
            partmp = zeros(5)
            partmp[rand(2:4)] = rand()*2-1
            if sum(partmp) == 0
                partmp[2] = .5
            end
            push!(Parasites,Parasite(partmp,0,[0]))
        end
    end

    #define arrays for storing population information from the simulation
    HostFit = zeros(Float64, Generations, NumHosts); #host population fitness
    ParFit = zeros(Float64, Generations, NumPars); #parasite population fitness
    ParBurden = zeros(Float64,Generations); 
    genParBurden = zeros(Float64,NumPars);
    PleioFitness = [] #fitness of pleiotropic hosts
    NonPleioFitness = [] #fitness of non-pleiotropic hosts
    PerPopPleio = zeros(Float64,Generations) # percent of the population that is pleiotropic 

    while gen <= Generations
        phosts = findPleioHosts((p -> p.Net1_Size).(Hosts),(p -> p.Network).(Hosts),500) #find all pleiotropic hosts in generation
        PerPopPleio[gen] = sum(phosts)

        #carry out host and parasite dynamics for the generation
        for i in 1:NumHosts
            if i <= NumPars
                #randomly generate an infection start time, eithe concurrent with dev or during adulthood 
                if Timing == "Synch"
                    ParIntro = rand(2:Adol)
                elseif Timing == "Asynch"
                    ParIntro = rand(Adol+1:LifeSpan - 50)
                end
                Parasites[i].InfTime = ParIntro
                Hosts[i].Concentrations,Parasites[i].Concentrations = OrgLife(Hosts[i],Parasites[i],Adol,LifeSpan,ParIntro,Credits) #rand infection occurs duing generic signaling
            else 
                Hosts[i].Concentrations,~ = OrgLife(Hosts[i], BlankParasite,Adol,LifeSpan,[],Credits)
            end
        end 
        #Calculate host and parasite Fitness
        HostToPop = []
        HostDead = 0
        for i in 1:NumHosts
            if i <= NumPars
                #calc the area under parasite infection
                ParCon = Parasites[i].Concentrations[Parasites[i].InfTime:end]
                ParArea = trapz(1:length(ParCon),ParCon)/trapz(1:length(ParCon),ones(length(ParCon)))

                #calc the area under immune effector 
                HostEffCon = Hosts[i].Concentrations[:,Hosts[i].Net1_Size]
                HostEffArea = trapz(1:length(HostEffCon),HostEffCon)/trapz(1:length(HostEffCon),ones(length(HostEffCon)))

                #calc match to dev signal, and kill hosts where parasites exceed deaththreshold
                GenMatch = NSMatch(Hosts[i].Concentrations,Hosts[i].Net1_Size, Adol)
                if ParArea > DeathThreshold
                    push!(HostToPop,i)
                    HostDead = HostDead + 1
                end

                ParFit[gen,i] = ParArea
                HostFit[gen,i],genParBurden[i] = OrgFitness(HostEffArea, OrgFitBalance, ParArea, GenMatch)
            else
                #no infection so norm area becomes 0, and post inf eq is the same as the pre-inf
                HostEffCon = Hosts[i].Concentrations[:,Hosts[i].Net1_Size]
                HostEffArea = trapz(1:length(HostEffCon),HostEffCon)/trapz(1:length(HostEffCon),ones(length(HostEffCon)))

                #calc match to dev signal
                GenMatch = NSMatch(Hosts[i].Concentrations,Hosts[i].Net1_Size, Adol)

                #calc host fitness when there is no infection
                HostFit[gen,i],~  = OrgFitness(HostEffArea, OrgFitBalance, 0, GenMatch)
            end
        end

        #populate pleio/non-pleio host fitness arrays
        tmp = []
        tmp2 = []
        for ii in 1:NumHosts
            if phosts[ii] == 1
                push!(tmp,HostFit[gen,ii])
            else
                push!(tmp2,HostFit[gen,ii])
            end
        end
        push!(PleioFitness, tmp)
        push!(NonPleioFitness, tmp2)
        ParBurden[gen] = mean(genParBurden)
        #if HostDead exceeds 30% of hosts then remove hosts from dead population at random
        while HostDead > NumHosts*DeathCoef
            rem = rand(1:HostDead)
            deleteat!(HostToPop,rem)
            HostDead -= 1 
        end

        #Cull the populations in a fitness weighted manner
        HostDeaths(HostDead,HostFit[gen,:],HostToPop,NumHosts,DeathCoef,Hosts)

        tmpParFit = ParFit[gen,:]
        ProgNum = zeros(NumPars)
        for i in 1:NumPars
            if tmpParFit[i] < .33
                ProgNum[i] = 1
            elseif .34 < tmpParFit[i] < .66
                ProgNum[i] = 2
            else
                ProgNum[i] = 3
            end 
        end
        
        p = sortperm(vec(tmpParFit), rev = true)
        Parasites = Parasites[p]
        ProgNum = ProgNum[p]
        Parasites = Parasites[1:NumPars - Int64(NumPars*DeathCoef)]          
        DeadHosts = NumHosts-length(Hosts)         

        #repopulate Hosts and parasitse allowing for evolution
        OrgEvolution(Hosts,HostFit[gen,:],DeadHosts)
        MeanSigs = floor(Statistics.mean([length(x.Network[1,:]) for x in Hosts]))-4 #minus four to account for detector and effector of both networks
        ParEvolution(Parasites, ProgNum, NumPars*DeathCoef, Int64(MeanSigs))
        
        #randomize order of host and parasite populations
        HostShuffle = Random.randperm(NumHosts)  
        Hosts = Hosts[HostShuffle]
        for i in 1:length(Hosts)
            Hosts[i].Concentrations = zeros(2, length(Hosts[i].Network[1,:]))
            Hosts[i].Concentrations[2,:] .+= .5
        end
        
        ParShuffle = Random.randperm(NumPars)
        Parasites = Parasites[ParShuffle]
        if gen%100 == 0
            print('\n',Run, ' ', gen,' ',(Dates.Time(Dates.now()) - Start_Time)/Nanosecond(1)* (1/1000000000))
        end
        
        gen = gen+1
    end
    return (p -> p.Network).(Hosts), collect(zip((p -> p.Net1_Size).(Hosts),(p -> p.Net2_Size).(Hosts))), HostFit, (p -> p.Network).(Parasites), ParFit, ParBurden, PerPopPleio, PleioFitness, NonPleioFitness
end
function HostDeaths(HostDead,HostFit,HostToPop,NumHosts,DeathCoef,Hosts)# kills specified number of hosts in a population
    culled = false
    Host=1
    while !culled
        if HostDead >= DeathCoef*NumHosts
            culled = true
        elseif Host == NumHosts
            culled = true
        end
        if !culled 
            Culpoint = culpoint(Host,HostFit,NumHosts)

            if HostFit[Host] < 1e-3
                push!(HostToPop,Host)
                HostDead += 1  
            elseif rand() > Culpoint
                push!(HostToPop,Host)
                HostDead += 1  
            end
            
            Host += 1
        end
    end

    while HostDead > NumHosts*DeathCoef
        rem = rand(1:HostDead)
        deleteat!(HostToPop,rem)
        HostDead -= 1 
    end
    HostToPop = sort(HostToPop, rev = true)
    for rem in HostToPop
        deleteat!(Hosts, rem)
    end
end
function OrgLife(Host, Para, Adol, LifeSpan, InfTime, Credits)# function that carries out simulations of host life
    #sim progresses until time steps match LifeSpan
    Net = InfectHost(Para, Host) #Host-parasite Network
    Net1Size = Host.Net1_Size
    Net2Size = Host.Net2_Size
    
    UseCoef = .01 #amount of protein deactivated for interacting with other proteins
    NumProts = length(Net[1,:])
    tmpConc = zeros(LifeSpan+5,NumProts) #+5 for 5 timesteps of 'burn-in' to allow the networks to erach Equilibrium
    #set initial protein concentrations
    tmpConc[1,1:Net1Size] .= .5
    tmpConc[1,Net1Size+2:end] .= .5

    #generate vector of regulatory coefficients and interaction number
    AllAct = Vector{Vector{Float64}}(undef,NumProts)
    AllUse = Vector{Int64}(undef,NumProts)
    for i in 1:NumProts
        AllAct[i] = Net[:,i]
        AllUse[i] = sum(Net[i,:].!= 0)
    end
    xs = range(0, 8pi, length= LifeSpan)
    NSsignal = (sin.(xs)./2).+ .5
    NSsignal[Adol+1:end] .= 0
    NSsignal = [zeros(5);NSsignal]
    #carry out calculations for LifeSpan number of time steps
    counter = 2
    reducedRange = [1:Net1Size;Net1Size+3:Net1Size+Net2Size+1]
    for i in 1:LifeSpan+5-1 #-1 because we don't need to simulate the timestep after te end of the host lifespan
        if i in InfTime
            tmpConc[i,Net1Size+1] = .5
        end
        tmpConc[i,end-(Net2Size-1)] = NSsignal[i]
        tmpCred = 0
        for prot in 1:NumProts #collect up and down regulating actions for each protein, sum them, and then add to the previous time steps active protein concentration
            InitProtConc = tmpConc[i,prot]
            ActingProts = AllAct[prot]
            tmpConc[counter,prot] = InitProtConc - (UseCoef*AllUse[prot]) 
            for j in 1:NumProts
                v = ActingProts[j]
                if v > 0 
                    tmpConc[counter,prot] = tmpConc[counter,prot] + (v*tmpConc[i,j])*(1-InitProtConc)
                 elseif v < 0
                    tmpConc[counter,prot] = tmpConc[counter,prot] + (v*tmpConc[i,j])*(InitProtConc)
                 end 
            end
            if tmpConc[counter,prot] > 1
                tmpConc[counter,prot] = 1
            elseif tmpConc[counter,prot] < 0
                tmpConc[counter,prot] = 0
            end
            if (prot in reducedRange) & (tmpConc[counter,prot] - tmpConc[i,prot] > 0) #skip dev signal and parasite
                tmpCred = tmpCred + (tmpConc[counter,prot] - tmpConc[i,prot])
            end
        end

        if tmpCred > Credits
            scale = Credits/tmpCred
            for prot in reducedRange
                tmp2 = tmpConc[counter,prot] - tmpConc[i,prot]
                if tmp2>0
                    tmpConc[counter,prot] = tmpConc[i,prot] + ((tmpConc[counter, prot] - tmpConc[i,prot])*scale)
                end
            end
        end

        counter += 1
    end 
    return tmpConc[6:end,:], tmpConc[6:end,Net1Size+1]
end
function InfectHost(Para, Host)# combines hosts and parasites to create infected hosts
    HostNet = Host.Network
    ParNet = Para.Network
    ParTarg = findfirst(x -> x != 0,ParNet)

    if (ParTarg in 2:Host.Net1_Size-1) || (ParTarg in Host.Net1_Size+2:Host.Net1_Size+Host.Net2_Size-1) #insures a signaling protein is being targeted
        tmp = vcat(HostNet[1:Host.Net1_Size,:],zeros(length(HostNet[1,:]))',HostNet[Host.Net1_Size+1:end,:])
        tmp[Host.Net1_Size+1,ParTarg] = ParNet[ParTarg]
        tmp[Host.Net1_Size+1,1] = 1

        parCol = zeros(length(tmp[:,1]))
        parCol[Host.Net1_Size + 1] = .8
        parCol[Host.Net1_Size] = -1
        tmp = hcat(tmp[:,1:Host.Net1_Size],parCol,tmp[:,Host.Net1_Size+1:end]) 
    else #if the parasite is not targeting one of the signaling proteins, then it does not get to effect hosts networks at all. Consider this a biological incompatibility
        tmp = vcat(HostNet[1:Host.Net1_Size,:],zeros(length(HostNet[1,:]))',HostNet[Host.Net1_Size+1:end,:])
        tmp[Host.Net1_Size+1,1] = 1

        parCol = zeros(length(tmp[:,1]))
        parCol[Host.Net1_Size + 1] = .8
        parCol[Host.Net1_Size] = -1
        tmp = hcat(tmp[:,1:Host.Net1_Size],parCol,tmp[:,Host.Net1_Size+1:end]) 
    end
    return tmp
end
function OrgEvolution(Hosts,HostFit,ToRep)# repopulates host population with chance for mutation in offspring
    NewOrgs = []
    while length(NewOrgs) < ToRep #run through survivors and potentially mutate
        RepHost = rand(1:length(Hosts))
        RepChance = count(x->x <= HostFit[RepHost],HostFit)/NumHosts
        RepCheck = rand()
        if (RepChance > RepCheck) 
            Dupe = deepcopy(Hosts[RepHost])
            ToMute = rand()

            if ToMute < 5e-3 #if the mut threshold is passed, go on to mutations
                Mutation = rand()
                if Mutation <= .25 #add edge to graph Lineage code: 1
                    Net1_len = Dupe.Net1_Size
                    Net2_len = Dupe.Net2_Size
                    Zers = findall(x->x==0,Dupe.Network)
                    #delete the zeros that would lead to detector/effector self regulation as well as D-E connections
                    deleteat!(Zers,findall(x->x == CartesianIndex(1,1), Zers)[1])
                    deleteat!(Zers,findall(x->x == CartesianIndex(1,Net1_len), Zers)[1])
                    deleteat!(Zers,findall(x->x == CartesianIndex(Net1_len,1), Zers)[1])
                    deleteat!(Zers,findall(x->x == CartesianIndex(Net1_len,Net1_len), Zers)[1])

                    #delete the zeros that would lead to detector/effector self regulation as well as D-E connections for d2 e2
                    deleteat!(Zers,findall(x->x == CartesianIndex(Net1_len+1,Net1_len+1), Zers)[1])
                    deleteat!(Zers,findall(x->x == CartesianIndex(Net1_len+1,Net1_len+Net2_len), Zers)[1])
                    deleteat!(Zers,findall(x->x == CartesianIndex(Net1_len+Net2_len,Net1_len+1), Zers)[1])
                    deleteat!(Zers,findall(x->x == CartesianIndex(Net1_len+Net2_len,Net1_len+Net2_len), Zers)[1])

                    #delete zeros that would allow d2<->e1 or d1<->e2
                    deleteat!(Zers,findall(x->x == CartesianIndex(1,Net1_len+Net2_len), Zers)[1])
                    deleteat!(Zers,findall(x->x == CartesianIndex(Net1_len+Net2_len,1), Zers)[1])
                    deleteat!(Zers,findall(x->x == CartesianIndex(Net1_len+1,Net1_len), Zers)[1])
                    deleteat!(Zers,findall(x->x == CartesianIndex(Net1_len,Net1_len+1), Zers)[1])
                    
                    if length(Zers) >0
                        AddEdge = Zers[rand(1:length(Zers))]
                        Dupe.Network[AddEdge[1],AddEdge[2]] = rand().*2 .-1
                    end
                    
                elseif (Mutation > .25) & (Mutation <= .5) #delete edge from graph
                    Ones = findall(x -> x!=0,Dupe.Network)
                    if length(Ones) >1
                        DelEdge = Ones[rand(1:length(Ones))]
                        Dupe.Network[DelEdge[1],DelEdge[2]] = 0
                    end 
                elseif (Mutation > .5) & (Mutation <= .8) #change coefficient by 10% randomly up or down
                    Coefs = findall(x -> x!=0, Dupe.Network)
                    if length(Coefs) > 0
                        ChanCoef = Coefs[rand(1:length(Coefs))]
                        delta = Dupe.Network[ChanCoef[1],ChanCoef[2]]*.1
                        tmp = rand()
                        if tmp > .5
                            Dupe.Network[ChanCoef[1],ChanCoef[2]] = Dupe.Network[ChanCoef[1],ChanCoef[2]]+delta
                        else
                            Dupe.Network[ChanCoef[1],ChanCoef[2]] = Dupe.Network[ChanCoef[1],ChanCoef[2]]-delta
                        end
                        if Dupe.Network[ChanCoef[1],ChanCoef[2]] > 1
                            Dupe.Network[ChanCoef[1],ChanCoef[2]] = 1
                        elseif Dupe.Network[ChanCoef[1],ChanCoef[2]] < -1
                            Dupe.Network[ChanCoef[1],ChanCoef[2]] = -1
                        end
                    end
                elseif (Mutation >.8) & (Mutation <= .9) #delete a protein from the network, reworked to insure that both networks maintain at least one signaling protein
                    if (Dupe.Net1_Size > 3) & (Dupe.Net2_Size > 3) #both bigger than 3 then choose from either net
                        ToDel = rand([2:Dupe.Net1_Size-1;Dupe.Net1_Size+2:Dupe.Net1_Size+Dupe.Net2_Size-1])
                    elseif (Dupe.Net1_Size > 3) & (Dupe.Net2_Size <= 3) #net1>3 net2<=3 then choose from 1
                        ToDel = rand(2:Dupe.Net1_Size-1)
                    elseif (Dupe.Net1_Size <= 3) & (Dupe.Net2_Size > 3) #net2>3 net1 <=3 then choose from net 2
                        ToDel = rand(Dupe.Net1_Size+2:Dupe.Net1_Size+Dupe.Net2_Size-1)
                    else #neither bigger so no deletion
                        ToDel = NaN
                    end
                    if ~isnan(ToDel)
                        Dupe.Network = Dupe.Network[1:end .!= ToDel,1:end .!= ToDel]
                        if  2 <= ToDel <= Dupe.Net1_Size-1
                            Dupe.Net1_Size -= 1
                        elseif Dupe.Net1_Size+2 <= ToDel <= Dupe.Net1_Size+Dupe.Net2_Size-1
                            Dupe.Net2_Size -= 1
                        end
                    end
                else #duplicate a protein in the network
                    toDup = rand([2:Dupe.Net1_Size-1;Dupe.Net1_Size+2:Dupe.Net1_Size+Dupe.Net2_Size-1])
                    tmpRow = Dupe.Network[toDup,:]
                    tmpRowStack = vcat(Dupe.Network[1:toDup-1,:],tmpRow',Dupe.Network[toDup:end,:])
                    tmpCol = tmpRowStack[:,toDup]
                    tmpColStack = hcat(tmpRowStack[:,1:toDup-1],tmpCol,tmpRowStack[:,toDup:end])
                    Dupe.Network = tmpColStack
                    
                    if  2 <= toDup <= Dupe.Net1_Size-1
                        Dupe.Net1_Size += 1
                    elseif Dupe.Net1_Size+2 <= toDup <= Dupe.Net1_Size+Dupe.Net2_Size-1
                        Dupe.Net2_Size += 1
                    end
                end
            end
            push!(NewOrgs,Dupe)
        end
    end
    for k in NewOrgs
        push!(Hosts,k)
    end
end 
function OrgFitness(EArea, OrgFitBalance, PArea, GenericMatch)# calculates organismal fitness
    if OrgFitBalance == "Balanced"
        ImmCoef = 2.5
        ParCoef = 2.5
        GenCoef = 2.5
    elseif OrgFitBalance == "GenHeavy"
        ImmCoef = .5
        ParCoef = .5
        GenCoef = 4.5
    elseif OrgFitBalance == "ImmHeavy"
        ImmCoef = 4.5
        ParCoef = 4.5
        GenCoef = .5
    end
    ImmCost = ImmCoef*EArea
    ParCost = ParCoef*PArea
    GenCost = GenCoef*abs(GenericMatch) #multiplies the absolute difference in Non Specific effector by 2 to scale to the same range of damage as immune effector 
    
    return exp(-(ImmCost + ParCost + GenCost )), ParCost
end
function ParEvolution(Parasites, ProgNum, ToRep, MeanSigs)# repopulates parasites with chance for mutation
    NewPars = []
    Rep = false
    par = 1
    while !Rep
        prog = ProgNum[par]
        for Surv in 1:prog
            Dupe = deepcopy(Parasites[par])
            if length(Dupe.Network) < 2+Int64(MeanSigs)
                tmp = zeros(2+MeanSigs)
                tmp[1:length(Dupe.Network)] = Dupe.Network
                Dupe.Network = tmp
            end
        
            if rand() < 1e-2
                Mutation = rand()
                if Mutation < .5
                    if MeanSigs == 1
                        prevTarg = findall(x->x !=0,Dupe.Network)[1]
                        Dupe.Network[2] = Dupe.Network[prevTarg]
                        if prevTarg != 2
                            Dupe.Network[prevTarg] = 0
                        end
                    else
                        prevTarg = findall(x->x !=0,Dupe.Network)[1]
                        NewTarg = Int64(rand(2:MeanSigs+1))
                        Dupe.Network[NewTarg] = Dupe.Network[prevTarg]
                        if NewTarg != prevTarg
                            Dupe.Network[prevTarg] = 0
                        end
                    end
                else
                    prevTarg = findall(x->x !=0,Dupe.Network)[1]
                    Dupe.Network[prevTarg] = rand()*2-1
                end
            end
            push!(NewPars,Dupe)
        end
        par += 1
        if length(NewPars) >= ToRep
            Rep = true
        end
    end
    for i in 1:Int64(ToRep)
        push!(Parasites, NewPars[i])
    end
end
function NSMatch(x,Net1Size,Adol)# function for determining how closely a host matches the dev input
    #Net1Size+2 because the array x includes parasite concentrations between Eff1 and Det2
    NSsignal = x[5:Adol-3,Net1Size+2] 
    #the 3 TS offset accounts for the time it takes to stimulate the network leading to the increase of effector
    NSeff = x[8:Adol,end]

    MeanDiff = mean(abs.(NSsignal.-NSeff))
    Correl = cor(NSsignal,NSeff)
    if isnan(Correl)
        Correl = 0
    end
    Score = (1-Correl) + (MeanDiff) #fitness improves by lowering score
 
    return Score/2.63 #the range for the score is [0,2.63] when correl is 1 and MeanDiff is 0, or Correl is -1 and MeanDiff is ~.63 (ie Effector is the -sin of input)
    #so divide by 2.63 to bring the range down to [0,1] which makes scaled fitness contributions easier to calculate.
end
function findPleioHosts(x,y,Runs)# function to find pleio hosts
    tmp = zeros(Runs) #records whether there was any pleiotropy
    for i in 1:Runs
        ImmToNS = count(y[i][1:x[i][1],x[i][1]+1:end] .!= 0)
        NSToImm = count(y[i][x[i][1]+1:end,1:x[i][1]] .!= 0)
        if (ImmToNS != 0) || (NSToImm != 0)
            tmp[i] = 1
        end
    end
    return tmp
end
#endregion
for cred in CredLimits #cycles through specified infections chances, conducts (Runs) simulations at a given level of infection 
    for batch in 1:10
        #collection of arrays for the saving of data
        Credits = cred
        inf = .5
        OrgFitBalance = "Balanced"
        # OrgFitBalance = "GenHeavy"
        # OrgFitBalance = "ImmHeavy"
        # Timing = "Synch"
        Timing = "Asynch"

        #collection of arrays for the saving of data
        Generations = 1000 #number of generations the sim runs for
        NumPars = Int64(NumHosts*inf); #how many parasites to create in each run
        HostFit = zeros(Float64, Runs,  Generations, NumHosts); #tracks host fitness across Runs, and each generation within a run
        ParasiteBurden = zeros(Float64, Runs, Generations); #tracks the individual fitness components of each host along with the lineage, allows mapping of the transistion through fitness space
        ParFit = zeros(Float64, Runs,  Generations, NumPars); #tracks parasites fitness across Runs, and each generation within a run
        
        PerPopPleio = zeros(Float64, Runs,  Generations); #tracks parasites fitness across Runs, and each generation within a run
        PleioHostFitness = zeros(Float64,Runs,Generations);
        NoPleioHostFitness = zeros(Float64,Runs,Generations);

        MostCommonNet =  [] #array to save the most common host in each simulation
        MostCommonNet_Size = []

        ParTargets = zeros(Float64, Runs,  Generations, NumPars); #what signaling protein the parasites are acting on

        Run = 1
        Start_Time = Dates.Time(Dates.now())
        while Run <= Runs #conducts Runs number of simulations
            #call simulation function
            print('\n',"******** TNTE " ,Timing," Credit limit ",Credits, " Fitness Effects ",OrgFitBalance," Batch ",batch," ********** ", (Dates.Time(Dates.now()) - Start_Time)/Nanosecond(1)* (1/1000000000))
            HostNetworks,NetSizes,HostFit[Run,:,:],ParNetworks,ParFit[Run,:,:],ParasiteBurden[Run,:],PerPopPleio[Run,:],tmpPleioFit,tmpNoPleioFit = Simulation(NumHosts,NumPars,Generations,Credits,OrgFitBalance,Timing,Start_Time,Run);
            for j in 1:length(tmpPleioFit)
                if size(tmpPleioFit[j])[1] == 0
                    tmpPleioFit[j] = NaN
                end
                if size(tmpNoPleioFit[j])[1] == 0
                    tmpNoPleioFit[j] = NaN
                end
            end
            PleioHostFitness[Run,:] = mean.(tmpPleioFit)
            NoPleioHostFitness[Run,:] = mean.(tmpNoPleioFit)

            #identify most common Host and save
            UnqNetCM = countmap(HostNetworks)
            Vals = collect(values(UnqNetCM))
            UnqNets = collect(keys(UnqNetCM))
            UnqNetSizes = [NetSizes[findfirst(x -> x == i, HostNetworks)] for i in UnqNets]
            MostCom = findall(x -> x == maximum(Vals),Vals)
            push!(MostCommonNet, UnqNets[MostCom])
            push!(MostCommonNet_Size, UnqNetSizes[MostCom])
            Run += 1 
            #print progress    
        end

        #saves data where each generation has been averaged (no individual host response, but significantly smaller file sizes)
        HostFit = reshape(mean(HostFit, dims = 3),(Runs,Generations))
        ParFit = reshape(mean(ParFit, dims = 3),(Runs,Generations))

        FileName = string(workDir,"/Data/Batches/",string(inf*100)[1:end-2],"_PercentInf_EC_Limit_",string(Credits*100)[1:end-2],"_FitEff_",OrgFitBalance,"_TNTE_",Timing,"_Batch_",batch,".jld2")
        save(FileName,"ParTarg",ParTargets,"HostFit",HostFit,"ParBurden",ParasiteBurden,"ParFit",ParFit,"PerPleio",PerPopPleio,"PleioHostFitness",PleioHostFitness,"NoPleioHostFitness",NoPleioHostFitness)
        save(string(FileName[1:end-5],"_Networks.jld2"),"Networks",MostCommonNet,"NetSizes", MostCommonNet_Size)
    end
end