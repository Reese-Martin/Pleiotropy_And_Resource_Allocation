# 
using Statistics
using FileIO
using HypothesisTests
# using CairoMakie #use CairoMakie for final figures because it can save as EPS
using GLMakie #use GLMakie for 'prototyping' ie not the final graphs because of the pop out window and iteration
# GLMakie.activate!(inline=false) # line necessary for GLMakie to generate popup window

#functions
function Jitter(a,Factor = .3)
    a .+= rand(-Factor:.001:Factor,length(a))
end
function AddMedianLines(x,Data,Runs, row, col, Color)
    Rows = Runs
    tmp = median(reshape(Data,(Rows,1)),dims = 1)
    lines!(ga[row,col],x,[tmp[1],tmp[1]], linewidth = 2, color = Color)
end 
function OE_KnockOuts(Hosts)
    Sizes = Hosts[1]
    Nets = Hosts[2]
    KOs = []#store each set of KOs 

    for i in 1:length(Sizes)
        tmp = [] #store each KO from the ith host
        tmp2 = [] #store the network size for each host
        for j in 2:Sizes[i][1][1] #range over signaling proteins in network, skips immune detector
            push!(tmp,Nets[i][1][setdiff(1:end,j),setdiff(1:end,j)])
            push!(tmp2,(Sizes[i][1][1]-1,Sizes[i][1][2]))
        end
        for j in Sizes[i][1][1]+2:Sizes[i][1][1]+Sizes[i][1][2] #range over signaling proteins in network, skips dev detector and the effector
            push!(tmp,Nets[i][1][setdiff(1:end,j),setdiff(1:end,j)])
            push!(tmp2,(Sizes[i][1][1],Sizes[i][1][2]-1))
        end
        push!(KOs,[tmp,tmp2])
    end
    return KOs
end
function OE_OrgLife(Net, Size, Adol, LifeSpan, InfTime, Credits)
    #sim progresses until time steps match LifeSpan
    Net = OE_InfectHost(Net) #Host-parasite Network
    Net1Size = Size[1]
    Net2Size = Size[2]
    
    UseCoef = .01 #amount of protein deactivated for interacting with other proteins
    NumProts = length(Net[1,:]) 
    tmpConc = zeros(LifeSpan,NumProts)
    #set initial protein concentrations 
    tmpConc[1,1:end-1] .= .5

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

    #carry out calculations for LifeSpan number of time steps
    counter = 2
    reducedRange = [1:Net1Size;Net1Size+2:Net1Size+Net2Size+1]
    for i in 1:LifeSpan-1
        if i in InfTime
            tmpConc[i,end] = .5
        end
        tmpConc[i,Net1Size+1] = NSsignal[i]
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
    return tmpConc[:,end-1]
end
function TE_OrgLife(Net, Size, Adol, LifeSpan, InfTime, Credits)
    #sim progresses until time steps match LifeSpan
    Net = TE_InfectHost(Net,Size) #Host-parasite Network
    Net1Size = Size[1]
    Net2Size = Size[2]
    
    UseCoef = .01 #amount of protein deactivated for interacting with other proteins
    NumProts = length(Net[1,:])
    tmpConc = zeros(LifeSpan,NumProts)
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
    
    #carry out calculations for LifeSpan number of time steps
    counter = 2
    reducedRange = [1:Net1Size;Net1Size+3:Net1Size+Net2Size+1]
    for i in 1:LifeSpan-1
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
    return tmpConc[:,Net1Size]
end
function TE_KnockOuts(Hosts)
    Sizes = Hosts[1]
    Nets = Hosts[2]
    KOs = []#store each set of KOs 

    for i in 1:length(Sizes)
        tmp = [] #store each KO from the ith host
        tmp2 = [] #store network sizes for each KO
        for j in 2:Sizes[i][1][1]-1 #range over signaling proteins in network, skips first and last entry
            push!(tmp,Nets[i][1][setdiff(1:end,j),setdiff(1:end,j)])
            push!(tmp2,(Sizes[i][1][1]-1,Sizes[i][1][2]))
        end
        push!(KOs,[tmp,tmp2])
    end
    return KOs
end
function TE_InfectHost(Net,Size) #combines hosts and parasites to create infected hosts
    tmp = vcat(Net[1:Size[1],:],zeros(length(Net[1,:]))',Net[Size[1]+1:end,:])
    tmp[Size[1]+1,1] = 1
    parCol = zeros(length(tmp[:,1]))
    parCol[Size[1] + 1] = .8
    parCol[Size[1]] = -1
    tmp= hcat(tmp[:,1:Size[1]],parCol,tmp[:,Size[1]+1:end]) 
    return tmp
end
function OE_InfectHost(Net) #combines hosts and parasites to create infected hosts
    tmp = vcat(Net,zeros(length(Net[1,:]))')
    tmp[end,1] = 1
    parCol = zeros(length(tmp[:,1]))
    parCol[end] = .8
    parCol[end-1] = -1
    tmp = hcat(tmp,parCol) 
    return tmp
end
function AddSigIndicator(x,Tests,row,col,pos)
    for j in 1:length(x)
        if x[j] < .05/Tests
            text!(ga[row,col],pos[1],pos[2],text = "*",fontsize = 60)
        end
    end
end
function TE_EffectorCons(NetSizes,Nets)
    tmp = zeros(length(NetSizes))
    for i in 1:length(NetSizes)
        tmp[i] = count(Nets[i][1:NetSizes[i][1],NetSizes[i][1]] .!= 0)
    end
    return tmp
end
function OE_EffectorCons(NetSizes,Nets)
    tmp = zeros(length(NetSizes))
    for i in 1:length(NetSizes)
        tmp[i] = count(Nets[i][:,end] .!= 0)
    end
    return tmp
end
#region load data
#topDir = "E:/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Two_Network_Model/"
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Two_Network_Model/"
Data = string(topDir,"Data/")
cd(Data)
Files = filter(x -> endswith(x, ".jld2"), readdir())
Files = filter(x -> !contains(x, "reverse"), Files)
TNTEFiles = filter(x->contains(x, "TNTE"), Files)
TNTEFiles = filter(x -> contains(x, "TNTE"), TNTEFiles)
AsynchTNTEFiles = filter(x -> contains(x, "Asynch"), TNTEFiles)
AsynchTNTEFiles = filter(x -> contains(x, "Net"), AsynchTNTEFiles)
AsynchTNTE_Net_Files = filter(x -> contains(x, "Nets"), AsynchTNTEFiles)[[2,3,1]]
AsynchTNTE_NetSizes_Files = filter(x -> contains(x, "NetSize"), AsynchTNTEFiles)[[2,3,1]]

SynchTNTEFiles = filter(x -> contains(x, "Synch"), TNTEFiles)
SynchTNTEFiles = filter(x -> contains(x, "Net"), SynchTNTEFiles)
SynchTNTE_Net_Files = filter(x -> contains(x, "Nets"), SynchTNTEFiles)[[2,3,1]]
SynchTNTE_NetSizes_Files = filter(x -> contains(x, "NetSize"), SynchTNTEFiles)[[2,3,1]]

TNOEFiles = filter(x -> contains(x, "TNOE"), Files)
AsynchTNOEFiles = filter(x -> contains(x, "Asynch"), TNOEFiles)
AsynchTNOEFiles = filter(x -> contains(x, "Net"), AsynchTNOEFiles)
AsynchTNOE_Net_Files = filter(x -> contains(x, "Nets"), AsynchTNOEFiles)[[2,3,1]]
AsynchTNOE_NetSizes_Files = filter(x -> contains(x, "NetSizes"), AsynchTNOEFiles)[[2,3,1]]

SynchTNOEFiles = filter(x -> contains(x, "Synch"), TNOEFiles)
SynchTNOEFiles = filter(x -> contains(x, "Net"), SynchTNOEFiles)
SynchTNOE_Net_Files = filter(x -> contains(x, "Nets"), SynchTNOEFiles)[[2,3,1]]
SynchTNOE_NetSizes_Files = filter(x -> contains(x, "NetSizes"), SynchTNOEFiles)[[2,3,1]]

SynchTNTE_Bal_Nets = []
SynchTNOE_Bal_Nets = []

AsynchTNTE_Bal_Nets = []
AsynchTNOE_Bal_Nets = []

for i in 1:3
    push!(SynchTNTE_Bal_Nets, [FileIO.load(string(Data,SynchTNTE_NetSizes_Files[i]))["NetSizes"],FileIO.load(string(Data,SynchTNTE_Net_Files[i]))["Nets"]])
    push!(SynchTNOE_Bal_Nets, [FileIO.load(string(Data,SynchTNOE_NetSizes_Files[i]))["NetSizes"],FileIO.load(string(Data,SynchTNOE_Net_Files[i]))["Nets"]])
    
    push!(AsynchTNTE_Bal_Nets, [FileIO.load(string(Data,AsynchTNTE_NetSizes_Files[i]))["NetSizes"],FileIO.load(string(Data,AsynchTNTE_Net_Files[i]))["Nets"]])
    push!(AsynchTNOE_Bal_Nets, [FileIO.load(string(Data,AsynchTNOE_NetSizes_Files[i]))["NetSizes"],FileIO.load(string(Data,AsynchTNOE_Net_Files[i]))["Nets"]])
end
#endregion

#region Generate KOs
SynchTNTE_Bal_KOs = []
SynchTNOE_Bal_KOs = []

AsynchTNTE_Bal_KOs = []
AsynchTNOE_Bal_KOs = []

for i in 1:3
    push!(SynchTNTE_Bal_KOs,TE_KnockOuts(SynchTNTE_Bal_Nets[i]))
    push!(SynchTNOE_Bal_KOs,OE_KnockOuts(SynchTNOE_Bal_Nets[i]))

    push!(AsynchTNTE_Bal_KOs,TE_KnockOuts(AsynchTNTE_Bal_Nets[i]))
    push!(AsynchTNOE_Bal_KOs,OE_KnockOuts(AsynchTNOE_Bal_Nets[i]))
end 
#endregion

#region Evaluate KOs
Adol = 50 
LifeSpan = 150

Credits = [.1,1,1]

SynchTE_Bal_Diffs = []
SynchOE_Bal_Diffs = []

AsynchTE_Bal_Diffs = []
AsynchOE_Bal_Diffs = []

avgSynchTE_Bal_Diffs = []
avgSynchOE_Bal_Diffs = []

avgAsynchTE_Bal_Diffs = []
avgAsynchOE_Bal_Diffs = []
for i in 1:3
    SynchTE_Bal_tmp = []
    SynchOE_Bal_tmp = []

    AsynchTE_Bal_tmp = []
    AsynchOE_Bal_tmp = []

    for j in 1:1000
        avgSynchTE_Bal_tmp = []
        avgSynchOE_Bal_tmp = []
    
        avgAsynchTE_Bal_tmp = []
        avgAsynchOE_Bal_tmp = []
        InfTime = 25
        SynchTNTE_Bal_IntactEffector = TE_OrgLife(SynchTNTE_Bal_Nets[i][2][j][1], SynchTNTE_Bal_Nets[i][1][j][1], Adol, LifeSpan, InfTime, Credits[i])
        for k in 1:length(SynchTNTE_Bal_KOs[i][j][1])
            SynchTNTE_Bal_KOEffector = TE_OrgLife(SynchTNTE_Bal_KOs[i][j][1][k],SynchTNTE_Bal_KOs[i][j][2][k], Adol, LifeSpan, InfTime, Credits[i])
            push!(SynchTE_Bal_tmp,mean(abs.(SynchTNTE_Bal_IntactEffector .- SynchTNTE_Bal_KOEffector)))
            push!(avgSynchTE_Bal_tmp,mean(abs.(SynchTNTE_Bal_IntactEffector .- SynchTNTE_Bal_KOEffector)))
        end
        push!(avgSynchTE_Bal_Diffs,mean(avgSynchTE_Bal_tmp))
        SynchTNOE_Bal_IntactEffector = OE_OrgLife(SynchTNOE_Bal_Nets[i][2][j][1], SynchTNOE_Bal_Nets[i][1][j][1], Adol, LifeSpan, InfTime, Credits[i])
        for k in 1:length(SynchTNOE_Bal_KOs[i][j][1])
            SynchTNOE_Bal_KOEffector = OE_OrgLife(SynchTNOE_Bal_KOs[i][j][1][k],SynchTNOE_Bal_KOs[i][j][2][k], Adol, LifeSpan, InfTime, Credits[i])
            push!(SynchOE_Bal_tmp,mean(abs.(SynchTNOE_Bal_IntactEffector .- SynchTNOE_Bal_KOEffector)))
            push!(avgSynchOE_Bal_tmp,mean(abs.(SynchTNOE_Bal_IntactEffector .- SynchTNOE_Bal_KOEffector)))
        end
        push!(avgSynchOE_Bal_Diffs,mean(avgSynchOE_Bal_tmp))
        InfTime = 75
        AsynchTNTE_Bal_IntactEffector = TE_OrgLife(AsynchTNTE_Bal_Nets[i][2][j][1], AsynchTNTE_Bal_Nets[i][1][j][1], Adol, LifeSpan, InfTime, Credits[i])
        for k in 1:length(AsynchTNTE_Bal_KOs[i][j][1])
            AsynchTNTE_Bal_KOEffector = TE_OrgLife(AsynchTNTE_Bal_KOs[i][j][1][k],AsynchTNTE_Bal_KOs[i][j][2][k], Adol, LifeSpan, InfTime, Credits[i])
            push!(AsynchTE_Bal_tmp,mean(abs.(AsynchTNTE_Bal_IntactEffector .- AsynchTNTE_Bal_KOEffector)))
            push!(avgAsynchTE_Bal_tmp,mean(abs.(AsynchTNTE_Bal_IntactEffector .- AsynchTNTE_Bal_KOEffector)))
        end
        push!(avgAsynchTE_Bal_Diffs,mean(avgAsynchTE_Bal_tmp))

        AsynchTNOE_Bal_IntactEffector = OE_OrgLife(AsynchTNOE_Bal_Nets[i][2][j][1], AsynchTNOE_Bal_Nets[i][1][j][1], Adol, LifeSpan, InfTime, Credits[i])
        for k in 1:length(AsynchTNOE_Bal_KOs[i][j][1])
            AsynchTNOE_Bal_KOEffector = OE_OrgLife(AsynchTNOE_Bal_KOs[i][j][1][k],AsynchTNOE_Bal_KOs[i][j][2][k], Adol, LifeSpan, InfTime, Credits[i])
            push!(AsynchOE_Bal_tmp,mean(abs.(AsynchTNOE_Bal_IntactEffector .- AsynchTNOE_Bal_KOEffector)))
            push!(avgAsynchOE_Bal_tmp,mean(abs.(AsynchTNOE_Bal_IntactEffector .- AsynchTNOE_Bal_KOEffector)))
        end
        push!(avgAsynchOE_Bal_Diffs,mean(avgAsynchOE_Bal_tmp))

    end
    push!(SynchTE_Bal_Diffs, SynchTE_Bal_tmp)
    push!(SynchOE_Bal_Diffs, SynchOE_Bal_tmp)

    push!(AsynchTE_Bal_Diffs, AsynchTE_Bal_tmp)
    push!(AsynchOE_Bal_Diffs, AsynchOE_Bal_tmp)
end
#endregion

#region Plot robustness for the balanced fitness case comparing synchronous and asycnhronous 4 groups per panel, 3 panels
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (nothing,nothing,0,1)
ga = f[1, 1] = GridLayout()
xlabels = ["Independent Effector", "\n Synchronous", "Shared Effector","Independent Effector", "\n Asynchronous", "Shared Effector"]
ga[1, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Scarce Resources \n Knockout Divergence",limits = lims,xgridvisible = false,ygridvisible = false)
ga[2, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Alternating Resources \n Knockout Divergence",limits = lims,xgridvisible = false,ygridvisible = false)
ga[3, 1] = Axis(f,xticks = ([1,1.125,1.25,1.5,1.625,1.75], xlabels),ylabel = "Plentiful Resources \n Knockout Divergence",limits = lims,xgridvisible = false,ygridvisible = false)
KruskalWallis = zeros(3)
Bal_MannU= zeros(3,6)
fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)
for i in 1:3
    SynchTNTE_Bal = Float64.(SynchTE_Bal_Diffs[i])
    SynchTNOE_Bal = Float64.(SynchOE_Bal_Diffs[i])
    AsynchTNTE_Bal = Float64.(AsynchTE_Bal_Diffs[i])
    AsynchTNOE_Bal = Float64.(AsynchOE_Bal_Diffs[i])

    KruskalWallis[i] = pvalue(KruskalWallisTest(SynchTNTE_Bal,SynchTNOE_Bal,AsynchTNTE_Bal,AsynchTNOE_Bal))
    if KruskalWallis[i] < .05
        group = [SynchTNTE_Bal,SynchTNOE_Bal,AsynchTNTE_Bal,AsynchTNOE_Bal]
        tests = [[1,4],[1,3],[1,2],[2,4],[2,3],[3,4]]
        for j in 1:6
            Bal_MannU[i,j] = pvalue(MannWhitneyUTest(group[tests[j][1]],group[tests[j][2]]))
        end
    end
    
    boxplot!(ga[i,1],vec(ones(length(SynchTNTE_Bal))), Float64.(SynchTNTE_Bal), color = :red, width = .25)
    boxplot!(ga[i,1],vec(ones(length(SynchTNOE_Bal))).+.25, Float64.(SynchTNOE_Bal), color = :black, width = .25)
    AddMedianLines([1.15,1.35],SynchTNOE_Bal,length(SynchTNOE_Bal),i,1,:red)

    boxplot!(ga[i,1],vec(ones(length(AsynchTNTE_Bal))).+.5, Float64.(AsynchTNTE_Bal), color = :red, width = .25)
    boxplot!(ga[i,1],vec(ones(length(AsynchTNOE_Bal))).+.75, Float64.(AsynchTNOE_Bal), color = :black, width = .25)
    AddMedianLines([1.65,1.85],AsynchTNOE_Bal,length(AsynchTNOE_Bal),i,1,:red)
end 
#  significance lettering was determined from Bal_MannU and added manually
text!(ga[1,1],string("A"),position = (1.05,.75))
text!(ga[1,1],string("B"),position = (1.3,.75))
text!(ga[1,1],string("C"),position = (1.55,.75))
text!(ga[1,1],string("D"),position = (1.8,.75))
text!(ga[2,1],string("A"),position = (1.05,.75))
text!(ga[2,1],string("B"),position = (1.3,.75))
text!(ga[2,1],string("C"),position = (1.55,.75))
text!(ga[2,1],string("C"),position = (1.8,.75))
text!(ga[3,1],string("A"),position = (1.05,.75))
text!(ga[3,1],string("B"),position = (1.3,.75))
text!(ga[3,1],string("C"),position = (1.55,.75))
text!(ga[3,1],string("D"),position = (1.8,.75))
cd(string(topDir,"/Images"))
# save(string("Evolved_Host_Robustness_Bal_Fit.png"), f, pt_per_unit = 1)
save("F3_Evolved_Host_Robustness_Bal_Fit.pdf", f)
#endregion

#region Figure S9 Plot robustness against the number of connections to the effector
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (.5,16,-.01,.7)
ga = f[1, 1] = GridLayout()
ga[1, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Synchronous Signaling \n Knockout Divergence",title = "Scarce Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[1, 2] = Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = "Plentiful Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[1, 3] = Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = "Varied Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[2, 1] = Axis(f,xticks = ([1,5,10,15],["1","5","10","15"]),ylabel = "Asynchronous Signaling \n Knockout Divergence",limits = lims,xgridvisible = false,ygridvisible = false)
ga[2, 2] = Axis(f,xticks = ([1,5,10,15],["1","5","10","15"]),xlabel = "Connections to Effector",yticklabelsvisible = false,yticksvisible = false,limits = lims,xgridvisible = false,ygridvisible = false)
ga[2, 3] = Axis(f,xticks = ([1,5,10,15],["1","5","10","15"]),yticklabelsvisible = false,yticksvisible = false,limits = lims,xgridvisible = false,ygridvisible = false)

Bal_TTests = zeros(2,3)
for i in 1:3
    SynchTNTE_Bal = Float64.(avgSynchTE_Bal_Diffs[1+(1000*(i-1)):(1000*i)])
    SynchTNTE_Bal_Count = TE_EffectorCons(getindex.(SynchTNTE_Bal_Nets[i][1],1),getindex.(SynchTNTE_Bal_Nets[i][2],1))
    SynchTNOE_Bal = Float64.(avgSynchOE_Bal_Diffs[1+(1000*(i-1)):(1000*i)])
    SynchTNOE_Bal_Count = OE_EffectorCons(getindex.(SynchTNOE_Bal_Nets[i][1],1),getindex.(SynchTNOE_Bal_Nets[i][2],1))
    AsynchTNTE_Bal = Float64.(avgAsynchTE_Bal_Diffs[1+(1000*(i-1)):(1000*i)])
    AsynchTNTE_Bal_Count = TE_EffectorCons(getindex.(AsynchTNTE_Bal_Nets[i][1],1),getindex.(AsynchTNTE_Bal_Nets[i][2],1))
    AsynchTNOE_Bal = Float64.(avgAsynchOE_Bal_Diffs[1+(1000*(i-1)):(1000*i)])
    AsynchTNOE_Bal_Count = OE_EffectorCons(getindex.(AsynchTNOE_Bal_Nets[i][1],1),getindex.(AsynchTNOE_Bal_Nets[i][2],1))

    scatter!(ga[1,i],SynchTNTE_Bal_Count, Float64.(SynchTNTE_Bal), color = :red)
    scatter!(ga[1,i],SynchTNOE_Bal_Count, Float64.(SynchTNOE_Bal), color = :black)

    scatter!(ga[2,i],AsynchTNTE_Bal_Count, Float64.(AsynchTNTE_Bal), color = :red)
    scatter!(ga[2,i],AsynchTNOE_Bal_Count, Float64.(AsynchTNOE_Bal), color = :black)
end 
fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)
cd(string(topDir,"/Images"))
# save(string("Evolved_Host_Robustness_Bal_Fit.png"), f, pt_per_unit = 1)
save("S10_Robustness_vs_effector_connections.pdf", f)
#endregion
