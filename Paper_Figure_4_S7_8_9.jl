# generate plots associated with Network features, prevalence of pleiotropy, size, connectivity, etc

using Statistics
using FileIO
using HypothesisTests
using CairoMakie #use CairoMakie for final figures because it can save as EPS
# using GLMakie #use GLMakie for 'prototyping' ie not the final graphs because of the pop out window and iteration
#region functions
function PerPleio(x,y,Runs) #returns the percent of pleiotropic conections in use 
    # x = netsizes, y = networks
    I2NS = zeros(Runs) #raw number of connections from I to NS
    NS2I = zeros(Runs) #raw number of connections from NS to I
    all = zeros(Runs)
    for i in 1:Runs
        I2NS[i] = count(y[i][1:x[i][1],x[i][1]+1:end] .!= 0)/((x[i][1]*x[i][2])-4) #minus 4 excludes potential D-E,D-D,E-E,E-D connections 
        NS2I[i] = count(y[i][x[i][1]+1:end,1:x[i][1]] .!= 0)/((x[i][1]*x[i][2])-4) #minus 4 excludes potential D-E,D-D,E-E,E-D connections 
        all[i] = (count(y[i][1:x[i][1],x[i][1]+1:end] .!= 0) + count(y[i][x[i][1]+1:end,1:x[i][1]] .!= 0))/(((x[i][1]*x[i][2])-4)*2)
    end
    return I2NS, NS2I, all
end
function PerPleio_OE(x,y,Runs) #returns the percent of pleiotropic conections in use 
    # x = netsizes, y = networks
    I2NS = zeros(Runs) #raw number of connections from I to NS
    NS2I = zeros(Runs) #raw number of connections from NS to I
    all = zeros(Runs)
    for i in 1:Runs
        I2NS[i] = count(y[i][1:x[i][1],x[i][1]+1:end-1] .!= 0)/((x[i][1]*x[i][2])-1) #minus 1 excludes potential D-D
        NS2I[i] = count(y[i][x[i][1]+1:end-1,1:x[i][1]] .!= 0)/((x[i][1]*x[i][2])-1) #minus 1 excludes potential D-D
        all[i] = (count(y[i][1:x[i][1],x[i][1]+1:end-1] .!= 0) + count(y[i][x[i][1]+1:end-1,1:x[i][1]] .!= 0))/(((x[i][1]*x[i][2])-1)*2)
    end
    return I2NS, NS2I, all
end
function Jitter(a,Factor = .3)
    a .+= rand(-Factor:.001:Factor,length(a))
end
function AddMedianLines(x,Data,Runs, row, col, Color)
    Rows = Runs
    tmp = median(reshape(Data,(Rows,1)),dims = 1)
    lines!(ga[row,col],x,[tmp[1],tmp[1]], linewidth = 2, color = Color)
end 
function AddSigIndicator(x,Tests,row,col,pos)
    for j in 1:length(x)
        if x[j] < .05/Tests
            text!(ga[row,col],pos[1],pos[2],text = "*",fontsize = 60)
        end
    end
end
function TE_Connectivity(NetSizes,Nets)
    tmp = zeros(length(NetSizes))
    for i in 1:length(Nets)
        a = sum(Nets[i][1:NetSizes[i][1],1:NetSizes[i][1]].!=0)
        b = sum(Nets[i][NetSizes[i][1]+1:NetSizes[i][1]+NetSizes[i][2],NetSizes[i][1]+1:NetSizes[i][1]+NetSizes[i][2]].!=0)
        tmp[i] = (a+b)/(NetSizes[i][1]*NetSizes[i][1]+NetSizes[i][2]*NetSizes[i][2]-8) #minus 8 to account for the 8 D<->E connections that cannot occur
    end
    return tmp
end
function OE_Connectivity(NetSizes,Nets)
    tmp = zeros(length(NetSizes))
    for i in 1:length(Nets)
        a = sum(Nets[i][1:NetSizes[i][1],1:NetSizes[i][1]].!=0)
        b = sum(Nets[i][NetSizes[i][1]+1:NetSizes[i][1]+NetSizes[i][2],NetSizes[i][1]+1:NetSizes[i][1]+NetSizes[i][2]].!=0)
        tmp[i] = (a+b)/(NetSizes[i][1]*NetSizes[i][1]+NetSizes[i][2]*NetSizes[i][2]) #No need to account for the D-E connections as they are not indluded in the slices taken
    end
    return tmp
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
#endregion
#region load data
#topDir = "E:/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Two_Network_Model/"
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Two_Network_Model/"
Data = string(topDir,"Data/")
cd(Data)
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Two_Network_Model/"
Data = string(topDir,"Data/")
cd(Data)
Files = filter(x->endswith(x, ".jld2"), readdir())
Files = filter(x->!contains(x, "reverse"), readdir())
TNTEFiles = filter(x -> contains(x, "TNTE"), Files)
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

SynchTNTE_Bal_NetSizes = []
SynchTNTE_Bal_Nets = []

SynchTNOE_Bal_NetSizes = []
SynchTNOE_Bal_Nets = []

AsynchTNTE_Bal_NetSizes = []
AsynchTNTE_Bal_Nets = []

AsynchTNOE_Bal_NetSizes = []
AsynchTNOE_Bal_Nets = []

for i in 1:3
    push!(SynchTNTE_Bal_NetSizes, getindex.(FileIO.load(string(Data,SynchTNTE_NetSizes_Files[i]))["NetSizes"],1))
    push!(SynchTNTE_Bal_Nets, getindex.(FileIO.load(string(Data,SynchTNTE_Net_Files[i]))["Nets"],1))
    push!(SynchTNOE_Bal_NetSizes, getindex.(FileIO.load(string(Data,SynchTNOE_NetSizes_Files[i]))["NetSizes"],1))
    push!(SynchTNOE_Bal_Nets, getindex.(FileIO.load(string(Data,SynchTNOE_Net_Files[i]))["Nets"],1))
    
    push!(AsynchTNTE_Bal_NetSizes, getindex.(FileIO.load(string(Data,AsynchTNTE_NetSizes_Files[i]))["NetSizes"],1))
    push!(AsynchTNTE_Bal_Nets, getindex.(FileIO.load(string(Data,AsynchTNTE_Net_Files[i]))["Nets"],1))
    push!(AsynchTNOE_Bal_NetSizes, getindex.(FileIO.load(string(Data,AsynchTNOE_NetSizes_Files[i]))["NetSizes"],1))
    push!(AsynchTNOE_Bal_Nets, getindex.(FileIO.load(string(Data,AsynchTNOE_Net_Files[i]))["Nets"],1))
end
#endregion

#region  figure 4 plot showing the percent of the most common hosts in a simulation that were pleiotropic (to any degree)
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1250, 1000))
lims = (nothing,nothing,0,1)
ylabels = ["50% Chance of Infection ","90% Chance of Infection"]
ga = f[1, 1] = GridLayout()
xlabels = ["Independent Effector", "\n Synchronous", "Shared Effector","Independent Effector", "\n Asynchronous", "Shared Effector"]
ga[1, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Scarce Resources", title = "Proportion of Hosts with upstream pleiotropy", limits = lims,xgridvisible = false,ygridvisible = false)
ga[2, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Alternating Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[3, 1] = Axis(f,xticks = ([1,1.125,1.25,1.5,1.625,1.75], xlabels),ylabel = "Plentiful Resources",limits = lims,xgridvisible = false,ygridvisible = false)
Bal_ChiSqTests = zeros(2,3)
fontsize_theme = Theme(fontsize = 24)
set_theme!(fontsize_theme)
for i in 1:3
    ~,~,SynchTNTE_Bal = PerPleio(SynchTNTE_Bal_NetSizes[i],SynchTNTE_Bal_Nets[i],1000)
    ~,~,SynchTNOE_Bal = PerPleio_OE(SynchTNOE_Bal_NetSizes[i],SynchTNOE_Bal_Nets[i],1000)

    ~,~,AsynchTNTE_Bal = PerPleio(AsynchTNTE_Bal_NetSizes[i],AsynchTNTE_Bal_Nets[i],1000)
    ~,~,AsynchTNOE_Bal = PerPleio_OE(AsynchTNOE_Bal_NetSizes[i],AsynchTNOE_Bal_Nets[i],1000)

    a = barplot!(ga[i,1], [1],[count(SynchTNTE_Bal.>0)/1000], color = :red, width = .25)
    b = barplot!(ga[i,1], [1.25],[count(SynchTNOE_Bal.>0)/1000], color = :black, width = .25)
    Bal_ChiSqTests[1,i] = pvalue(ChisqTest([count(SynchTNTE_Bal.>0) count(SynchTNTE_Bal.==0) ; count(SynchTNOE_Bal.>0) count(SynchTNOE_Bal.==0)]))
    
    barplot!(ga[i,1], [1.5],[count(AsynchTNTE_Bal.>0)/1000], color = :red, width = .25)
    barplot!(ga[i,1], [1.75],[count(AsynchTNOE_Bal.>0)/1000], color = :black, width = .25)
    Bal_ChiSqTests[2,i] = pvalue(ChisqTest([count(AsynchTNTE_Bal.>0) count(AsynchTNTE_Bal.==0) ; count(AsynchTNOE_Bal.>0) count(AsynchTNOE_Bal.==0)]))

    AddSigIndicator(Bal_ChiSqTests[1,i],6,i,1,[1.12,.7])
    AddSigIndicator(Bal_ChiSqTests[2,i],6,i,1,[1.62,.7])
end
cd(string(topDir,"/Images"))
save(string("f4_prevalence_of_Pleiotropy_Bal_Fit.pdf"), f)
#endregion

#region NetworkSize plots
fontsize_theme = Theme(fontsize = 28)
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (nothing,nothing,4,25)
ga = f[1, 1] = GridLayout()
xlabels = ["Independent Effector", "\n Synchronous", "Shared Effector","Independent Effector", "\n Asynchronous", "Shared Effector"]
ga[1, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Scarce Resources \n Network Size",limits = lims,xgridvisible = false,ygridvisible = false)
ga[2, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Alternating Resources \n Network Size",limits = lims,xgridvisible = false,ygridvisible = false)
ga[3, 1] = Axis(f,xticks = ([1,1.125,1.25,1.5,1.625,1.75], xlabels),ylabel = "Plentiful Resources \n Network Size",limits = lims,xgridvisible = false,ygridvisible = false)
KruskalWallis = zeros(3)
Bal_SignRank = zeros(3,6)
for i in 1:3
    SynchTNTE_Bal = sum.(SynchTNTE_Bal_NetSizes[i]).-1
    SynchTNOE_Bal = sum.(SynchTNOE_Bal_NetSizes[i]).+1

    AsynchTNTE_Bal = sum.(AsynchTNTE_Bal_NetSizes[i]).-1#subtract 1 to account for the extra protein these networks start with
    AsynchTNOE_Bal = sum.(AsynchTNOE_Bal_NetSizes[i]).+1#add 1 to account for the shared effector which is otherwise not included in netsize
    
    KruskalWallis[i] = pvalue(KruskalWallisTest(SynchTNTE_Bal,SynchTNOE_Bal,AsynchTNTE_Bal,AsynchTNOE_Bal))
    if KruskalWallis[i] < .05
        group = [SynchTNTE_Bal,SynchTNOE_Bal,AsynchTNTE_Bal,AsynchTNOE_Bal]
        tests = [[1,4],[1,3],[1,2],[2,4],[2,3],[3,4]]
        for j in 1:6
            Bal_SignRank[i,j] = pvalue(SignedRankTest(group[tests[j][1]],group[tests[j][2]]))
        end
    end

    boxplot!(ga[i,1],vec(ones(length(SynchTNTE_Bal))), Float64.(SynchTNTE_Bal), color = :red, width = .25)
    boxplot!(ga[i,1],vec(ones(length(SynchTNOE_Bal))).+.25, Float64.(SynchTNOE_Bal), color = :black, width = .25)
    AddMedianLines([1.15,1.35],SynchTNOE_Bal,length(SynchTNOE_Bal),i,1,:red)

    boxplot!(ga[i,1],vec(ones(length(AsynchTNTE_Bal))).+.5, Float64.(AsynchTNTE_Bal), color = :red, width = .25)
    boxplot!(ga[i,1],vec(ones(length(AsynchTNOE_Bal))).+.75, Float64.(AsynchTNOE_Bal), color = :black, width = .25)
    AddMedianLines([1.65,1.85],AsynchTNOE_Bal,length(AsynchTNOE_Bal),i,1,:red)    
end 
text!(ga[1,1],string("A"),position = (1.05,21.75))
text!(ga[1,1],string("A"),position = (1.3,21.75))
text!(ga[1,1],string("B"),position = (1.55,21.75))
text!(ga[1,1],string("C"),position = (1.8,21.75))

text!(ga[2,1],string("A"),position = (1.05,21.75))
text!(ga[2,1],string("B"),position = (1.3,21.75))
text!(ga[2,1],string("C"),position = (1.55,21.75))
text!(ga[2,1],string("D"),position = (1.8,21.75))

text!(ga[3,1],string("A"),position = (1.05,21.75))
text!(ga[3,1],string("B"),position = (1.3,21.75))
text!(ga[3,1],string("C"),position = (1.55,21.75))
text!(ga[3,1],string("D"),position = (1.8,21.75))
set_theme!(fontsize_theme)
cd(string(topDir,"/Images"))
save(string("S7_Network_Sizes_Bal_Fit.pdf"), f, pt_per_unit = 1)
#endregion
#region Connectivity plots
fontsize_theme = Theme(fontsize = 28)
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (nothing,nothing,0,1)
ga = f[1, 1] = GridLayout()
xlabels = ["Independent Effector", "\n Synchronous", "Shared Effector","Independent Effector", "\n Asynchronous", "Shared Effector"]
ga[1, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Scarce Resources \n Network Connectivity",limits = lims,xgridvisible = false,ygridvisible = false)
ga[2, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Alternating Resources \n Network Connectivity",limits = lims,xgridvisible = false,ygridvisible = false)
ga[3, 1] = Axis(f,xticks = ([1,1.125,1.25,1.5,1.625,1.75], xlabels),ylabel = "Plentiful Resources \n Network Connectivity",limits = lims,xgridvisible = false,ygridvisible = false)
KruskalWallis = zeros(3)
Bal_SignRank = zeros(3,6)
for i in 1:3
    SynchTNTE_Bal = TE_Connectivity(SynchTNTE_Bal_NetSizes[i],SynchTNTE_Bal_Nets[i])
    SynchTNOE_Bal = OE_Connectivity(SynchTNOE_Bal_NetSizes[i],SynchTNOE_Bal_Nets[i])

    AsynchTNTE_Bal = TE_Connectivity(AsynchTNTE_Bal_NetSizes[i],AsynchTNTE_Bal_Nets[i])
    AsynchTNOE_Bal = OE_Connectivity(AsynchTNOE_Bal_NetSizes[i],AsynchTNOE_Bal_Nets[i])
    
    KruskalWallis[i] = pvalue(KruskalWallisTest(SynchTNTE_Bal,SynchTNOE_Bal,AsynchTNTE_Bal,AsynchTNOE_Bal))
    if KruskalWallis[i] < .05
        group = [SynchTNTE_Bal,SynchTNOE_Bal,AsynchTNTE_Bal,AsynchTNOE_Bal]
        tests = [[1,4],[1,3],[1,2],[2,4],[2,3],[3,4]]
        for j in 1:6
            Bal_SignRank[i,j] = pvalue(SignedRankTest(group[tests[j][1]],group[tests[j][2]]))
        end
    end

    boxplot!(ga[i,1],vec(ones(length(SynchTNTE_Bal))), Float64.(SynchTNTE_Bal), color = :red, width = .25)
    boxplot!(ga[i,1],vec(ones(length(SynchTNOE_Bal))).+.25, Float64.(SynchTNOE_Bal), color = :black, width = .25)
    AddMedianLines([1.15,1.35],SynchTNOE_Bal,length(SynchTNOE_Bal),i,1,:red)

    boxplot!(ga[i,1],vec(ones(length(AsynchTNTE_Bal))).+.5, Float64.(AsynchTNTE_Bal), color = :red, width = .25)
    boxplot!(ga[i,1],vec(ones(length(AsynchTNOE_Bal))).+.75, Float64.(AsynchTNOE_Bal), color = :black, width = .25)
    AddMedianLines([1.65,1.85],AsynchTNOE_Bal,length(AsynchTNOE_Bal),i,1,:red)    
end 
text!(ga[1,1],string("A"),position = (1.05,.75))
text!(ga[1,1],string("B"),position = (1.3,.75))
text!(ga[1,1],string("C"),position = (1.55,.75))
text!(ga[1,1],string("A"),position = (1.8,.75))

text!(ga[2,1],string("A"),position = (1.05,.75))
text!(ga[2,1],string("B"),position = (1.3,.75))
text!(ga[2,1],string("C"),position = (1.55,.75))
text!(ga[2,1],string("D"),position = (1.8,.75))

text!(ga[3,1],string("A"),position = (1.05,.75))
text!(ga[3,1],string("B"),position = (1.3,.75))
text!(ga[3,1],string("A"),position = (1.55,.75))
text!(ga[3,1],string("C"),position = (1.8,.75))
set_theme!(fontsize_theme)
cd(string(topDir,"/Images"))
save(string("S8_Network_Connectivity_Bal_Fit_alt.pdf"), f)
#endregion

#region effector connectivity plot
fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (nothing,nothing,0,15)
ga = f[1, 1] = GridLayout()
xlabels = ["Independent Effector", "\n Synchronous", "Shared Effector","Independent Effector", "\n Asynchronous", "Shared Effector"]
ga[1, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Scarce Resources \n Effector Connections",limits = lims,xgridvisible = false,ygridvisible = false)
ga[2, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Alternating Resources \n Effector Connections",limits = lims,xgridvisible = false,ygridvisible = false)
ga[3, 1] = Axis(f,xticks = ([1,1.125,1.25,1.5,1.625,1.75], xlabels),ylabel = "Plentiful Resources \n Effector Connections",limits = lims,xgridvisible = false,ygridvisible = false)
KruskalWallis = zeros(3)
Bal_SignRank = zeros(3,6)
for i in 1:3
    SynchTNTE_Bal = TE_EffectorCons(SynchTNTE_Bal_NetSizes[i],SynchTNTE_Bal_Nets[i])
    SynchTNOE_Bal = OE_EffectorCons(SynchTNOE_Bal_NetSizes[i],SynchTNOE_Bal_Nets[i])

    AsynchTNTE_Bal = TE_EffectorCons(AsynchTNTE_Bal_NetSizes[i],AsynchTNTE_Bal_Nets[i])
    AsynchTNOE_Bal = OE_EffectorCons(AsynchTNOE_Bal_NetSizes[i],AsynchTNOE_Bal_Nets[i])
    
    KruskalWallis[i] = pvalue(KruskalWallisTest(SynchTNTE_Bal,SynchTNOE_Bal,AsynchTNTE_Bal,AsynchTNOE_Bal))
    if KruskalWallis[i] < .05
        group = [SynchTNTE_Bal,SynchTNOE_Bal,AsynchTNTE_Bal,AsynchTNOE_Bal]
        tests = [[1,4],[1,3],[1,2],[2,4],[2,3],[3,4]]
        for j in 1:6
            Bal_SignRank[i,j] = pvalue(SignedRankTest(group[tests[j][1]],group[tests[j][2]]))
        end
    end

    boxplot!(ga[i,1],vec(ones(length(SynchTNTE_Bal))), Float64.(SynchTNTE_Bal), color = :red, width = .25)
    boxplot!(ga[i,1],vec(ones(length(SynchTNOE_Bal))).+.25, Float64.(SynchTNOE_Bal), color = :black, width = .25)
    AddMedianLines([1.15,1.35],SynchTNOE_Bal,length(SynchTNOE_Bal),i,1,:red)

    boxplot!(ga[i,1],vec(ones(length(AsynchTNTE_Bal))).+.5, Float64.(AsynchTNTE_Bal), color = :red, width = .25)
    boxplot!(ga[i,1],vec(ones(length(AsynchTNOE_Bal))).+.75, Float64.(AsynchTNOE_Bal), color = :black, width = .25)
    AddMedianLines([1.65,1.85],AsynchTNOE_Bal,length(AsynchTNOE_Bal),i,1,:red)  
end 
text!(ga[1,1],string("A"),position = (1.05,12.75))
text!(ga[1,1],string("B"),position = (1.3,12.75))
text!(ga[1,1],string("B"),position = (1.55,12.75))
text!(ga[1,1],string("C"),position = (1.8,12.75))

text!(ga[2,1],string("A"),position = (1.05,12.75))
text!(ga[2,1],string("B"),position = (1.3,12.75))
text!(ga[2,1],string("C"),position = (1.55,12.75))
text!(ga[2,1],string("D"),position = (1.8,12.75))

text!(ga[3,1],string("A"),position = (1.05,12.75))
text!(ga[3,1],string("B"),position = (1.3,12.75))
text!(ga[3,1],string("C"),position = (1.55,12.75))
text!(ga[3,1],string("D"),position = (1.8,12.75))
set_theme!(fontsize_theme)
cd(string(topDir,"/Images"))
save(string("S9_Effector_Connectivity.pdf"), f)
#endregion
