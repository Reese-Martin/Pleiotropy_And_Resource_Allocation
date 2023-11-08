# this script creates figure 1, the comparative fitness plot, and s11 the plot showing early fitness traces
using Statistics
using StatsBase
using FileIO
using HypothesisTests
using CairoMakie #use CairoMakie for final figures because it can save as EPS
#using GLMakie #use GLMakie for 'prototyping' ie not the final graphs because of the pop out window and iteration
using GeometryBasics
#functions
function Jitter(a,Factor = .3)
    a .+= rand(-Factor:.001:Factor,length(a))
end
function AddMedianLines(x,Data,Runs, row, col, Color)
    Rows = Runs
    tmp = median(reshape(Data,(Rows,1)),dims = 1)
    lines!(ga[row,col],x,[tmp[1],tmp[1]], linewidth = 2, color = Color)
end 
function AddGeoMLines(x,Data,Runs, row, col, Color)
    Rows = Runs
    tmp = geomean(reshape(Data,(Rows,1)))
    lines!(ga[row,col],x,[tmp[1],tmp[1]], linewidth = 2, color = Color)
end 
function AddSigIndicator(x,Tests,row,col,pos)
    for j in 1:length(x)
        if x[j] < .05/Tests
            text!(ga[row,col],pos[1],pos[2],text = "*",fontsize = 60)
        end
    end
end
#region load data
# change to directory that holds both a data folder and an Images folder to read data and save output
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Two_Network_Model/"
Data = string(topDir,"Data/")
cd(Data)
# so that only the necessary is included
Files = filter(x->endswith(x, ".jld2"), readdir())
Files = filter(x->!contains(x, "reverse"), Files)

TNTEFiles = filter(x -> contains(x, "TNTE"), Files)
AsynchTNTEFiles = filter(x -> contains(x, "Asynch"), TNTEFiles)
AsynchTNTEFiles = filter(x -> !contains(x, "Net"), AsynchTNTEFiles)
AsynchTNTEFiles = filter(x -> !contains(x, "vs"), AsynchTNTEFiles)
AsynchTNTEFiles = filter(x -> !contains(x, "ParFit"), AsynchTNTEFiles)
AsynchTNTEBalFiles = filter(x -> contains(x, "Balanced"), AsynchTNTEFiles)[[2,3,1]] #grab the variable resource condition 50 for comparison

SynchTNTEFiles = filter(x -> contains(x, "Synch"), TNTEFiles)
SynchTNTEFiles = filter(x -> !contains(x, "Net"), SynchTNTEFiles)
SynchTNTEFiles = filter(x -> !contains(x, "vs"), SynchTNTEFiles)
SynchTNTEFiles = filter(x -> !contains(x, "ParFit"), SynchTNTEFiles)
SynchTNTEBalFiles = filter(x -> contains(x, "Balanced"), SynchTNTEFiles)[[2,3,1]] #grab the variable resource condition 50 for comparison

TNOEFiles = filter(x -> contains(x, "TNOE"), Files)
AsynchTNOEFiles = filter(x -> contains(x, "Asynch"), TNOEFiles)
AsynchTNOEFiles = filter(x -> !contains(x, "Net"), AsynchTNOEFiles)
AsynchTNOEFiles = filter(x -> !contains(x, "vs"), AsynchTNOEFiles)
AsynchTNOEFiles = filter(x -> !contains(x, "ParFit"), AsynchTNOEFiles)
AsynchTNOEBalFiles = filter(x -> contains(x, "Balanced"), AsynchTNOEFiles)[[2,3,1]] #grab the variable resource condition 50 for comparison

SynchTNOEFiles = filter(x -> contains(x, "Synch"), TNOEFiles)
SynchTNOEFiles = filter(x -> !contains(x, "Net"), SynchTNOEFiles)
SynchTNOEFiles = filter(x -> !contains(x, "vs"), SynchTNOEFiles)
SynchTNOEFiles = filter(x -> !contains(x, "ParFit"), SynchTNOEFiles)
SynchTNOEBalFiles = filter(x -> contains(x, "Balanced"), SynchTNOEFiles)[[2,3,1]] #grab the variable resource condition 50 for comparison

AsynchTNTE_Bal_Fit = []
AsynchTNOE_Bal_Fit = []

SynchTNTE_Bal_Fit = []
SynchTNOE_Bal_Fit = []

for i in 1:3
    AsynchTNOE_Bal_Tmp = FileIO.load(string(Data,AsynchTNOEBalFiles[i]))
    AsynchTNTE_Bal_Tmp = FileIO.load(string(Data,AsynchTNTEBalFiles[i]))
    
    push!(AsynchTNTE_Bal_Fit, AsynchTNTE_Bal_Tmp["HostFit"])
    push!(AsynchTNOE_Bal_Fit,AsynchTNOE_Bal_Tmp["HostFit"])

    SynchTNTE_Bal_Tmp = FileIO.load(string(Data,SynchTNTEBalFiles[i]))
    SynchTNOE_Bal_Tmp = FileIO.load(string(Data,SynchTNOEBalFiles[i]))
    
    push!(SynchTNTE_Bal_Fit, SynchTNTE_Bal_Tmp["HostFit"])
    push!(SynchTNOE_Bal_Fit, SynchTNOE_Bal_Tmp["HostFit"])
end
#endregion

#region Figure 1 Host fitness in the final generation of independent evolution 4 groups per panel, 3 panels
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (nothing,nothing,0,.85)
ga = f[1, 1] = GridLayout()
xlabels = ["Independent Effector", "\n Synchronous", "Shared Effector","Independent Effector", "\n Asynchronous", "Shared Effector"]
ga[1, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Scarce Resources \n Absolute Fitness",limits = lims,xgridvisible = false,ygridvisible = false)
ga[2, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Alternating Resources  \n Absolute Fitness",limits = lims,xgridvisible = false,ygridvisible = false)
ga[3, 1] = Axis(f,xticks = ([1,1.25,1.5,2,2.25,2.5], xlabels),ylabel = "Plentiful Resources  \n Absolute Fitness",limits = lims,xgridvisible = false,ygridvisible = false)
p_big = decompose(Point2f, Circle(Point2f(0), .4))
p_small = decompose(Point2f, Circle(Point2f(0), 0.3))
mark = Polygon(p_big, [p_small])
KruskalWallis = zeros(3)
Bal_SignRank = zeros(3,6)
Bal_Synch_Variance = zeros(2,3)
Bal_Asynch_Variance = zeros(2,3)
Bal_FTest = zeros(2,3)
for i in 1:3
    SynchTNOE_Bal = SynchTNOE_Bal_Fit[i][:,end]
    SynchTNTE_Bal = SynchTNTE_Bal_Fit[i][:,end]

    AsynchTNOE_Bal = AsynchTNOE_Bal_Fit[i][:,end]
    AsynchTNTE_Bal = AsynchTNTE_Bal_Fit[i][:,end]    
   
    Bal_Synch_Variance[1,i] = var(SynchTNTE_Bal)
    Bal_Synch_Variance[2,i] = var(SynchTNOE_Bal)
    Bal_Asynch_Variance[1,i] = var(AsynchTNTE_Bal)
    Bal_Asynch_Variance[2,i] = var(AsynchTNOE_Bal)
    Bal_FTest[1,i] = pvalue(VarianceFTest(SynchTNTE_Bal,SynchTNOE_Bal)) 
    Bal_FTest[2,i] = pvalue(VarianceFTest(AsynchTNTE_Bal,AsynchTNOE_Bal)) 

    KruskalWallis[i] = pvalue(KruskalWallisTest(SynchTNTE_Bal,SynchTNOE_Bal,AsynchTNTE_Bal,AsynchTNOE_Bal))

    scatter!(ga[i,1],vec(Jitter(ones(length(SynchTNTE_Bal)), .2)), Float64.(SynchTNTE_Bal), color = :red)
    AddMedianLines([.8,1.2],SynchTNTE_Bal,length(SynchTNTE_Bal),i,1,:black)
    scatter!(ga[i,1],vec(Jitter(ones(length(SynchTNOE_Bal)), .2)).+.5, Float64.(SynchTNOE_Bal), color = :black)
    AddMedianLines([1.3,1.7],SynchTNOE_Bal,length(SynchTNOE_Bal),i,1,:red)

    scatter!(ga[i,1],vec(Jitter(ones(length(AsynchTNTE_Bal)), .2)).+1, marker = mark, Float64.(AsynchTNTE_Bal), color = :red)
    AddMedianLines([1.8,2.2],AsynchTNTE_Bal,length(AsynchTNTE_Bal),i,1,:black)
    scatter!(ga[i,1],vec(Jitter(ones(length(AsynchTNOE_Bal)), .2)).+1.5, marker = mark, Float64.(AsynchTNOE_Bal), color = :black)
    AddMedianLines([2.3,2.7],AsynchTNOE_Bal,length(AsynchTNOE_Bal),i,1,:red)
    if KruskalWallis[i] < .05
        group = [SynchTNTE_Bal,SynchTNOE_Bal,AsynchTNTE_Bal,AsynchTNOE_Bal]
        tests = [[1,4],[1,3],[1,2],[2,4],[2,3],[3,4]]
        for j in 1:6
            Bal_SignRank[i,j] = pvalue(SignedRankTest(group[tests[j][1]],group[tests[j][2]]))
        end
    end
end 
text!(ga[1,1],string("A"),position = (1,.75))
text!(ga[1,1],string("B"),position = (1.5,.75))
text!(ga[1,1],string("A"),position = (2,.75))
text!(ga[1,1],string("C"),position = (2.5,.75))
text!(ga[2,1],string("A"),position = (1,.75))
text!(ga[2,1],string("B"),position = (1.5,.75))
text!(ga[2,1],string("A"),position = (2,.75))
text!(ga[2,1],string("C"),position = (2.5,.75))
text!(ga[3,1],string("A"),position = (1,.75))
text!(ga[3,1],string("B"),position = (1.5,.75))
text!(ga[3,1],string("C"),position = (2,.75))
text!(ga[3,1],string("D"),position = (2.5,.75))
fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)
cd(string(topDir,"/Images"))
CairoMakie.save("F1_Final_Gen_HostFit.pdf", f)
#endregion

#region Figure S11 Host fitness across 100 generations
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1500, 1000))
lims = (nothing,nothing,0,1)
ga = f[1, 1] = GridLayout()
ga[1, 1] = Axis(f,xticklabelsvisible = false,xticksvisible = false,ylabel = "Synchronous Signaling \n Absolute Fitness", title = "Scarce Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[1, 2] = Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = "Alternating Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[1, 3] = Axis(f,xticklabelsvisible = false,xticksvisible = false,yticklabelsvisible = false,yticksvisible = false,title = "Plentiful Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[2, 1] = Axis(f,ylabel = "Asynchronous Signaling \n Absolute Fitness",limits = lims,xgridvisible = false,ygridvisible = false)
ga[2, 2] = Axis(f,xlabel = "Generations",yticklabelsvisible = false,yticksvisible = false,limits = lims,xgridvisible = false,ygridvisible = false)
ga[2, 3] = Axis(f,yticklabelsvisible = false,yticksvisible = false,limits = lims,xgridvisible = false,ygridvisible = false)

for i in 1:3
    SynchTNTE_Bal = median(SynchTNTE_Bal_Fit[i][:,1:100],dims = 1)
    SynchTNOE_Bal = median(SynchTNOE_Bal_Fit[i][:,1:100],dims = 1)

    AsynchTNTE_Bal = median(AsynchTNTE_Bal_Fit[i][:,1:100],dims = 1)
    AsynchTNOE_Bal = median(AsynchTNOE_Bal_Fit[i][:,1:100],dims = 1)

    a = lines!(ga[1,i],vec(SynchTNTE_Bal), linewidth = 3, color = :red)
    b = lines!(ga[1,i],vec(SynchTNOE_Bal), linewidth = 3, color = :black)

    lines!(ga[2,i],vec(AsynchTNTE_Bal), linewidth = 3, color = :red)
    lines!(ga[2,i],vec(AsynchTNOE_Bal), linewidth = 3, color = :black)
    labels = ["Independent Effector", "Shared Effector"]
    title = "Host Type"
    # elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    group_size = [a,b]
    Legend(f[1,2], group_size, labels, title)
end 
fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)
cd(string(topDir,"/Images"))
CairoMakie.save("S11_HostFit_100.pdf", f)
#endregion
