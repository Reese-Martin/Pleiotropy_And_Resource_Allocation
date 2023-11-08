# this file is the final file for figure 1 in the Two net paper, will trim excess code 
# Based on Comparative EC fitness v2
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
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Two_Network_Model/"
Data = string(topDir,"Data/")
cd(Data)
# so that only the necessary is included
Files = filter(x->contains(x, "reverse"), readdir())
TNTEFiles_reverse = filter(x -> contains(x, "TNTE"), Files)
AsynchTNTEFiles_reverse = filter(x -> contains(x, "Asynch"), TNTEFiles_reverse)
AsynchTNTEFiles_reverse = filter(x -> !contains(x, "Net"), AsynchTNTEFiles_reverse)
AsynchTNTEFiles_reverse = filter(x -> !contains(x, "vs"), AsynchTNTEFiles_reverse)
AsynchTNTEFiles_reverse = filter(x -> !contains(x, "ParFit"), AsynchTNTEFiles_reverse)
AsynchTNTEBalFiles_reverse = filter(x -> contains(x, "Balanced"), AsynchTNTEFiles_reverse)[[2,1,3]] #grab the variable resource condition 50 for comparison

SynchTNTEFiles_reverse = filter(x -> contains(x, "Synch"), TNTEFiles_reverse)
SynchTNTEFiles_reverse = filter(x -> !contains(x, "Net"), SynchTNTEFiles_reverse)
SynchTNTEFiles_reverse = filter(x -> !contains(x, "vs"), SynchTNTEFiles_reverse)
SynchTNTEFiles_reverse = filter(x -> !contains(x, "ParFit"), SynchTNTEFiles_reverse)
SynchTNTEBalFiles_reverse = filter(x -> contains(x, "Balanced"), SynchTNTEFiles_reverse)[[2,1,3]] #grab the variable resource condition 50 for comparison

TNOEFiles_reverse = filter(x -> contains(x, "TNOE"), Files)
AsynchTNOEFiles_reverse = filter(x -> contains(x, "Asynch"), TNOEFiles_reverse)
AsynchTNOEFiles_reverse = filter(x -> !contains(x, "Net"), AsynchTNOEFiles_reverse)
AsynchTNOEFiles_reverse = filter(x -> !contains(x, "vs"), AsynchTNOEFiles_reverse)
AsynchTNOEFiles_reverse = filter(x -> !contains(x, "ParFit"), AsynchTNOEFiles_reverse)
AsynchTNOEBalFiles_reverse = filter(x -> contains(x, "Balanced"), AsynchTNOEFiles_reverse)[[2,1,3]] #grab the variable resource condition 50 for comparison

SynchTNOEFiles_reverse = filter(x -> contains(x, "Synch"), TNOEFiles_reverse)
SynchTNOEFiles_reverse = filter(x -> !contains(x, "Net"), SynchTNOEFiles_reverse)
SynchTNOEFiles_reverse = filter(x -> !contains(x, "vs"), SynchTNOEFiles_reverse)
SynchTNOEFiles_reverse = filter(x -> !contains(x, "ParFit"), SynchTNOEFiles_reverse)
SynchTNOEBalFiles_reverse = filter(x -> contains(x, "Balanced"), SynchTNOEFiles_reverse)[[2,1,3]] #grab the variable resource condition 50 for comparison

AsynchTNTE_Bal_Fit_reverse = []
AsynchTNOE_Bal_Fit_reverse = []

SynchTNTE_Bal_Fit_reverse = []
SynchTNOE_Bal_Fit_reverse = []

for i in 1:3
    AsynchTNOE_Bal_Tmp = FileIO.load(string(Data,AsynchTNOEBalFiles_reverse[i]))
    AsynchTNTE_Bal_Tmp = FileIO.load(string(Data,AsynchTNTEBalFiles_reverse[i]))
    
    push!(AsynchTNTE_Bal_Fit_reverse, AsynchTNTE_Bal_Tmp["HostFit"])
    push!(AsynchTNOE_Bal_Fit_reverse,AsynchTNOE_Bal_Tmp["HostFit"])

    SynchTNTE_Bal_Tmp = FileIO.load(string(Data,SynchTNTEBalFiles_reverse[i]))
    SynchTNOE_Bal_Tmp = FileIO.load(string(Data,SynchTNOEBalFiles_reverse[i]))
    
    push!(SynchTNTE_Bal_Fit_reverse, SynchTNTE_Bal_Tmp["HostFit"])
    push!(SynchTNOE_Bal_Fit_reverse, SynchTNOE_Bal_Tmp["HostFit"])
end
#endregion

#region reverse absolute fitness calculations
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
    SynchTNOE_Bal = SynchTNOE_Bal_Fit_reverse[i][:,end]
    SynchTNTE_Bal = SynchTNTE_Bal_Fit_reverse[i][:,end]

    AsynchTNOE_Bal = AsynchTNOE_Bal_Fit_reverse[i][:,end]
    AsynchTNTE_Bal = AsynchTNTE_Bal_Fit_reverse[i][:,end]    
   
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
text!(ga[1,1],string("C"),position = (2,.75))
text!(ga[1,1],string("B"),position = (2.5,.75))
text!(ga[2,1],string("A"),position = (1,.75))
text!(ga[2,1],string("A"),position = (1.5,.75))
text!(ga[2,1],string("A"),position = (2,.75))
text!(ga[2,1],string("B"),position = (2.5,.75))
text!(ga[3,1],string("A"),position = (1,.75))
text!(ga[3,1],string("B"),position = (1.5,.75))
text!(ga[3,1],string("C"),position = (2,.75))
text!(ga[3,1],string("D"),position = (2.5,.75))
fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)
cd(string(topDir,"/Images"))
# save(string("Example_Plots.png"), f, pt_per_unit = 1)
# save(string("Final_Gen_HostFit_BalFitContributions.png"), f, pt_per_unit = 1)
# CairoMakie.save("F1_Final_Gen_HostFit.pdf", f)
CairoMakie.save("S2_absolute_fitness_reverse.pdf", f)
#endregion