# Figure 2 is the quasi binomial regression and proportion of population plots
using StatsBase
using TypedTables
using Statistics
using FileIO
using HypothesisTests
using CairoMakie #use CairoMakie for final figures because it can save as EPS
# using GLMakie #use GLMakie for 'prototyping' ie not the final graphs because of the pop out window and iteration
using GLM
using QuasiGLM
using DataFrames
using GeometryBasics
#region functions
function Extend_Final_Per!(x) # fix the competitive simulation matrices so that any fixation events show as bein fixed for the entirety of the competition (even playing field for comparison)
    for i in 1:size(x,1)
        if maximum(x[i,:]) == 1
            x[i,findfirst(x -> x == 1,x[i,:]):end] .=1
        end
    end
end 
function Jitter(a,Factor = .3)
    a .+= rand(-Factor:.001:Factor,length(a))
end
function AddMedianLines(x,Data,Runs, row,col)
    Rows = Runs
    tmp = median(reshape(Data,(Rows,1)),dims = 1)
    lines!(ga[row,col],x,[tmp[1],tmp[1]], linewidth = 5, color = :Black)
end
function AddGeoMeanLines(x,Data,Runs, row,col)
    Rows = Runs
    tmp = geomean(reshape(Data,(Rows,1)),dims = 1)
    lines!(ga[row,col],x,[tmp[1],tmp[1]], linewidth = 5, color = :Black)
end  
function NanMean(x,dim)
    tmp = zeros(size(x,dim))
    for i in 1:size(x,dim)
        if dim == 1
            tmp[i] = mean(x[i,(!).(isnan.(x[i,:]))])
        elseif dim == 1
            tmp[i] = mean(x[(!).(isnan.(x[:,i])),i])
        end
    end
    return tmp
end
function AddSigIndicator(x,Tests,row,col,pos)
    for j in 1:length(x)
        if x[j] < .05/Tests
            text!(ga[row,col],pos[1],pos[2],text = "*",fontsize = 60)
        end
    end
end
#endregion

#region load data 
topDir = "/Users/martinra4/Lab Work/OneDrive - Vanderbilt/Vanderbilt/Tate Lab/Two_Network_Model/"
Data = string(topDir,"Data/")
cd(Data)
Files = filter(x->endswith(x, ".jld2"), readdir())
Files = filter(x -> contains(x, "vs"), Files)

AsynchFiles = filter(x -> contains(x, "Asynch"), Files)
AsynchFiles = filter(x -> !contains(x, "Networks"), AsynchFiles)

SynchFiles = filter(x -> contains(x, "Synch"), Files)
SynchFiles = filter(x -> !contains(x, "Networks"), SynchFiles)

Synch_OE_Fit_Files = filter(x -> contains(x, "HostFitOE"), SynchFiles)[[5,6,4,8,9,7,2,3,1]]#[[2,1,3]]#[[5,6,4,2,3,1,8,9,7]]
Synch_TE_Fit_Files = filter(x -> contains(x, "HostFitTE"), SynchFiles)[[5,6,4,8,9,7,2,3,1]]#[[2,1,3]]#[[5,6,4,2,3,1,8,9,7]]
Synch_Per_Pop_TE_Files = filter(x -> contains(x, "PerPopTE"), SynchFiles)[[5,6,4,8,9,7,2,3,1]]#[[2,1,3]]#[[5,6,4,2,3,1,8,9,7]]

Asynch_OE_Fit_Files = filter(x -> contains(x, "HostFitOE"), AsynchFiles)[[5,6,4,8,9,7,2,3,1]]#[[2,1,3]]#[[5,6,4,2,3,1,8,9,7]]
Asynch_TE_Fit_Files = filter(x -> contains(x, "HostFitTE"), AsynchFiles)[[5,6,4,8,9,7,2,3,1]]#[[2,1,3]]#[[5,6,4,2,3,1,8,9,7]]
Asynch_Per_Pop_TE_Files = filter(x -> contains(x, "PerPopTE"), AsynchFiles)[[5,6,4,8,9,7,2,3,1]]#[[2,1,3]]#[[5,6,4,2,3,1,8,9,7]]

SynchBal_Per_TE = []
SynchBal_Host_Fit_TE = [] 
SynchBal_Host_Fit_OE = []

AsynchBal_Per_TE= []
AsynchBal_Host_Fit_TE = [] 
AsynchBal_Host_Fit_OE = []


for i in 1:9
    Synch_OE_Tmp = FileIO.load(string(Data,Synch_OE_Fit_Files[i]))
    Synch_TE_Tmp = FileIO.load(string(Data,Synch_TE_Fit_Files[i]))
    Synch_PPT_Tmp = FileIO.load(string(Data,Synch_Per_Pop_TE_Files[i]))
    push!(SynchBal_Per_TE, Synch_PPT_Tmp["PerPopTE"])
    push!(SynchBal_Host_Fit_TE, Synch_TE_Tmp["HostFitTE"])
    push!(SynchBal_Host_Fit_OE, Synch_OE_Tmp["HostFitOE"])

    Asynch_OE_Tmp = FileIO.load(string(Data,Asynch_OE_Fit_Files[i]))
    Asynch_TE_Tmp = FileIO.load(string(Data,Asynch_TE_Fit_Files[i]))
    Asynch_PPT_Tmp = FileIO.load(string(Data,Asynch_Per_Pop_TE_Files[i]))
    push!(AsynchBal_Per_TE, Asynch_PPT_Tmp["PerPopTE"])
    push!(AsynchBal_Host_Fit_TE, Asynch_TE_Tmp["HostFitTE"])
    push!(AsynchBal_Host_Fit_OE, Asynch_OE_Tmp["HostFitOE"])
end

#in event of a fixation, extend the final percnt of the population through the remaining generations
Extend_Final_Per!.(SynchBal_Per_TE)
Extend_Final_Per!.(AsynchBal_Per_TE)
# discard simulations in which there was a draw
Synch_Final = []
Asynch_Final = []
for i in 1:9
    push!(Synch_Final,1 .- SynchBal_Per_TE[i][:,end])
    push!(Asynch_Final,1 .- AsynchBal_Per_TE[i][:,end])
end
#endregion

#region generate dataframe for regression analysis
Synch_Final = reduce(vcat,Synch_Final)
Asynch_Final = reduce(vcat,Asynch_Final)
Synch = fill(1,length(Synch_Final)) #binary, 1 = synch, 0 = asynch
Asynch = fill(0,length(Synch_Final))
EQGens = reduce(vcat,repeat([fill(250,1000),fill(500,1000),fill(1000,1000)],3)) #number of burn in generations
Resources = reduce(vcat,[fill(.1,3000),fill(.5,3000),fill(1,3000)]) #scarce, alternating, plentiful

#DF containing resutls and conditions that lead to those results
DF = DataFrame(Outcomes = vcat(Synch_Final,Asynch_Final), SigTime = vcat(Synch,Asynch), EQgens = vcat(EQGens,EQGens), Resources = vcat(Resources,Resources))
#fit model
BinModel =  glm(@formula(Outcomes ~ SigTime + EQgens + Resources), DF, Binomial(), LogitLink())
AdjustQuasiGLM(BinModel, DF; level=0.95) #run function to adjust model to quasi-Binomial data

#generate test table 
synch_timing = fill(1,9)
asynch_timing = fill(0,9)
gens = repeat([250,500,1000],3)
resources = reduce(vcat,[fill(.1,3),fill(.5,3),fill(1,3)])
Synch_Table = Table(SigTime = synch_timing, EQgens = gens, Resources = resources)
Asycnh_Table = Table(SigTime = asynch_timing, EQgens = gens, Resources = resources)

synch_Preds = predict(BinModel,Synch_Table) # 1- to account for the data being in a Indep winner format rather than a Shared winner format
asynch_Preds = predict(BinModel,Asycnh_Table) # 1- to account for the data being in a Indep winner format rather than a Shared winner format
#endregion

#region line plots of competitions won by One Effector hosts
# Plot comparing winning percentages for the synch and asynch conditions when fitness is balanced
fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (1950, 1000))
lims = (.7,3.3,-.05,1.05)
ga = f[1, 1] = GridLayout()
ga[1, 1] = Axis(f,xticks = ([1,2,3],["250", "500", "1000"]),ylabel = "Proportion of population \n with shared effectors", title = "Scarce Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[1, 2] = Axis(f,xticks = ([1,2,3],["250", "500", "1000"]),xlabel = "Generations Prior to Competition", title = "Alternating Resources",yticklabelsvisible = false,yticksvisible = false,limits = lims,xgridvisible = false,ygridvisible = false)
ga[1, 3] = Axis(f,xticks = ([1,2,3],["250", "500", "1000"]), title = "Plentiful Resources",yticklabelsvisible = false,yticksvisible = false,limits = lims,xgridvisible = false,ygridvisible = false)
p_big = decompose(Point2f, Circle(Point2f(0), .4))
p_small = decompose(Point2f, Circle(Point2f(0), 0.3))
mark = Polygon(p_big, [p_small])
for i in 1:3 
    a = lines!(ga[1,i],[.8,1],[synch_Preds[3*(i-1)+1],synch_Preds[3*(i-1)+1]], color = :black, linewidth = 7)
    lines!(ga[1,i],[1.8,2],[synch_Preds[3*(i-1)+2],synch_Preds[3*(i-1)+2]], color = :black, linewidth = 7)
    lines!(ga[1,i],[2.8,3],[synch_Preds[3*(i-1)+3],synch_Preds[3*(i-1)+3]], color = :black, linewidth = 7)

    tmp = countmap(SynchBal_Per_TE[(3*(i-1)+1)][:,end])
    tmp2 = countmap(SynchBal_Per_TE[(3*(i-1)+2)][:,end])
    tmp3 = countmap(SynchBal_Per_TE[(3*(i-1)+3)][:,end])

    scatter!(ga[1,i],ones(length(tmp)) .- .1,1 .- collect(keys(tmp)), color = :black, marker = :circle, markersize = ceil.(collect(values(tmp))./20).+4,strokewidth = 1)
    scatter!(ga[1,i],ones(length(tmp2)).+.9,1 .- collect(keys(tmp2)), color = :black, marker = :circle, markersize = ceil.(collect(values(tmp2))./20).+4, strokewidth = 1)
    scatter!(ga[1,i],ones(length(tmp3)).+1.9,1 .- collect(keys(tmp3)), color = :black, marker = :circle, markersize = ceil.(collect(values(tmp3))./20).+4,strokewidth = 1)

    b = lines!(ga[1,i],[1,1.2],[asynch_Preds[3*(i-1)+1],asynch_Preds[3*(i-1)+1]], color = :grey, linewidth = 7)
    lines!(ga[1,i],[2,2.2],[asynch_Preds[3*(i-1)+2],asynch_Preds[3*(i-1)+2]], color = :grey, linewidth = 7)
    lines!(ga[1,i],[3,3.2],[asynch_Preds[3*(i-1)+3],asynch_Preds[3*(i-1)+3]], color = :grey, linewidth = 7)

    tmp4 = countmap(AsynchBal_Per_TE[(3*(i-1)+1)][:,end])
    tmp5 = countmap(AsynchBal_Per_TE[(3*(i-1)+2)][:,end])
    tmp6 = countmap(AsynchBal_Per_TE[(3*(i-1)+3)][:,end])

    scatter!(ga[1,i],ones(length(tmp4)) .+ .1,1 .- collect(keys(tmp4)), color = :black, marker = mark, markersize = ceil.(collect(values(tmp4))./20).+4,strokewidth = 1)
    scatter!(ga[1,i],ones(length(tmp5)).+1.1,1 .- collect(keys(tmp5)), color = :black, marker = mark, markersize = ceil.(collect(values(tmp5))./20).+4,strokewidth = 1)
    scatter!(ga[1,i],ones(length(tmp6)).+2.1,1 .- collect(keys(tmp6)), color = :black, marker = mark, markersize = ceil.(collect(values(tmp6))./20).+4,strokewidth = 1)

    labels = ["Synchronous", "Asynchronous","Synch. Win %", "Asynch. Win %"]
    title = "Signaling Timing"

    group_size = [MarkerElement(marker = :circle, color = :black),MarkerElement(marker = mark, color = :black),a,b]
    Legend(f[1,2], group_size, labels, title)
end
# set 50% mark using dashed lines
lines!(ga[1,1],[.7,3.3],[.5,.5],linestyle = :dash)
lines!(ga[1,2],[.7,3.3],[.5,.5],linestyle = :dash)
lines!(ga[1,3],[.7,3.3],[.5,.5],linestyle = :dash)

cd(string(topDir,"/Images"))
save(string("F2_Probability_of_OE_Winning.pdf"), f)
#endregion