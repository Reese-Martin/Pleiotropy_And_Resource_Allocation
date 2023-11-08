# v4 refactors plots to show asynch and synch in the same plots for better viewing and removes reference to the non-standard Effector use cases (save these for supplement)
# if looking for the other effector use cases check v3

using Statistics
using FileIO
using HypothesisTests
using CairoMakie #use CairoMakie for final figures because it can save as EPS
# using GLMakie #use GLMakie for 'prototyping' ie not the final graphs because of the pop out window and iteration

#comp net 1,2 indicates the networks being compared to the standard networks
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
#endregion

#region mean difference in Immune cost between TE and OE 
# immune cost is sheet 1 of the given Host fit array (immune cost for Bal_Host_Fit_TE[1] is [:,:,1])
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (2000, 1000))
lims = (.9,3.1,-.5,.15)
colors = [:Blue, :Orange]
ga = f[1, 1] = GridLayout()
ga[1, 1] = Axis(f,xticks = ([1,2,3],["250", "500", "1000"]),ylabel = "Shared Effector Relative \n Immune Fitness",title = "Scarce Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[1, 2] = Axis(f,xticks = ([1,2,3],["250", "500", "1000"]),xlabel = "Generations Prior to Competition",yticklabelsvisible = false,yticksvisible = false,title = "Alternating Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[1, 3] = Axis(f,xticks = ([1,2,3],["250", "500", "1000"]),yticklabelsvisible = false,yticksvisible = false,title = "Plentiful Resources",limits = lims,xgridvisible = false,ygridvisible = false)
for i in 1:3
    tmp1 = NanMean(SynchBal_Host_Fit_TE[(3*(i-1)+1)][:,:,1].-SynchBal_Host_Fit_OE[(3*(i-1)+1)][:,:,1],1)
    tmp2 = NanMean(SynchBal_Host_Fit_TE[(3*(i-1)+2)][:,:,1].-SynchBal_Host_Fit_OE[(3*(i-1)+2)][:,:,1],1)
    tmp3 = NanMean(SynchBal_Host_Fit_TE[(3*(i-1)+3)][:,:,1].-SynchBal_Host_Fit_OE[(3*(i-1)+3)][:,:,1],1)
    a = lines!(ga[1,i],[1,2,3],mean.([tmp1,tmp2,tmp3]), color = :black, linewidth = 7)
    errorbars!(ga[1,i],[1,2,3],mean.([tmp1,tmp2,tmp3]),(std.([tmp1,tmp2,tmp3])./sqrt(length(tmp1))).*1.96,whiskerwidth = 15, color = :Black) #multiply by 1.96 to convert SEM to 95% Confidence Interval
    
    tmp4 = NanMean(AsynchBal_Host_Fit_TE[(3*(i-1)+1)][:,:,1].-AsynchBal_Host_Fit_OE[(3*(i-1)+1)][:,:,1],1)
    tmp5 = NanMean(AsynchBal_Host_Fit_TE[(3*(i-1)+2)][:,:,1].-AsynchBal_Host_Fit_OE[(3*(i-1)+2)][:,:,1],1)
    tmp6 = NanMean(AsynchBal_Host_Fit_TE[(3*(i-1)+3)][:,:,1].-AsynchBal_Host_Fit_OE[(3*(i-1)+3)][:,:,1],1)
    b = lines!(ga[1,i],[1,2,3],mean.([tmp4,tmp5,tmp6]), color = :black, linestyle = :dot, linewidth = 7)
    errorbars!(ga[1,i],[1,2,3],mean.([tmp4,tmp5,tmp6]),(std.([tmp4,tmp5,tmp6])./sqrt(length(tmp3))).*1.96,whiskerwidth = 15, color = :black)
    labels = ["Synchronous", "Asynchronous"]
    title = "Signaling Timing"
    # elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    Legend(f[1,2], [a,b], labels, title)
end
lines!(ga[1,1],[1,3],[0,0], color = :black, linewidth = 2, linestyle = :dash)
lines!(ga[1,2],[1,3],[0,0], color = :black, linewidth = 2, linestyle = :dash)
lines!(ga[1,3],[1,3],[0,0], color = :black, linewidth = 2, linestyle = :dash)

fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)

cd(string(topDir,"/Images"))
# save(string("Competition_Immune_Cost_Difference_Bal_Fit.png"), f, pt_per_unit = 1)
save("S4_Competition_Imm_Cost_Difference_Bal_Fit.pdf", f)
#endregion

#region mean difference in par cost between TE and OE
# parasite cost is sheet 2 of the given Host fit array (parasite cost for Bal_Host_Fit_TE[1] is [:,:,2])
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (2000, 1000))
lims = (.9,3.1,-.2,.2)
colors = [:Blue, :Orange]
ga = f[1, 1] = GridLayout()
ga[1, 1] = Axis(f,xticks = ([1,2,3],["250", "500", "1000"]),ylabel = "Shared Effector Relative \n Parasitic Fitness",title = "Scarce Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[1, 2] = Axis(f,xticks = ([1,2,3],["250", "500", "1000"]),xlabel = "Generations Prior to Competition",yticklabelsvisible = false,yticksvisible = false,title = "Alternating Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[1, 3] = Axis(f,xticks = ([1,2,3],["250", "500", "1000"]),yticklabelsvisible = false,yticksvisible = false,title = "Plentiful Resources",limits = lims,xgridvisible = false,ygridvisible = false)
for i in 1:3
    tmp1 = NanMean(SynchBal_Host_Fit_TE[(3*(i-1)+1)][:,:,2].-SynchBal_Host_Fit_OE[(3*(i-1)+1)][:,:,2],1)
    tmp2 = NanMean(SynchBal_Host_Fit_TE[(3*(i-1)+2)][:,:,2].-SynchBal_Host_Fit_OE[(3*(i-1)+2)][:,:,2],1)
    tmp3 = NanMean(SynchBal_Host_Fit_TE[(3*(i-1)+3)][:,:,2].-SynchBal_Host_Fit_OE[(3*(i-1)+3)][:,:,2],1)
    a = lines!(ga[1,i],[1,2,3],mean.([tmp1,tmp2,tmp3]), color = :black, linewidth = 7)
    errorbars!(ga[1,i],[1,2,3],mean.([tmp1,tmp2,tmp3]),(std.([tmp1,tmp2,tmp3])./sqrt(length(tmp1))).*1.96,whiskerwidth = 15, color = :black)
    
    tmp4 = NanMean(AsynchBal_Host_Fit_TE[(3*(i-1)+1)][:,:,2].-AsynchBal_Host_Fit_OE[(3*(i-1)+1)][:,:,2],1)
    tmp5 = NanMean(AsynchBal_Host_Fit_TE[(3*(i-1)+2)][:,:,2].-AsynchBal_Host_Fit_OE[(3*(i-1)+2)][:,:,2],1)
    tmp6 = NanMean(AsynchBal_Host_Fit_TE[(3*(i-1)+3)][:,:,2].-AsynchBal_Host_Fit_OE[(3*(i-1)+3)][:,:,2],1)
    b = lines!(ga[1,i],[1,2,3],mean.([tmp4,tmp5,tmp6]), color = :black, linewidth = 7, linestyle = :dot)
    errorbars!(ga[1,i],[1,2,3],mean.([tmp4,tmp5,tmp6]),(std.([tmp4,tmp5,tmp6])./sqrt(length(tmp3))).*1.96,whiskerwidth = 15, color = :black)
    labels = ["Synchronous", "Asynchronous"]
    title = "Signaling Timing"
    # elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    Legend(f[1,2], [a,b], labels, title)
end
lines!(ga[1,1],[1,3],[0,0], color = :black, linewidth = 2, linestyle = :dash)
lines!(ga[1,2],[1,3],[0,0], color = :black, linewidth = 2, linestyle = :dash)
lines!(ga[1,3],[1,3],[0,0], color = :black, linewidth = 2, linestyle = :dash)

fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)
cd(string(topDir,"/Images"))
save(string("S5_Competition_Parasite_Cost_Difference_Bal_Fit.pdf"), f)
#endregion

#region mean difference in Dev cost between TE and OE
# Dev cost is sheet 3 of the given Host fit array (immune cost for Bal_Host_Fit_TE[1] is [:,:,3])
f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),
    resolution = (2000, 1000))
lims = (.9,3.1,-.05,.4)
colors = [:Blue, :Orange]
ga = f[1, 1] = GridLayout()
ga[1, 1] = Axis(f,xticks = ([1,2,3],["250", "500", "1000"]),ylabel = "Shared Effector Relative \n Developmental Fitness",title = "Scarce Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[1, 2] = Axis(f,xticks = ([1,2,3],["250", "500", "1000"]),xlabel = "Generations Prior to Competition",yticklabelsvisible = false,yticksvisible = false,title = "Alternating Resources",limits = lims,xgridvisible = false,ygridvisible = false)
ga[1, 3] = Axis(f,xticks = ([1,2,3],["250", "500", "1000"]),yticklabelsvisible = false,yticksvisible = false,title = "Plentiful Resources",limits = lims,xgridvisible = false,ygridvisible = false)
for i in 1:3
    tmp1 = NanMean(SynchBal_Host_Fit_TE[(3*(i-1)+1)][:,:,3].-SynchBal_Host_Fit_OE[(3*(i-1)+1)][:,:,3],1)
    tmp2 = NanMean(SynchBal_Host_Fit_TE[(3*(i-1)+2)][:,:,3].-SynchBal_Host_Fit_OE[(3*(i-1)+2)][:,:,3],1)
    tmp3 = NanMean(SynchBal_Host_Fit_TE[(3*(i-1)+3)][:,:,3].-SynchBal_Host_Fit_OE[(3*(i-1)+3)][:,:,3],1)
    a = lines!(ga[1,i],[1,2,3],mean.([tmp1,tmp2,tmp3]), color = :black, linewidth = 7)
    errorbars!(ga[1,i],[1,2,3],mean.([tmp1,tmp2,tmp3]),(std.([tmp1,tmp2,tmp3])./sqrt(length(tmp1))).*1.96,whiskerwidth = 15, color = :black)
    
    tmp4 = NanMean(AsynchBal_Host_Fit_TE[(3*(i-1)+1)][:,:,3].-AsynchBal_Host_Fit_OE[(3*(i-1)+1)][:,:,3],1)
    tmp5 = NanMean(AsynchBal_Host_Fit_TE[(3*(i-1)+2)][:,:,3].-AsynchBal_Host_Fit_OE[(3*(i-1)+2)][:,:,3],1)
    tmp6 = NanMean(AsynchBal_Host_Fit_TE[(3*(i-1)+3)][:,:,3].-AsynchBal_Host_Fit_OE[(3*(i-1)+3)][:,:,3],1)
    b = lines!(ga[1,i],[1,2,3],mean.([tmp4,tmp5,tmp6]), color = :black, linewidth = 7, linestyle = :dot)
    errorbars!(ga[1,i],[1,2,3],mean.([tmp4,tmp5,tmp6]),(std.([tmp4,tmp5,tmp6])./sqrt(length(tmp3))).*1.96,whiskerwidth = 15, color = :black)
    labels = ["Synchronous", "Asynchronous"]
    title = "Signaling Timing"
    # elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
    Legend(f[1,2], [a,b], labels, title)
end
lines!(ga[1,1],[1,3],[0,0], color = :black, linewidth = 2, linestyle = :dash)
lines!(ga[1,2],[1,3],[0,0], color = :black, linewidth = 2, linestyle = :dash)
lines!(ga[1,3],[1,3],[0,0], color = :black, linewidth = 2, linestyle = :dash)

fontsize_theme = Theme(fontsize = 28)
set_theme!(fontsize_theme)
cd(string(topDir,"/Images"))
# save(string("Competition_Dev_Cost_Difference_Bal_Fit.png"), f, pt_per_unit = 1)
save("S6_Competition_Dev_Cost_Difference_Bal_Fit.pdf", f)
#endregion