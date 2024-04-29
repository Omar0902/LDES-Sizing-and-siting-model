using Revise
using PowerSystems
using PowerSimulations
using Dates
using Logging
using DataFrames
using CSV
using StorageSystemsSimulations
using InfrastructureSystems
logger = configure_logging(console_level=Logging.Info)
const PSI = PowerSimulations
const PSY = PowerSystems
const SS = StorageSystemsSimulations
const IS = InfrastructureSystems
using TimeSeries
using Plots
using JuMP
using HiGHS
using Random

sys_path = (@__DIR__)*"/../systems_data/RTS_system_PV.json"

sys = System(sys_path)

include((@__DIR__)*"/template_functions.jl")
sol =#FILL THIS IN #get_results_dict((@__DIR__)*"/RTS_sims", "RTS_PV", data = data, key = key, sys = sys)

sol["dual_flow"] = hcat(sol["dual_flow_line"], sol["dual_flow_tt"][:, 2:end])

end_pt = 8400
metrics = get_metrics_RTS(sol, end_pt, sys = sys)
rds = collect(get_components(RenewableDispatch, sys))
buses = collect(get_components(Bus, sys))
lines = collect(get_components(Line, sys))
tts = collect(get_components(TapTransformer, sys))
line_names = [l.name for l in lines]
tt_names = [tt.name for tt in tts]
bus_list = [buses[i].name for i in 1:length(buses)]

datadf = DataFrame("buses" => bus_list)
datadf[!, "VRE"] = zeros(length(bus_list))
datadf[!, "VRE_int"] = zeros(length(bus_list))
datadf[!, "load"] = zeros(length(bus_list))
datadf[!, "netload"] = zeros(length(bus_list))
datadf[!, "int_cap"] = zeros(length(bus_list))
datadf[!, "100lines"] = zeros(length(bus_list))
datadf[!, "LMP_avg"] = zeros(length(bus_list))
datadf[!, "LMP_std"] = zeros(length(bus_list))
datadf[!, "avg_int_cap"] = zeros(length(bus_list))
datadf[!, "line_cons"] = [String[] for i in 1:length(bus_list)]

bus_map = Dict()
for i in 1:length(bus_list)
    bus_map[bus_list[i]] = i
end

lmp_buses = names(metrics["lmp_df"])
for i in 1:length(lmp_buses)
    datadf[bus_map[lmp_buses[i]], "LMP_avg"] = metrics["lmp_avg"][i]
    datadf[bus_map[lmp_buses[i]], "LMP_std"] = metrics["lmp_std"][i]
end

rd = sol["rd"]
rd_param = sol["rd_parameter"]
load = sol["load_parameter"]

for i in 2:size(rd_param, 2)
    rd_name = names(rd_param)[i]
    rd_comp = get_component(RenewableDispatch, sys, rd_name)
    bus_name = rd_comp.bus.name
    datadf[bus_map[bus_name], "VRE"] += sum(rd_param[:, i])
    datadf[bus_map[bus_name], "VRE_int"] += sum(rd[:, i])
end

for i in 2:size(load, 2)
    load_name = names(load)[i]
    load_comp = get_component(PowerLoad, sys, load_name)
    bus_name = load_comp.bus.name
    datadf[bus_map[bus_name], "load"] += abs.(sum(load[:, i]))
end

datadf[!, "netload"] = datadf[!, "load"] - datadf[!, "VRE"]

transmission_mat = zeros(length(buses), length(buses))

for l in lines
    rating = l.rate .* 100
    arc = l.arc
    bus1 = arc.from.name
    bus2 = arc.to.name
    datadf[bus_map[bus1], "int_cap"] += rating
    datadf[bus_map[bus2], "int_cap"] += rating
    
    fapv_norm = sol["fapv_norm"][!, l.name]
    datadf[bus_map[bus1], "100lines"] += sum(fapv_norm .>= 1)
    datadf[bus_map[bus2], "100lines"] += sum(fapv_norm .>= 1)
    push!(datadf[bus_map[bus1], "line_cons"], l.name)
    push!(datadf[bus_map[bus2], "line_cons"], l.name)

    bus_idx1 = bus_map[bus1]
    bus_idx2 = bus_map[bus2]

    transmission_mat[bus_idx1, bus_idx2] += rating
    transmission_mat[bus_idx2, bus_idx1] += rating
end

for tt in tts
    rating = tt.rate .* 100
    arc = tt.arc
    bus1 = arc.from.name
    bus2 = arc.to.name
    datadf[bus_map[bus1], "int_cap"] += rating
    datadf[bus_map[bus2], "int_cap"] += rating
    
    fapv_norm = sol["fapv_tt_norm"][!, tt.name]
    datadf[bus_map[bus1], "100lines"] += sum(fapv_norm .>= 1)
    datadf[bus_map[bus2], "100lines"] += sum(fapv_norm .>= 1)
    push!(datadf[bus_map[bus1], "line_cons"], tt.name)
    push!(datadf[bus_map[bus2], "line_cons"], tt.name)
    
    bus_idx1 = bus_map[bus1]
    bus_idx2 = bus_map[bus2]

    transmission_mat[bus_idx1, bus_idx2] += rating
    transmission_mat[bus_idx2, bus_idx1] += rating
end

datadf[:, "100lines"] ./= length.(datadf[:, "line_cons"])

for j in 1:length(bus_list)
    flow_vec = zeros(8400)
    total_cap = [0.]
    line_con_list = datadf[j, "line_cons"]
    for i in line_con_list
        if i in line_names
            flow_vec .+= sol["fapv_abs"][:, i]
            line_comp = get_component(Line, sys, i)
            total_cap[1] += line_comp.rate .* 100
        elseif i in tt_names
            flow_vec .+= sol["fapv_tt_abs"][:, i]
            line_comp = get_component(TapTransformer, sys, i)
            total_cap[1] += line_comp.rate .* 100
        end
    end
    datadf[j, "avg_int_cap"] = mean(flow_vec ./ total_cap[1])
end


# create the topological metrics

inv_mat = zeros(size(transmission_mat))
for i in 1:size(inv_mat, 1)
    for j in 1:size(inv_mat, 2)
        if transmission_mat[i, j] != 0
            inv_mat[i, j] = 1 / transmission_mat[i, j]
        end
    end
end

#=
using Graphs
g = SimpleGraph(length(buses))

for l in lines
    bus1_idx = bus_map[l.arc.from.name]
    bus2_idx = bus_map[l.arc.to.name]
    add_edge!(g, bus1_idx, bus2_idx)
end

for tt in tts
    bus1_idx = bus_map[tt.arc.from.name]
    bus2_idx = bus_map[tt.arc.to.name]
    add_edge!(g, bus1_idx, bus2_idx)
end

abc = closeness_centrality(g, inv_mat)

datadf[!, "Closeness_Centrality"] .= abc

=#

#= Optionally plot with PlasmoData
using PlasmoData, PlasmoDataPlots

dg = DataGraph()

for i in 1:length(buses)
    add_node!(dg, i)
end

for l in lines
    bus1_idx = bus_map[l.arc.from.name]
    bus2_idx = bus_map[l.arc.to.name]
    add_edge!(dg, bus1_idx, bus2_idx)
end

for tt in tts
    bus1_idx = bus_map[tt.arc.from.name]
    bus2_idx = bus_map[tt.arc.to.name]
    add_edge!(dg, bus1_idx, bus2_idx)
end

for i in 1:length(buses)
    for j in 1:length(buses)
        if transmission_mat[i,j] > 0
            add_edge_data!(dg, i, j, transmission_mat[i,j], "transmission_cap")
        end
    end
end


for i in 1:length(buses)
    add_node_data!(dg, i, datadf[i, "Closeness_Centrality"], "cc")
    add_node_data!(dg, i, datadf[i, "LMP_avg"], "LMP_avg")
    add_node_data!(dg, i, datadf[i, "LMP_std"], "LMP_std")
    add_node_data!(dg, i, datadf[i, "100lines"], "100lines")
    add_node_data!(dg, i, datadf[i, "load"], "load")
    add_node_data!(dg, i, (datadf[i, "VRE"] - datadf[i, "VRE_int"]), "curtailment")
    add_node_data!(dg, i, datadf[i, "netload"], "netload")
end

plot_graph(dg, node_z = get_node_data(dg, "cc"), nodecolor = :grays, linewidth = 3, nodesize = 10, legend = true, rev = true, save_fig = true, fig_name = (@__DIR__)*"/cc_metric.png")
plot_graph(dg, node_z = get_node_data(dg, "LMP_avg"), nodecolor = :grays, linewidth = 3, nodesize = 10, legend = true, rev = true, save_fig = true, fig_name = (@__DIR__)*"/lmp_avg.png")
plot_graph(dg, node_z = get_node_data(dg, "LMP_std"), nodecolor = :grays, linewidth = 3, nodesize = 10, legend = true, rev = true)
plot_graph(dg, node_z = get_node_data(dg, "100lines"), nodecolor = :grays, linewidth = 3, nodesize = 10, legend = true, rev = true)
plot_graph(dg, node_z = get_node_data(dg, "load"), nodecolor = :grays, linewidth = 3, nodesize = 10, legend = true, rev = true)
plot_graph(dg, node_z = get_node_data(dg, "curtailment"), nodecolor = :grays, linewidth = 3, nodesize = 10, legend = true, rev = true)
plot_graph(dg, node_z = get_node_data(dg, "netload"), nodecolor = :grays, linewidth = 3, nodesize = 10, legend = true, rev = true)



