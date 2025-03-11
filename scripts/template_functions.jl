using PowerNetworkMatrices
const PNM = PowerNetworkMatrices
using OrderedCollections
using Statistics

# The scripts below are designed such that a user can load the required data by calling
# data = get_results_dict(path_to_file, name_of_subfolder, data = data, key = key, sys = sys)
# for example: data = get_results_dict((@__DIR__)*"/5bus_sims/move_wind", "5bus_Wind_ld1_wind1_PTDF", data = data, key = key, sys = sys)
# Metrics (such as production cost) can then be computed by calling
# metrics = get_metrics(data, num_time_points)

function abs_flow_variable(data_dict)

    if "fapv" in keys(data_dict)
        fapv_copy = copy(data_dict["fapv"])
        fapv_copy[:, 2:end] .= abs.(fapv_copy[:, 2:end])
        data_dict["fapv_abs"] = fapv_copy

        fapv_norm_copy = copy(data_dict["fapv_abs"])
        for i in 2:size(fapv_norm_copy, 2)

            max_value = maximum(fapv_norm_copy[:, i])
            fapv_norm_copy[:, i] ./= max_value
        end
        data_dict["fapv_norm"] = fapv_norm_copy
    end
end


function abs_flow_variable_RTS(data_dict, sys)
    sys_to_use = sys

    if "fapv" in keys(data_dict)
        fapv_copy = copy(data_dict["fapv"])
        fapv_copy[:, 2:end] .= abs.(fapv_copy[:, 2:end])
        data_dict["fapv_abs"] = fapv_copy

        fapv_norm_copy = copy(data_dict["fapv_abs"])
        for i in 2:size(fapv_norm_copy, 2)
            line_name = names(fapv_norm_copy)[i]
            l = get_component(Line, sys_to_use, line_name)

            max_value = l.rate * 100
            fapv_norm_copy[:, i] ./= max_value
        end
        data_dict["fapv_norm"] = fapv_norm_copy
    end
    if "fapv_tt" in keys(data_dict)
        fapv_tt_copy = copy(data_dict["fapv_tt"])
        fapv_tt_copy[:, 2:end] .= abs.(fapv_tt_copy[:, 2:end])
        data_dict["fapv_tt_abs"] = fapv_tt_copy

        fapv_tt_norm_copy = copy(data_dict["fapv_tt_abs"])
        for i in 2:size(fapv_tt_norm_copy, 2)
            tt_name = names(fapv_tt_norm_copy)[i]
            tt = get_component(TapTransformer, sys_to_use, tt_name)
            max_value = tt.rate * 100
            fapv_tt_norm_copy[:, i] ./= max_value
        end
        data_dict["fapv_tt_norm"] = fapv_tt_norm_copy
    end
end

function psi_ptdf_lmps(cp_duals, flow_duals, sys)
    ptdf = PNM.PTDF(sys)
    λ = Matrix{Float64}(cp_duals[:, propertynames(cp_duals) .!= :DateTime])

    μ = Matrix{Float64}(flow_duals[:, PNM.get_branch_ax(ptdf)])

    buses = get_components(Bus, sys)
    lmps = OrderedDict()

    bus_dict = Dict()
    bus_nums = ptdf.axes[1]
    for (i, num) in enumerate(bus_nums)
        bus_dict[num] = i
    end

    for bus in buses
        bus_idx = bus_dict[get_number(bus)]
        lmps[get_name(bus)] = μ * ptdf[:, bus_idx]
    end
    lmp = λ .+ DataFrames.DataFrame(lmps)
    return lmp[!, sort(propertynames(lmp))]
end

function psi_ptdf_lmps(data_dict, sys)
    cp_duals = data_dict["cp_duals"]

    λ = Matrix{Float64}(cp_duals[:, propertynames(cp_duals) .!= :DateTime])

    flow_duals = data_dict["flow_duals"]

    μ = Matrix{Float64}(flow_duals[:, PNM.get_branch_ax(ptdf)])

    buses = get_components(Bus, sys)
    lmps = OrderedDict()
    for bus in buses
        lmps[get_name(bus)] = μ * ptdf[:, get_number(bus)]
    end
    lmp = λ .+ DataFrames.DataFrame(lmps)
    return lmp[!, sort(propertynames(lmp))]
end

function get_results_dict(path1, path2; key = [], data = [], sys = nothing)

    sim_results = SimulationResults(path1, path2)
    uc_results = get_decision_problem_results(sim_results, "UC")

    data_dict = Dict()

    for (i, d) in enumerate(data)
        if occursin("Parameter", d)
            g = read_realized_parameter
        elseif occursin("Expression", d)
            g = read_realized_expression
        elseif occursin("Constraint", d)
            g = read_realized_dual
        else
            g = read_realized_variable
        end

        try
            data_results = g(uc_results, d)
        catch e
            continue
        else
            data_results = g(uc_results, d)
            data_dict[key[i]] = data_results
        end
    end

    if isnothing(sys)
        abs_flow_variable(data_dict)
    else
        abs_flow_variable_RTS(data_dict, sys)
    end

    return data_dict
end



function get_pcm(data_dict, end_pt = 8400)
    key_list = [i for i in keys(data_dict) if occursin("pcm", i)]
    pcm_val = [0.]

    for i in key_list
        data = data_dict[i]
        pcm_val[1] += sum(Matrix(data[1:end_pt, 2:end]))
    end

    return pcm_val[1]
end

function get_avg_data(data_vector)
    if length(data_vector) % 24 !=  0
        error("data_vector is wrong length")
    end

    daily_vector = zeros(Int(length(data_vector) / 24))
    hourly_vector = zeros(24)

    for i in 1:length(daily_vector)
        daily_vector[i] = mean(data_vector[((i - 1) * 24 + 1):(i * 24)])
    end

    sigmas = Float64[]
    for i in 1:24
        idx_vec = [i + (j - 1) * 24 for j in 1:length(daily_vector)]
        hourly_vector[i] = mean(data_vector[idx_vec])
        sigma = std(data_vector[idx_vec])
        push!(sigmas, sigma)
    end

    return daily_vector, hourly_vector, sigmas
end

function get_avg_week_data(data_vector)

    weekly_vector = zeros(Int(floor(length(data_vector) / (7 * 24))))

    for i in 1:length(weekly_vector)
        weekly_vector[i] = mean(data_vector[((i - 1) * (24 * 7) + 1):(i * 24 * 7)])
    end

    return weekly_vector
end

PV_rd_dict = Dict("node_a" => ["SolarPV3"], "node_b" => String[], "node_c" => ["SolarPV1"], "node_d" => ["SolarPV2"], "node_e" => ["Wind"])
Wind_rd_dict = Dict("node_a" => String[], "node_b" => ["Wind2"], "node_c" => ["SolarPV1"], "node_d" => String[], "node_e" => ["Wind"])

function get_nodal_net_load(data_dict, rd_dict)
    load_data = data_dict["load_parameter"]
    rd_data = data_dict["rd_parameter"]

    nodal_net_load = DataFrame()
    nodal_net_load[!, "DateTime"] = load_data[:, "DateTime"]

    buses = ["node_a", "node_b", "node_c", "node_d", "node_e"]

    for b in buses
        if b in names(load_data)
            nodal_load = load_data[:, b]
        else
            nodal_load = zeros(size(nodal_net_load, 1))
        end

        nodal_rd = zeros(size(nodal_net_load, 1))

        for i in rd_dict[b]
            nodal_rd .+= rd_data[:, rd_dict[b]]
        end

        nodal_net_load[!, b] = (nodal_load .* -1) - nodal_rd
    end

    return nodal_net_load
end

function get_metrics(data_dict, end_pt, sys)

    metrics_dict = Dict()

    sat_lines = sum(Matrix(data_dict["fapv_norm"][1:end_pt, 2:end]) .== 1, dims = 1)
    tot_sat_lines = sum(sat_lines)

    lines_75 = sum(Matrix(data_dict["fapv_norm"][1:end_pt, 2:end]) .>= .75, dims = 1)
    tot_lines_75 = sum(lines_75)

    lines_50 = sum(Matrix(data_dict["fapv_norm"][1:end_pt, 2:end]) .>= .50, dims = 1)
    tot_lines_50 = sum(lines_50)

    lines_10 = sum(Matrix(data_dict["fapv_norm"][1:end_pt, 2:end]) .<= .10, dims = 1)
    tot_lines_10 = sum(lines_10)

    U90_lines = sum(Matrix(data_dict["fapv_norm"][1:end_pt, 2:end]) .>= 0.9, dims = 1)
    U90_tot = sum(U90_lines)

    metrics_dict["sat"] = sat_lines
    metrics_dict["sat_tot"] = tot_sat_lines
    metrics_dict["75"] = lines_75
    metrics_dict["75_tot"] = tot_lines_75
    metrics_dict["50"] = lines_50
    metrics_dict["50_tot"] = tot_lines_50
    metrics_dict["10"] = lines_10
    metrics_dict["10_tot"] = tot_lines_10
    metrics_dict["90"] = U90_lines
    metrics_dict["90_tot"] = U90_tot

    average_pu = mean(Matrix(data_dict["fapv_norm"][1:end_pt, 2:end]), dims = 2)[:]
    pu_daily_avg, pu_hourly_avg, pu_hourly_sigmas = get_avg_data(average_pu)
    pu_weekly_avg = get_avg_week_data(average_pu)

    metrics_dict["pu_avg"] = average_pu
    metrics_dict["pu_yearly_avg"] = mean(average_pu)
    metrics_dict["pu_hourly_sigma"] = pu_hourly_sigmas
    metrics_dict["pu_daily_avg"] = pu_daily_avg
    metrics_dict["pu_hourly_avg"] = pu_hourly_avg
    metrics_dict["pu_weekly_avg"] = pu_weekly_avg
    metrics_dict["pu_sigma"] = std(average_pu)

    lines = [i for i in names(data_dict["fapv_norm"][:, 2:end])]

    for l in lines
        pu_data = data_dict["fapv_norm"][1:end_pt, l]
        d_avg, h_avg, sig = get_avg_data(pu_data)
        metrics_dict[l*"d"] = d_avg
        metrics_dict[l*"h"] = h_avg
        metrics_dict[l*"y"] = mean(pu_data)
        metrics_dict[l*"h_sigma"] = sig
        metrics_dict[l*"sigma"] = std(pu_data)
    end

    pc_val = get_pcm(data_dict, end_pt)
    metrics_dict["pc"] = pc_val

    rd_param = data_dict["rd_parameter"]
    rd_param = abs.(sum(Matrix(rd_param[1:end_pt, 2:end]), dims = 2))

    rd_realized = data_dict["rd"]
    rd_realized = abs.(sum(Matrix(rd_realized[1:end_pt, 2:end]), dims = 2))

    curtailment = rd_param .- rd_realized
    curt_d, curt_h, curt_h_sigma = get_avg_data(curtailment)
    curt_w = get_avg_week_data(curtailment)

    tot_curtailment = sum(rd_param) - sum(rd_realized)

    metrics_dict["vre_int"] = rd_realized
    metrics_dict["curt"] = tot_curtailment
    metrics_dict["curt_d"] = curt_d
    metrics_dict["curt_h"] = curt_h
    metrics_dict["curt_w"] = curt_w
    metrics_dict["ts_tot"] = sum(Matrix(data_dict["ts"][1:end_pt, 2:end]))

    average_ts = mean(Matrix(data_dict["ts"][1:end_pt, 2:end]), dims = 2)[:]
    sum_ts = sum(Matrix(data_dict["ts"][1:end_pt, 2:end]), dims = 2)[:]
    ts_daily_avg, ts_hourly_avg, ts_hourly_sigma = get_avg_data(average_ts)
    ts_weekly_avg = get_avg_week_data(average_ts)

    ts_max_individual = []
    for i in 2:size(data_dict["ts"], 2)
        push!(ts_max_individual, maximum(data_dict["ts"][1:end_pt, i]))
    end

    metrics_dict["ts_avg"] = average_ts
    metrics_dict["ts_sigma"] = std(sum_ts)
    metrics_dict["ts_max"] = maximum(sum_ts)
    metrics_dict["ts_hourly_sigma"] = ts_hourly_sigma
    metrics_dict["ts_max_ind"] = ts_max_individual
    metrics_dict["ts_d"] = ts_daily_avg
    metrics_dict["ts_h"] = ts_hourly_avg
    metrics_dict["ts_w"] = ts_weekly_avg
    metrics_dict["ts_starts"] = sum(Matrix(data_dict["start_ts"][1:end_pt, 2:end]))
    metrics_dict["ts_stops"] = sum(Matrix(data_dict["stop_ts"][1:end_pt, 2:end]))

    if "ev_gb" in keys(data_dict)
        ev_vec = data_dict["ev_gb"][1:end_pt, 2]

        metrics_dict["ev_max"] = sum(ev_vec .== 5700)
        metrics_dict["ev_80"] = sum(ev_vec .>= 5700 * 0.8)
        metrics_dict["ev_50"] = sum(ev_vec .>= 5700 * 0.5)

        d_avg, h_avg, sigma = get_avg_data(ev_vec)
        metrics_dict["ev_d"] = d_avg
        metrics_dict["ev_h"] = h_avg
    end

    bus_set = ["node_a", "node_b", "node_c", "node_d", "node_e"]
    line_set = ["line_ab", "line_ad", "line_ae", "line_bc", "line_cd", "line_de"]

    bus_connections = [1 1 1 0 0 0; 1 0 0 1 0 0 ; 0 0 0 1 1 0; 0 1 0 0 1 1; 0 0 1 0 0 1]
    bus_connections_adj = [1 1 1 0 0 0; -1 0 0 1 0 0 ; 0 0 0 -1 1 0; 0 -1 0 0 -1 1; 0 0 -1 0 0 -1]
    line_caps = [400, 400, 400, 400, 400, 240]

    flow_amounts = DataFrame()
    flow_amounts[!, "DateTime"] = data_dict["fapv_abs"][1:end_pt, "DateTime"]
    bus_caps = [1200, 800, 800, 1040, 640]

    available_import = DataFrame()
    available_export = DataFrame()
    available_import[!, "DateTime"] = data_dict["fapv_abs"][1:end_pt, "DateTime"]
    available_export[!, "DateTime"] = data_dict["fapv_abs"][1:end_pt, "DateTime"]

    for (i, bus) in enumerate(bus_set)
        flows = zeros(end_pt)
        av_imp = zeros(end_pt)
        av_exp = zeros(end_pt)

        for (j, line) in enumerate(line_set)
            fapv_abs = data_dict["fapv_abs"][1:end_pt, line]
            fapv = data_dict["fapv"][1:end_pt, line]

            flows[:] .+= fapv_abs .* bus_connections[i, j]


            if bus_connections_adj[i, j] != 0
                direction = sign(bus_connections_adj[i, j])
                for k in 1:end_pt
                    if fapv[k] * direction > 0
                        av_imp[k] += line_caps[j] - fapv[k] * direction
                    else
                        av_exp[k] += line_caps[j] + fapv[k] * direction
                    end
                end
            end
        end
        flow_amounts[!, bus] = flows
        available_import[!, bus] = av_imp
        available_export[!, bus] = av_exp
    end
    flow_amounts_norm = copy(flow_amounts)
    for i in 2:6
        flow_amounts_norm[:, i] ./= bus_caps[i - 1]
    end

    metrics_dict["av_import"] = available_import
    metrics_dict["av_export"] = available_export

    metrics_dict["nodal_flows"] = flow_amounts
    metrics_dict["nodal_flows_norm"] = flow_amounts_norm
    metrics_dict["nodal_avg"] = mean(Matrix(flow_amounts[:, 2:end]), dims = 1)
    metrics_dict["nodal_sigma"] = std(Matrix(flow_amounts[:, 2:end]), dims = 1)
    metrics_dict["nodal_avg_norm"] = mean(Matrix(flow_amounts_norm[:, 2:end]), dims = 1)
    metrics_dict["nodal_sigma_norm"] = std(Matrix(flow_amounts_norm[:, 2:end]), dims = 1)
    metrics_dict["nodal_flows_max"] = maximum(Matrix(flow_amounts[:, 2:end]), dims = 1)
    metrics_dict["nodal_flows_min"] = minimum(Matrix(flow_amounts[:, 2:end]), dims = 1)
    metrics_dict["nodal_flows_max_norm"] = maximum(Matrix(flow_amounts_norm[:, 2:end]), dims = 1)
    metrics_dict["nodal_flows_min_norm"] = minimum(Matrix(flow_amounts_norm[:, 2:end]), dims = 1)

    if sys == "PV"
        sys_path = (@__DIR__)*"/../systems_data/dataset_60p_system/5bus_system_PV.json"
        sys = System(sys_path)
    elseif sys == "Wind"
        sys_path = (@__DIR__)*"/../systems_data/dataset_60p_system/5bus_system_Wind_caseB.json"
        sys = System(sys_path)
    else
        sys_to_use = sys
    end

    if "dual_cp" in keys(data_dict) && "dual_flow_line" in keys(data_dict)
        cp_duals = data_dict["dual_cp"]
        flow_duals = data_dict["dual_flow_line"]

        LMP_df = psi_ptdf_lmps(cp_duals, flow_duals, sys_to_use)
        metrics_dict["lmp_df"] = LMP_df
        metrics_dict["lmp_avg"] = mean(Matrix(LMP_df[:, :]), dims = 1)
        metrics_dict["lmp_std"] = std(Matrix(LMP_df[:, :]), dims = 1)
        metrics_dict["lmp_med"] = median(Matrix(LMP_df[:, :]), dims = 1)
    end

    return metrics_dict
end


function get_metrics_RTS(data_dict, end_pt, sys)

    metrics_dict = Dict()

    sat_lines = sum(Matrix(data_dict["fapv_norm"][1:end_pt, 2:end]) .== 1, dims = 1)
    tot_sat_lines = sum(sat_lines)

    lines_75 = sum(Matrix(data_dict["fapv_norm"][1:end_pt, 2:end]) .>= .75, dims = 1)
    tot_lines_75 = sum(lines_75)

    lines_50 = sum(Matrix(data_dict["fapv_norm"][1:end_pt, 2:end]) .>= .50, dims = 1)
    tot_lines_50 = sum(lines_50)

    lines_10 = sum(Matrix(data_dict["fapv_norm"][1:end_pt, 2:end]) .<= .10, dims = 1)
    tot_lines_10 = sum(lines_10)

    U90_lines = sum(Matrix(data_dict["fapv_norm"][1:end_pt, 2:end]) .>= 0.9, dims = 1)
    U90_tot = sum(U90_lines)

    metrics_dict["sat"] = sat_lines
    metrics_dict["sat_tot"] = tot_sat_lines
    metrics_dict["75"] = lines_75
    metrics_dict["75_tot"] = tot_lines_75
    metrics_dict["50"] = lines_50
    metrics_dict["50_tot"] = tot_lines_50
    metrics_dict["10"] = lines_10
    metrics_dict["10_tot"] = tot_lines_10
    metrics_dict["90"] = U90_lines
    metrics_dict["90_tot"] = U90_tot

    average_pu = mean(Matrix(data_dict["fapv_norm"][1:end_pt, 2:end]), dims = 2)[:]
    pu_daily_avg, pu_hourly_avg, pu_hourly_sigmas = get_avg_data(average_pu)
    pu_weekly_avg = get_avg_week_data(average_pu)

    metrics_dict["pu_avg"] = average_pu
    metrics_dict["pu_yearly_avg"] = mean(average_pu)
    metrics_dict["pu_hourly_sigma"] = pu_hourly_sigmas
    metrics_dict["pu_daily_avg"] = pu_daily_avg
    metrics_dict["pu_hourly_avg"] = pu_hourly_avg
    metrics_dict["pu_weekly_avg"] = pu_weekly_avg
    metrics_dict["pu_sigma"] = std(average_pu)

    lines = [i for i in names(data_dict["fapv_norm"][:, 2:end])]

    for l in lines
        pu_data = data_dict["fapv_norm"][1:end_pt, l]
        d_avg, h_avg, sig = get_avg_data(pu_data)
        metrics_dict[l*"d"] = d_avg
        metrics_dict[l*"h"] = h_avg
        metrics_dict[l*"y"] = mean(pu_data)
        metrics_dict[l*"sigma"] = std(pu_data)
    end

    pc_val = get_pcm(data_dict, end_pt)
    metrics_dict["pc"] = pc_val

    rd_param = data_dict["rd_parameter"]
    rd_param = abs.(sum(Matrix(rd_param[1:end_pt, 2:end]), dims = 2))

    rd_realized = data_dict["rd"]
    rd_realized = abs.(sum(Matrix(rd_realized[1:end_pt, 2:end]), dims = 2))

    curtailment = rd_param .- rd_realized
    curt_d, curt_h, curt_h_sigma = get_avg_data(curtailment)
    curt_w = get_avg_week_data(curtailment)

    tot_curtailment = sum(rd_param) - sum(rd_realized)

    metrics_dict["curt"] = tot_curtailment
    metrics_dict["curt_d"] = curt_d
    metrics_dict["curt_h"] = curt_h
    metrics_dict["curt_w"] = curt_w
    metrics_dict["ts_tot"] = sum(Matrix(data_dict["ts"][1:end_pt, 2:end]))

    average_ts = mean(Matrix(data_dict["ts"][1:end_pt, 2:end]), dims = 2)[:]
    sum_ts = sum(Matrix(data_dict["ts"][1:end_pt, 2:end]), dims = 2)[:]
    ts_daily_avg, ts_hourly_avg, ts_hourly_sigma = get_avg_data(average_ts)
    ts_weekly_avg = get_avg_week_data(average_ts)

    ts_max_individual = []
    for i in 2:size(data_dict["ts"], 2)
        push!(ts_max_individual, maximum(data_dict["ts"][1:end_pt, i]))
    end

    metrics_dict["ts_avg"] = average_ts
    metrics_dict["ts_sigma"] = std(sum_ts)
    metrics_dict["ts_max"] = maximum(sum_ts)
    metrics_dict["ts_hourly_sigma"] = ts_hourly_sigma
    metrics_dict["ts_max_ind"] = ts_max_individual
    metrics_dict["ts_d"] = ts_daily_avg
    metrics_dict["ts_h"] = ts_hourly_avg
    metrics_dict["ts_w"] = ts_weekly_avg
    metrics_dict["ts_starts"] = sum(Matrix(data_dict["start_ts"][1:end_pt, 2:end]))
    metrics_dict["ts_stops"] = sum(Matrix(data_dict["stop_ts"][1:end_pt, 2:end]))

    if "ev_gb" in keys(data_dict)
        ev_vec = data_dict["ev_gb"][1:end_pt, 2]

        metrics_dict["ev_max"] = sum(ev_vec .== 5700)
        metrics_dict["ev_80"] = sum(ev_vec .>= 5700 * 0.8)
        metrics_dict["ev_50"] = sum(ev_vec .>= 5700 * 0.5)

        d_avg, h_avg, sigma = get_avg_data(ev_vec)
        metrics_dict["ev_d"] = d_avg
        metrics_dict["ev_h"] = h_avg
    end

    if "dual_cp" in keys(data_dict) && "dual_flow" in keys(data_dict)
        sys_to_use = sys

        cp_duals = data_dict["dual_cp"]
        flow_duals = data_dict["dual_flow"]

        LMP_df = psi_ptdf_lmps(cp_duals, flow_duals, sys_to_use)
        metrics_dict["lmp_df"] = LMP_df
        metrics_dict["lmp_avg"] = mean(Matrix(LMP_df[:, :]), dims = 1)
        metrics_dict["lmp_std"] = std(Matrix(LMP_df[:, :]), dims = 1)
        metrics_dict["lmp_med"] = median(Matrix(LMP_df[:, :]), dims = 1)
    end

    return metrics_dict
end


function plot_data(
    data_list,
    data_entry,
    iters,
    labels;
    title = "",
    second_axis = false,
    second_axis_data = nothing,
    second_axis_iters = nothing,
    second_axis_label = nothing,
    second_axis_legend = :none,
    save_fig = false,
    fig_name = (@__DIR__)*"/fig_name.png",
    linewidth = 1,
    legend = :topleft,
    xlabel = nothing,
    ylabel = nothing
)
    plt = plot(iters, data_list[1][data_entry], label = labels[1], linewidth = linewidth, legend = legend, xlabel = xlabel, ylabel = ylabel)
    for i in 2:length(data_list)
        plot!(plt, iters, data_list[2][data_netry], label = labels[i], linewidth = linewidth)
    end

    if title != ""
        title!(plt, title)
    end

    if second_axis
        plot(twinx(), second_axis_iters, second_axis_data, ylabel = second_axis_label, legend = second_axis_legend)
    end

    if save_fig
        savefig(fig_name)
    end

    display(plt)
    return plt
end


function print_vals(metric, data_set; names = ["Bus 1: ", "Bus 2: ", "Bus 3: ", "Bus 4: ", "Bus 5: "])
    println()
    println("Values for metric ", metric, " are:")
    for i in 1:length(data_set)
        println(names[i], "   ", data_set[i][metric])
    end
    println()
end

function get_scalar_metrics(
    data_dict,
    end_pt = 8400;
    interconnect_cap = [1200, 800, 800, 1040, 640],
    nnd = [2.333, 2.5, 2.5, 2.38, 3.0],
    cc = [320, 266.7, 240, 282.4, 218.2],
    wind = true
)

    if wind
        rd_map = Dict("node_a" => String[], "node_b" => ["Wind2"], "node_c" => ["SolarPV1"], "node_d" => String[], "node_e" => ["Wind"])
    else
        rd_map = Dict("node_a" => ["SolarPV3"], "node_b" => String[], "node_c" => ["SolarPV1"], "node_d" => ["SolarPV2"], "node_e" => ["Wind"])
    end

    ts_map = Dict("node_a" => String[], "node_b" => String[], "node_c" => ["Solitude"], "node_d" => String[], "node_e" => ["Brighton"])


    rd = data_dict["rd_parameter"]
    load = data_dict["load_parameter"]
    rd_realized = data_dict["rd"]
    ts = data_dict["ts"]

    buses = ["node_a", "node_b", "node_c", "node_d", "node_e"]
    # want to return metrics by each node/bus
    bus_load = DataFrame()
    bus_load[!, "DateTime"] = rd[1:end_pt, "DateTime"]
    bus_vre = copy(bus_load)
    bus_net_load = copy(bus_load)
    bus_vre_used = copy(bus_load)
    bus_ts  = copy(bus_load)

    metrics = Dict()

    net_load_per_ic = DataFrame()
    tot_load_per_ic = DataFrame()
    vre_per_ic = DataFrame()
    vre_per_avg_load = DataFrame()
    vre_per_avg_net_load = DataFrame()
    net_load_per_nnd = DataFrame()
    tot_load_per_nnd = DataFrame()
    vre_per_nnd = DataFrame()
    net_load_per_cc = DataFrame()
    tot_load_per_cc = DataFrame()
    vre_per_cc = DataFrame()
    std_net_load_per_ic = DataFrame()
    std_net_load_per_nnd = DataFrame()
    std_net_load_per_cc = DataFrame()
    std_load_df = DataFrame()
    std_load_per_cc = DataFrame()
    load_df = DataFrame()
    vre_df = DataFrame()
    net_load_df = DataFrame()
    std_net_load_df = DataFrame()
    std_vre_df = DataFrame()
    tot_prod_df = DataFrame()

    for (i, b) in enumerate(buses)
        if length(rd_map[b]) > 0
            bus_vre[!, b] = rd[1:end_pt, rd_map[b][1]]
            bus_vre_used[!, b] = rd_realized[1:end_pt, rd_map[b][1]]
        else
            bus_vre[!, b] = zeros(Float64, end_pt)
            bus_vre_used[!, b] = zeros(Float64, end_pt)
        end

        if b in names(load)
            bus_load[!, b] = load[1:end_pt, b]
        else
            bus_load[!, b] = zeros(Float64, end_pt)
        end

        if length(ts_map[b]) > 0
            bus_ts[!, b] = ts[1:end_pt, ts_map[b][1]]
        else
            bus_ts[!, b] = zeros(Float64, end_pt)
        end

        bus_net_load[!, b] = bus_load[:, b] .- bus_vre[:, b]

        load_vec = bus_load[:, b]
        net_load_vec = bus_net_load[:, b]
        vre_vec = bus_vre[:, b]
        ts_vec = bus_ts[:, b]
        vre_used_vec = bus_vre[:, b]

        avg_load = mean(load_vec)
        avg_net_load = mean(net_load_vec)
        avg_vre = mean(vre_vec)
        sum_load = sum(load_vec)
        sum_net_load = sum(net_load_vec)
        sum_vre = sum(vre_vec)
        std_net_load = std(net_load_vec)
        std_load = std(load_vec)
        std_vre = std(vre_vec)

        net_load_per_ic[!, b] = [sum_net_load / interconnect_cap[i]]
        tot_load_per_ic[!, b] = [sum_load / interconnect_cap[i]]
        vre_per_ic[!, b] = [sum_vre / interconnect_cap[i]]
        if sum_load == 0
            vre_per_avg_load[!, b] = [0.]
        else
            vre_per_avg_load[!, b] = [sum_vre ./ avg_load]
        end
        if avg_net_load == 0
            vre_per_avg_net_load[!, b] = [0.]
        else
            vre_per_avg_net_load[!, b] = [sum_vre ./ avg_net_load]
        end
        net_load_per_nnd[!, b] = [sum_net_load / nnd[i]]
        tot_load_per_nnd[!, b] = [sum_load / nnd[i]]
        vre_per_nnd[!, b] = [sum_vre / nnd[i]]
        net_load_per_cc[!, b] = [sum_net_load / cc[i]]
        tot_load_per_cc[!, b] = [sum_load / cc[i]]
        vre_per_cc[!, b] = [sum_vre / cc[i]]
        std_net_load_per_ic[!, b] = [std_net_load / interconnect_cap[i]]
        std_net_load_per_nnd[!, b] = [std_net_load / nnd[i]]
        std_net_load_per_cc[!, b] = [std_net_load / cc[i]]
        std_load_df[!, b] = [std_load]
        std_load_per_cc[!, b] = [std_load / cc[i]]
        load_df[!, b] = [avg_load]
        vre_df[!, b] = [avg_vre]
        net_load_df[!, b] = [avg_net_load]
        std_net_load_df[!, b] = [std_net_load]
        std_vre_df[!, b] = [std_vre]
        tot_prod_df[!, b] = [mean(ts_vec .+ vre_used_vec)]
    end


    metrics["load"] = load_df
    metrics["vre"] = vre_df
    metrics["net_load"] = net_load_df
    metrics["std_net_load"] = std_net_load_df
    metrics["std_vre"] = std_vre_df
    metrics["net_load_per_ic"] = net_load_per_ic
    metrics["load_per_ic"] = tot_load_per_ic
    metrics["vre_per_ic"] = vre_per_ic
    metrics["vre_per_avg_load"] = vre_per_avg_load
    metrics["vre_per_avg_net_load"] = vre_per_avg_net_load
    metrics["net_load_per_nnd"] = net_load_per_nnd
    metrics["load_per_nnd"] = tot_load_per_nnd
    metrics["vre_per_nnd"] = vre_per_nnd
    metrics["net_load_per_cc"] = net_load_per_cc
    metrics["load_per_cc"] = tot_load_per_cc
    metrics["vre_per_cc"] = vre_per_cc
    metrics["std_net_load_per_ic"] = std_net_load_per_ic
    metrics["std_net_load_per_nnd"] = std_net_load_per_nnd
    metrics["std_net_load_per_cc"] = std_net_load_per_cc
    metrics["std_load"] = std_load_df
    metrics["std_load_per_cc"] = std_load_per_cc
    metrics["tot_prod"] = tot_prod_df

    return metrics
end

lines = ["line_ab", "line_ad", "line_ae", "line_bc", "line_cd", "line_de"]

data = ["FlowActivePowerVariable__Line",
        "FlowActivePowerVariable__TapTransformer",
        "ActivePowerTimeSeriesParameter__PowerLoad",
        "ActivePowerTimeSeriesParameter__RenewableDispatch",
        "ActivePowerVariable__ThermalStandard",
        "StartVariable__ThermalStandard",
        "StopVariable__ThermalStandard",
        "OnVariable__ThermalStandard",
        "ActivePowerVariable__RenewableDispatch",
        "ActivePowerOutVariable__BatteryEMS",
        "ActivePowerInVariable__BatteryEMS",
        "EnergyVariable__BatteryEMS",
        "SystemBalanceSlackUp__System",
        "SystemBalanceSlackDown__System",
        "ProductionCostExpression__ThermalStandard",
        "ProductionCostExpression__RenewableDispatch",
        "ProductionCostExpression__BatteryEMS",
        "ProductionCostExpression__GenericBattery",
        "ProductionCostExpression__ThermalMultiStart",
        "ActivePowerOutVariable__GenericBattery",
        "ActivePowerInVariable__GenericBattery",
        "EnergyVariable__GenericBattery",
        "CopperPlateBalanceConstraint__System",
        "NetworkFlowConstraint__Line",
        "NetworkFlowConstraint__TapTransformer"
]
key = ["fapv",
        "fapv_tt",
        "load_parameter",
        "rd_parameter",
        "ts",
        "start_ts",
        "stop_ts",
        "on_ts",
        "rd",
        "out_ems",
        "in_ems",
        "ev_ems",
        "slackup",
        "slackdown",
        "pcm_ts",
        "pcm_rd",
        "pcm_ems",
        "pcm_gb",
        "pcm_tms",
        "out_gb",
        "in_gb",
        "ev_gb",
        "dual_cp",
        "dual_flow_line",
        "dual_flow_tt"
]
