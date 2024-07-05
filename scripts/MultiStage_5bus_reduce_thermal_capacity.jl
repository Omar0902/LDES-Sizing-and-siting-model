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
using JuMP
using HiGHS
using Xpress
using Random

Random.seed!(10)
include((@__DIR__)*"/simulation_utils.jl")

### Parsing Args
sys_name = (@__DIR__)*"/../systems_data/5bus_system_Wind_caseB.json"

interval = 24 # simulation step interval in Hours
num_periods = 1
horizon = 15*24 # total time periods in hours, including the day-ahead time
steps = 350 # number of steps in the simulation
battery = true
form = "StorageDispatch"
network_formulation = "StandardPTDFModel"
output_dir = (@__DIR__)*"/5bus_sims/reduce_thermal_capacity"

solver = optimizer_with_attributes(
    Xpress.Optimizer, "MIPRELSTOP" => 1e-5, "OUTPUTLOG" => 0, "MAXTIME" => 1000, "THREADS" => 208
)

if !ispath(output_dir)
    mkpath(output_dir)
end

template_uc = get_template_uc(network_formulation, "StorageDispatch")


# function takes sim_name, bus (for LDES location), and ts_restriction (amount of reduction in thermal standard)
function run_new_simulation(sim_name, bus, ts_restriction)
    sys_UC = System(sys_name)

    # Reduce the thermal generation capacity for generator on bus 3
    ts = get_component(ThermalStandard, sys_UC, "Solitude")
    ts.active_power_limits = (min = ts.active_power_limits.min, max = ts_restriction)
    ts.reactive_power_limits = (min = ts.reactive_power_limits.min, max = ts_restriction)

    # Uncomment the lines below (and comment the lines above) to apply the restriction to all generators
    #=
    ts = collect(get_components(ThermalStandard, sys_UC))

    for i in 1:length(ts)
        ts[i].active_power_limits = (min = ts[i].active_power_limits.min, max = ts_restriction)
        ts[i].reactive_power_limits = (min = ts[i].reactive_power_limits.min, max = ts_restriction)
    end
    =#

    # make sure to set the system base AFTER changing the thermal capacity
    set_units_base_system!(sys_UC, PSY.UnitSystem.SYSTEM_BASE)

    # Move LDES
    ldes = get_component(GenericBattery, sys_UC, "5bus_60_long_duration")

    new_bus = get_component(Bus, sys_UC, bus)

    new_bus_num = new_bus.number
    ldes.bus = new_bus
    ldes.name = "LDES_bus$new_bus_num"

    # Add forecast errors for simulations
    Random.seed!(10)
    add_single_time_series_forecast_error!(sys_UC, horizon, Hour(interval), 0.05)

    # Define models and simulations
    models = SimulationModels(
        decision_models=[
            DecisionModel(template_uc,
                sys_UC,
                name="UC",
                optimizer=solver,
                initialize_model=false,
                optimizer_solve_log_print=false,
                direct_mode_optimizer=true,
                check_numerical_bounds=false,
                warm_start=true,
            ),
        ],
    )

    sequence = SimulationSequence(
        models=models,
        ini_cond_chronology=InterProblemChronology(),
    )

    sim = Simulation(
        name="$(sim_name)",
        steps=steps,
        models=models,
        sequence=sequence,
        simulation_folder=output_dir,
        # initial_time=DateTime("2024-01-01T00:00:00"),
    )

    build!(sim, console_level=Logging.Info, file_level=Logging.Info)
    model = get_simulation_model(sim, :UC)
    if occursin("RTS", sys_name)
        add_must_run_constraint!(model)
    end

    exec_time = @elapsed execute!(sim, enable_progress_bar=true, cache_size_mib = 512, min_cache_flush_size_mib = 100)


    results = SimulationResults(sim);
    results_uc = get_decision_problem_results(results, "UC");
    export_results_csv(results_uc, "UC", joinpath(results.path, "results"))
end

# loop through buses with different battery location
#
# run_new_simulation(sim_name, new_bus_num; transmission_line = [], transmission_factor = 0.5)
a = 1
lines = ["line_ab", "line_ad", "line_ae", "line_bc", "line_cd", "line_de"]
tfs = [0.75, 0.5, 0]
sim_name = "5bus_PV_battery_PTDF"
buses = ["node_a", "node_b", "node_c", "node_d", "node_e"]

ts_restrictions = [0.95, 0.9, 0.85]

for ts_restriction in ts_restrictions
    for (i, ldes_bus) in enumerate(buses)
        run_new_simulation("5bus_Wind_ld$(i)_ts$(ts_restriction)_PTDF", ldes_bus, ts_restriction)
    end
end
