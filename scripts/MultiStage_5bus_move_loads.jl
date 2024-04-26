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
include(joinpath(pwd(), "simulation_scripts/simulation_utils.jl"))

### Parsing Args
sys_name = (@__DIR__)*"/../systems_data/dataset_60p_system/5bus_system_Wind.json"

interval = 24 # simulation step interval in Hours
num_periods = 1
horizon = 15*24 # total time periods in hours, including the day-ahead time
steps = 350 # number of steps in the simulation
battery = true 
form = "StorageDispatch"
network_formulation = "StandardPTDFModel"
output_dir = (@__DIR__)*"/5bus_sims/move_loads"

solver = optimizer_with_attributes(
    Xpress.Optimizer, "MIPRELSTOP" => 1e-5, "OUTPUTLOG" => 0, "MAXTIME" => 1000, "THREADS" => 208
)

if !ispath(output_dir)
    mkpath(output_dir)
end

template_uc = get_template_uc(network_formulation, "StorageDispatch")

# function takes sim_name, bus (for LDES location), and load_bus (for load location)
function run_new_simulation(sim_name, bus, load_bus)
    sys_UC = System(sys_name)

    set_units_base_system!(sys_UC, PSY.UnitSystem.SYSTEM_BASE)

    # Move LDES
    ldes = get_component(GenericBattery, sys_UC, "5bus_60_long_duration")
    remove_component!(sys_UC, ldes)
    
    new_bus = get_component(Bus, sys_UC, bus)

    new_bus_num = new_bus.number
    ldes.bus = new_bus
    ldes.name = "LDES_bus$new_bus_num"

    add_component!(sys_UC, ldes)



    loadd = get_component(PowerLoad, sys_UC, "node_d")
    tsd = get_time_series(SingleTimeSeries, loadd, "max_active_power")

    # if the bus to move to has a load, swap it with bus four; otherwise, move the 
    # load from bus four to that node
    if load_bus in ["node_b", "node_c"]
        load_to_switch = get_component(PowerLoad, sys_UC, load_bus)
        ts_switch = get_time_series(SingleTimeSeries, load_to_switch, "max_active_power")
        
        remove_time_series!(sys_UC, SingleTimeSeries, loadd, "max_active_power")
        remove_time_series!(sys_UC, SingleTimeSeries, load_to_switch, "max_active_power")

        loadd.max_active_power = 3.0
        load_to_switch.max_active_power = 4.0

        add_time_series!(sys_UC, load_to_switch, tsd)
        add_time_series!(sys_UC, loadd, ts_switch)

    elseif load_bus == "node_d"
        a = 1
    else
        new_bus = get_component(Bus, sys_UC, load_bus)
        new_bus_num = new_bus.number

        remove_component!(sys_UC, loadd)
        loadd.bus = new_bus
        loadd.name = load_bus
        add_component!(sys_UC, loadd)
        add_time_series!(sys_UC, loadd, tsd)
    end
    
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

buses = ["node_a", "node_b", "node_c", "node_d", "node_e"]

# Run simulations
for (ii, load_bus) in enumerate(buses)
    for (i, ldes_bus) in enumerate(buses)
        run_new_simulation("5bus_Wind_ld$(i)_load$(ii)_PTDF", ldes_bus, load_bus)
    end
end