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
include((@__DIR)*"/simulation_utils.jl")

sim_name = "RTS_PV"

sys_name = (@__DIR__)*"/../systems_data/dataset_60p_system/RTS_system_PV.json"

interval = 24
num_periods = 1
horizon = 4 * 24
steps = 350
battery = true 
form = "StorageDispatch"
network_formulation = "StandardPTDFModel"
output_dir = (@__DIR__)*"/RTS_sims"

solver = optimizer_with_attributes(
    Xpress.Optimizer, "MIPRELSTOP" => 1e-5, "OUTPUTLOG" => 0, "MAXTIME" => 2000, "THREADS" => Threads.nthreads()
)

### Simulation Setup
if !ispath(output_dir)
    mkpath(output_dir)
end

template_uc = get_template_uc(network_formulation, "StorageDispatch")

sys_UC = PSY.System(sys_name); 

ldes = collect(get_components(GenericBattery, sys_UC))

# Uncomment this section to move LDES to other bus
#=
ldes = first(get_components(GenericBattery, sys_UC)) #(there is only one LDES in RTS system
new_bus = get_component(Bus, sys_UC, "Chifa")
ldes.bus = new_bus
=#


if occursin("RTS", sys_name)
    convert_must_run_units!(sys_UC)
end

set_units_base_system!(sys_UC, PSY.UnitSystem.SYSTEM_BASE)

Random.seed!(10)
add_single_time_series_forecast_error!(sys_UC, horizon, Hour(interval), 0.05)
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
set_system!(results_uc, sys_UC);
export_results_csv(results_uc, "UC", joinpath(results.path, "results"))
