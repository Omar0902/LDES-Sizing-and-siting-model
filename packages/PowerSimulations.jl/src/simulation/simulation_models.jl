"""
Stores the OperationProblem definitions to be used in the simulation. When creating the
SimulationModels, the order in which the models are created determines the order on which
the simulation is executed.
"""
mutable struct SimulationModels
    decision_models::Vector{<:DecisionModel}
    emulation_model::Union{Nothing, EmulationModel}

    function SimulationModels(
        decision_models::Vector,
        emulation_model::Union{Nothing, EmulationModel} = nothing,
    )
        all_names = [get_name(x) for x in decision_models]
        emulation_model !== nothing && push!(all_names, get_name(emulation_model))
        model_count =
            if emulation_model === nothing
                length(decision_models)
            else
                length(decision_models) + 1
            end
        if length(Set(all_names)) != model_count
            error("All model names must be unique: $all_names")
        end

        return new(decision_models, emulation_model)
    end
end

function SimulationModels(
    decision_models::DecisionModel,
    emulation_model::Union{Nothing, EmulationModel} = nothing,
)
    return SimulationModels([decision_models], emulation_model)
end

function SimulationModels(;
    decision_models,
    emulation_model::Union{Nothing, EmulationModel} = nothing,
)
    return SimulationModels(decision_models, emulation_model)
end

function get_simulation_model(models::SimulationModels, name::Symbol)
    for model in models.decision_models
        if get_name(model) == name
            return model
        end
    end
    em = models.emulation_model
    if em !== nothing
        if get_name(em) == name
            return em
        end
    end

    error("Model $name not stored in SimulationModels")
end

function get_simulation_model(models::SimulationModels, index::Int)
    n_decision_models = length(get_decision_models(models))
    if index == n_decision_models + 1
        return models.emulation_model
    elseif index <= n_decision_models
        return get_decision_models(models)[index]
    else
        error("Model number $index is invalid")
    end
end

get_decision_models(models::SimulationModels) = models.decision_models
get_emulation_model(models::SimulationModels) = models.emulation_model

function determine_horizons!(models::SimulationModels)
    horizons = OrderedDict{Symbol, Int}()
    for model in models.decision_models
        container = get_optimization_container(model)
        settings = get_settings(container)
        horizon = get_horizon(settings)
        if horizon == UNSET_HORIZON
            sys = get_system(model)
            horizon = PSY.get_forecast_horizon(sys)
            set_horizon!(settings, horizon)
        end
        horizons[get_name(model)] = horizon
    end
    em = models.emulation_model
    if em !== nothing
        horizons[get_name(em)] = 1
    end
    return horizons
end

function determine_intervals(models::SimulationModels)
    intervals = OrderedDict{Symbol, Dates.Millisecond}()
    for model in models.decision_models
        system = get_system(model)
        interval = PSY.get_forecast_interval(system)
        if interval == Dates.Millisecond(0)
            throw(IS.InvalidValue("Interval of model $(get_name(model)) not set correctly"))
        end
        intervals[get_name(model)] = IS.time_period_conversion(interval)
    end
    em = models.emulation_model
    if em !== nothing
        emulator_system = get_system(em)
        emulator_interval = PSY.get_time_series_resolution(emulator_system)
        intervals[get_name(em)] = IS.time_period_conversion(emulator_interval)
    end
    return intervals
end

function determine_resolutions(models::SimulationModels)
    resolutions = OrderedDict{Symbol, Dates.Millisecond}()
    for model in models.decision_models
        system = get_system(model)
        resolution = PSY.get_time_series_resolution(system)
        if resolution == Dates.Millisecond(0)
            throw(
                IS.InvalidValue("Resolution of model $(get_name(model)) not set correctly"),
            )
        end
        resolutions[get_name(model)] = IS.time_period_conversion(resolution)
    end
    em = models.emulation_model
    if em !== nothing
        emulator_system = get_system(em)
        emulator_resolution = PSY.get_time_series_resolution(emulator_system)
        resolutions[get_name(em)] = IS.time_period_conversion(emulator_resolution)
    end
    return resolutions
end

function initialize_simulation_internals!(models::SimulationModels, uuid::Base.UUID)
    for (ix, model) in enumerate(get_decision_models(models))
        info = SimulationInfo(ix, uuid)
        set_simulation_info!(model, info)
    end
    em = get_emulation_model(models)
    if em !== nothing
        ix = length(get_decision_models(models)) + 1
        info = SimulationInfo(ix, uuid)
        set_simulation_info!(em, info)
    end
    return
end

function get_model_names(models::SimulationModels)
    all_names = get_name.(get_decision_models(models))
    em = get_emulation_model(models)
    if em !== nothing
        push!(all_names, get_name(em))
    end
    return all_names
end
