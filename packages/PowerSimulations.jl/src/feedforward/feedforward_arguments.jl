function add_feedforward_arguments!(
    container::OptimizationContainer,
    model::DeviceModel,
    devices::IS.FlattenIteratorWrapper{V},
) where {V <: PSY.Component}
    for ff in get_feedforwards(model)
        @debug "arguments" ff V _group = LOG_GROUP_FEEDFORWARDS_CONSTRUCTION
        _add_feedforward_arguments!(container, model, devices, ff)
    end
    return
end

function add_feedforward_arguments!(
    container::OptimizationContainer,
    model::ServiceModel,
    service::V,
) where {V <: PSY.AbstractReserve}
    for ff in get_feedforwards(model)
        @debug "arguments" ff V _group = LOG_GROUP_FEEDFORWARDS_CONSTRUCTION
        contributing_devices = get_contributing_devices(model)
        _add_feedforward_arguments!(container, model, contributing_devices, ff)
    end
    return
end

function _add_feedforward_arguments!(
    container::OptimizationContainer,
    model::DeviceModel,
    devices::IS.FlattenIteratorWrapper{T},
    ff::AbstractAffectFeedforward,
) where {T <: PSY.Component}
    parameter_type = get_default_parameter_type(ff, T)
    add_parameters!(container, parameter_type, ff, model, devices)
    return
end

function _add_feedforward_arguments!(
    container::OptimizationContainer,
    model::ServiceModel{SR},
    contributing_devices::Vector{T},
    ff::AbstractAffectFeedforward,
) where {T <: PSY.Component, SR <: PSY.AbstractReserve}
    parameter_type = get_default_parameter_type(ff, SR)
    add_parameters!(container, parameter_type, ff, model, contributing_devices)
    return
end

function _add_feedforward_arguments!(
    container::OptimizationContainer,
    model::DeviceModel,
    devices::IS.FlattenIteratorWrapper{T},
    ff::SemiContinuousFeedforward,
) where {T <: PSY.Component}
    parameter_type = get_default_parameter_type(ff, T)
    add_parameters!(container, parameter_type, ff, model, devices)
    add_to_expression!(
        container,
        ActivePowerRangeExpressionUB,
        parameter_type(),
        devices,
        model,
    )
    add_to_expression!(
        container,
        ActivePowerRangeExpressionLB,
        parameter_type(),
        devices,
        model,
    )
    return
end

function _add_feedforward_arguments!(
    container::OptimizationContainer,
    model::DeviceModel,
    devices::IS.FlattenIteratorWrapper{T},
    ff::EnergyTargetFeedforward,
) where {T <: PSY.Component}
    parameter_type = get_default_parameter_type(ff, T)
    add_parameters!(container, parameter_type, ff, model, devices)
    # Enabling this FF requires the addition of an extra variable
    add_variables!(container, EnergyShortageVariable, devices, get_formulation(model)())
    return
end
