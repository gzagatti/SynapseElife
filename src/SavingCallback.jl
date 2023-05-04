# the objective is to have the savings callback save exactly at the same time as the integrator
mutable struct _SavingAffect{SaveFunc,tType,savevalType,saveatType,saveatCacheType}
    save_func::SaveFunc
    saved_values::SavedValues{tType,savevalType}
    saveat::saveatType
    saveat_cache::saveatCacheType
    save_everystep::Bool
    save_start::Bool
    save_end::Bool
    saveiter::Int
end

function (affect!::_SavingAffect)(integrator, force_save = false)
    just_saved = false
    # see OrdinaryDiffEq.jl -> integrator_utils.jl, function savevalues!
    while !isempty(affect!.saveat) &&
        integrator.tdir * first(affect!.saveat) <= integrator.tdir * integrator.t # Perform saveat
        affect!.saveiter += 1
        curt = pop!(affect!.saveat) # current time
        if curt != integrator.t # If <t, interpolate
            if integrator isa SciMLBase.AbstractODEIntegrator
                # Expand lazy dense for interpolation
                DiffEqBase.addsteps!(integrator)
            end
            if !DiffEqBase.isinplace(integrator.sol.prob)
                curu = integrator(curt)
            else
                curu = first(get_tmp_cache(integrator))
                integrator(curu, curt) # inplace since save_func allocates
            end
            copyat_or_push!(affect!.saved_values.t, affect!.saveiter, curt)
            copyat_or_push!(
                affect!.saved_values.saveval,
                affect!.saveiter,
                affect!.save_func(curu, curt, integrator),
                Val{false},
            )
        else # ==t, just save
            just_saved = true
            copyat_or_push!(affect!.saved_values.t, affect!.saveiter, integrator.t)
            copyat_or_push!(
                affect!.saved_values.saveval,
                affect!.saveiter,
                affect!.save_func(integrator.u, integrator.t, integrator),
                Val{false},
            )
        end
    end
    if !just_saved && affect!.save_everystep ||
       force_save ||
       (
           affect!.save_end &&
           affect!.saved_values.t[affect!.saveiter] != integrator.t &&
           integrator.t == integrator.sol.prob.tspan[end]
       )
        affect!.saveiter += 1
        copyat_or_push!(affect!.saved_values.t, affect!.saveiter, integrator.t)
        copyat_or_push!(
            affect!.saved_values.saveval,
            affect!.saveiter,
            affect!.save_func(integrator.u, integrator.t, integrator),
            Val{false},
        )
    end
    u_modified!(integrator, false)
end

function _saving_initialize!(cb, u, t, integrator)
    integrator.p.xd .= integrator.p.xd0
    cb.affect!.saveat = deepcopy(integrator.opts.saveat)
    cb.affect!.save_everystep = integrator.opts.save_everystep
    cb.affect!.save_start = integrator.opts.save_start
    cb.affect!.save_end = integrator.opts.save_end
    cb.affect!.saveiter = 0
    cb.affect!.save_start && cb.affect!(integrator, true)
end

function _SavingCallback(save_func, saved_values::SavedValues; save_modified = true)
    # adapted from DiffEqCallbacks.jl/src/saving.jl
    saveat_internal =
        DataStructures.BinaryHeap{eltype(saved_values.t)}(DataStructures.FasterForward())
    affect! = _SavingAffect(
        save_func,
        saved_values,
        saveat_internal,
        nothing,
        false,
        false,
        false,
        0,
    )
    condition = if save_modified
        function (u, t, integrator)
            if integrator.u_modified
                push!(affect!.saveat, t)
            end

            return true
        end
    else
        function (u, t, integrator)
            return true
        end
    end
    DiscreteCallback(
        condition,
        affect!;
        initialize = _saving_initialize!,
        save_positions = (false, false),
    )
end
