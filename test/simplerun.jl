# using Revise, Plots
using Test
using Synapse
using OrdinaryDiffEq, Sundials, LSODA # solvers
using PiecewiseDeterministicMarkovProcesses, JumpProcesses, DiffEqCallbacks
const PDMP = PiecewiseDeterministicMarkovProcesses

##### Parameters
p_synapse = SynapseParams(t_end = 1000.0);
glu = 0.0;
events_sorted_times = [500.0];
is_pre_or_post_event = [true];
events_bap = events_sorted_times[is_pre_or_post_event.==false];
bap_by_epsp = Float64[];
t1 = 0.0;
t2 = 500.0;
nu = buildTransitionMatrix();

##### Initial conditions
xc0 = initial_conditions_continuous_temp(p_synapse);
xd0 = initial_conditions_discrete(p_synapse);

##### Jump problem
jsave_positions = (false, true);
jsaveat = 1 / p_synapse.sampling_rate;

# jalgos = (Tsit5(), Tsit5());
# jalgos = (TRBDF2(), TRBDF2());
# jalgos = (lsoda(), lsoda());
# jalgos = (CVODE_BDF(), CVODE_BDF());
jalgos = (AutoTsit5(Rosenbrock23()), AutoTsit5(Rosenbrock23()));

jumps = J_synapse(p_synapse, nu);
p = (xd0 = copy(xd0), xd = copy(xd0), Glu = p_synapse.glu_amp * glu, p_synapse = p_synapse)
oprob = ODEProblem(
    (du, u, p, t) -> F_synapse(du, u, p.xd, p.p_synapse, t, events_bap, bap_by_epsp),
    xc0,
    (t1, t2),
    p,
);
xdsol = SavedValues(typeof(t1), typeof(xd0));
cb = Synapse._SavingCallback((u, t, integrator) -> integrator.p.xd[:], xdsol);
dep_graph = buildRxDependencyGraph(nu);

# Coevolve
coagg = CoevolveSynced();

coprob = JumpProblem(
    oprob,
    coagg,
    jumps;
    dep_graph = dep_graph,
    save_positions = jsave_positions,
    callback = cb,
    saveat = jsaveat,
);
cosol = @time (
    xcsol = solve(coprob, jalgos[1]; saveat = jsaveat, save_everystep = false),
    xdsol = xdsol,
);

@test length(cosol.xcsol.t) == length(cosol.xdsol.t)

@info "Integrator" cosol.xcsol.stats
# with tweaked JumpProcesses.jl
# total_jumps = (coprob.discrete_jump_aggregation.rejections + coprob.discrete_jump_aggregation.jumps);
# rejections = coprob.discrete_jump_aggregation.rejections[total_jumps .> 0];
# total_jumps = total_jumps[total_jumps .> 0];
# rejection_rate = sum(rejections ./ total_jumps) ./ length(total_jumps);
# @info "Rejection rate" rejection_rate

coresult = @time evolveSynapse(
    xc0,
    xd0,
    p_synapse,
    events_sorted_times,    # external events
    is_pre_or_post_event,   # pre or post?
    bap_by_epsp,
    [true],
    nu,
    jalgos,
    coagg;
    save_positions = jsave_positions,
    saveat = jsaveat,
    save_everystep = false,
);

@test length(coresult.t) == length(coresult.XD)

# @test sum(coresult.XD[end]) == 278

# # Replicate PDMP stepper
# jumps = J_synapse(p_synapse, nu);
# p = (xd0 = copy(xd0), xd = copy(xd0), Glu = p_synapse.glu_amp * glu, p_synapse = p_synapse)
# dep_graph = buildRxDependencyGraph(nu);
# tstops = collect(t1:0.05:t2);

# function mystepper(xc0, xd0, jumps, p, dep_graph, t1, t2, tstops)

#   sol = (xc = VectorOfArray([xc0]), xd = VectorOfArray([xd0]), t = [0.])
#   p.xd .= p.xd0
#   oprob = ODEProblem((du, u, p, t) -> G_synapse(du, u, p.xd, p.p_synapse, t, events_bap, bap_by_epsp), xc0, (t1, t2), p; reltol=1e-7, abstol=1e-9, tstops=tstops);
#   coagg = CoevolveSynced();
#   coprob = JumpProblem(oprob, coagg, jumps; dep_graph, save_positions = (false, false), reltol=1e-7, abstol=1e-9, tstops=tstops);
#   integrator = init(coprob, jalgos[1]; save_everystep=false, advance_to_tstop=true)

#   while (integrator.t < t2)
#     step!(integrator)
#     if integrator.t in tstops
#       push!(sol.xc, integrator.u)
#       push!(sol.xd, integrator.p.xd)
#       push!(sol.t, integrator.t)
#     end
#   end

#   return integrator

# end


# Profile.init(delay=1e-9)
# Profile.clear()
# coresult = @profile evolveSynapse(
# 	xc0,
# 	xd0,
# 	p_synapse,
# 	events_sorted_times,    # external events
# 	is_pre_or_post_event,   # pre or post?
# 	bap_by_epsp,
# 	[true],
# 	nu,
# 	jalgos,
# 	coagg;
# 	save_positions = jsave_positions,
# );
# ProfileSVG.save("./coevolve-synced-profile.svg",
# 	title="evolveSynpase CoevolveSynced", maxdepth=110)
# Profile.clear()


# too slow on this branch
# # Direct
# diagg = Direct();

# diprob = JumpProblem(oprob, diagg, jumps; dep_graph = dep_graph, save_positions = jsave_positions, callback=cb, saveat=jsaveat);
# disol = (xcsol = solve(diprob, jalgos[1]; save_everystep = false), xdsol = xdsol);
# disol = @time (xcsol = solve(diprob, jalgos[1]; save_everystep = false), xdsol = xdsol);

# # unstable
# diresult = @time evolveSynapse(
# 	xc0,
# 	xd0,
# 	p_synapse,
# 	events_sorted_times,    # external events
# 	is_pre_or_post_event,   # pre or post?
# 	bap_by_epsp,
# 	[true],
# 	nu,
# 	jalgos,
# 	diagg;
# 	save_positions = jsave_positions,
# 	save_everystep = false,
# );

#### PDMP problem
pdmpsave_positions = (false, true);

pdmpagg = nothing;

# pdmpalgos = (CHV(:lsoda), CHV(:lsoda));
# pdmpalgos = (CHV(CVODE_BDF()), CHV(CVODE_BDF()));
pdmpalgos = (CHV(AutoTsit5(Rosenbrock23())), CHV(AutoTsit5(Rosenbrock23())));

pdmpprob = PDMP.PDMPProblem(
    (dxc, xc, xd, p, t) -> F_synapse(dxc, xc, xd, p, t, events_bap, bap_by_epsp),
    (rate, xc, xd, p, t, sum_rate) -> R_synapse(rate, xc, xd, p, t, sum_rate, glu),
    nu,
    xc0,
    xd0,
    p_synapse,
    (t1, t2);
    Ncache = 12,
);
pdmpsol = @time solve(pdmpprob, pdmpalgos[1]);

pdmpresult = @time evolveSynapse(
    xc0,
    xd0,
    p_synapse,
    events_sorted_times,    # external events
    is_pre_or_post_event,   # pre or post?
    bap_by_epsp,
    [true],
    nu,
    pdmpalgos,
    pdmpagg;
    save_positions = pdmpsave_positions,
    save_everystep = false,
);

# @test sum(pdmpresult.XD[end])-maximum(pdmpresult.XD[end]) == 278

# @test ~isnothing(result);

##### Plots
using Plots

plot(coresult.t, coresult.XD[1, :], label = "CoevolveSynced");
# plot!(diresult.t, diresult.XD[1, :], label="Direct")
plot!(pdmpresult.t, pdmpresult.XD[1, :], label = "PDMP");
title!("N_ampa")

plot(coresult.t, coresult.XC[1, :], label = "CoevolveSynced");
# plot!(diresult.t, diresult.XC[1, :], label="Direct")
plot!(pdmpresult.t, pdmpresult.XC[1, :], label = "PDMP");
title!("Vsp")

# # plot the discrete variables
# Synapse.plot_discrete(result.t, result.XC, result.XD)
#
# # plot specific variable
# Synapse.plot_variable(result.t, result.XC, result.XD, :Vsp; xlim = (480., 520))
