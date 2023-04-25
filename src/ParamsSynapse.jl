"""
$(TYPEDEF)

Presynaptic parameters

- Firing events are processed separately from the main simulation (at
  `src/OnlyStp.jl`) it takes the firing structure as input from the function
  [`firingPattern`](@ref).
- Using the presynaptic stimulation times, the vesicle depletion and AP induced
  by EPSP are estimated, however one can use a tag to deactivate it (covering
  sub-threshold EPSP cases) as in  [`dataProtocol`](@ref).
- The presynaptic part of our model is phenomenological, for instance, the
  variable `Soma` in `src/OnlyStp.jl:38` was made to represent the voltage and
  can accumulate for faster frequencies but has an abstract unit.

# Equations

Based on D. Sterratt et al book; [`Principles of Computational Modelling in Neuroscience`](https://www.compneuroprinciples.org/) page 183

```math
\\begin{aligned}
raterec  &= (`R_0` - `R`) ⋅ τ_R ⋅ `rrp` \\\\
raterrp  &= (`D_0` - `D`) ⋅ τ_D ⋅ `rec` \\\\
rateref  &= (`R_0` - `R`) ⋅ ref_dt \\\\
\\end{aligned}
```

# Fields

$(FIELDS)
"""
@with_kw struct PreSynapseParams
    "recovery constant of pre calcium decay function."
    τ_rec::Float64 = 20000
    "fraction of decay constant of pre calcium decay f."
    δ_ca::Float64 = 0.0004
    "decay time constant of pre calcium."
    τ_pre::Float64 = 20
    "decay time constant for AP induced by EPSP."
    τ_V::Float64 = 40
    "delay to EPSPs onset and evoked AP."
    δ_delay_AP::Float64 = 15.0
    "initial conditions ready releaseble pool."
    D_0::Int64 = 25
    "initial conditions recovery pool."
    R_0::Int64 = 30
    "rate for `D -> R`."
    τ_R::Float64 = 5000
    "rate for `R -> D`."
    τ_D::Float64 = 45000
    "rate for `infinite reservoir -> R`."
    τ_R_ref::Float64 = 40000
    "sigmoid parameter for release probability."
    s::Float64 = 2.0
    "sigmoid parameter for release probability."
    h::Float64 = 0.7 # this value changes given Ca external concentration
    "sampling rate for plotting / printing."
    sampling_rate::Float64 = 1.0
end

"""
$(TYPEDEF)

Postsynaptic parameters.

# Units
- time: ms
- length: µm (area µm^2, volume µm^3)
- voltage: mV
- current: pA
- conductance: nS
- resistance: GOhm
- capacitance: pF ( ms.pA/mV = ms.nS = ms/GOhm)
- concentration: µM

# Fields

$(FIELDS)

# References

[bartol2015] Bartol, T.M., Keller, D.X., Kinney, J.P., Bajaj, C.L., Harris, K.M., Sejnowski,
T.J., Kennedy, M.B., 2015. Computational reconstitution of spine calcium
transients from individual proteins. Frontiers in synaptic neuroscience 7, 17.

[liu1999] Liu, G., Choi, S., Tsien, R.W., 1999. Variability of neurotransmitter
concentration and nonsaturation of postsynaptic ampa receptors at synapses in
hippocampal cultures and slices. Neuron 22, 395–409.

[magee1995] Magee, J.C., Johnston, D., 1995. Characterization of single
voltage-gated na+ and ca2+ channels in apical dendrites of rat ca1 pyramidal
neurons. The Journal of physiology 487, 67–90.

[maylie2004] Maylie, J., Bond, C.T., Herson, P.S., Lee, W.S., Adelman, J.P., 2004. 
Small conductance ca2+-activated k+ channels and calmodulin. The Journal
of physiology 554, 255–261.

[quintana2005] Quintana, A.R., Wang, D., Forbes, J.E., Waxham, M.N., 2005.
Kinetics of calmodulin binding to calcineurin. Biochemical and biophysical
research communications 334, 674–680.

[tigaret2016] Tigaret, C.M., Olivo, V., Sadowski, J.H., Ashby, M.C., Mellor, J.R., 2016.
Coordinated activation of distinct ca 2+ sources and metabotropic glutamate
receptors encodes hebbian synaptic plasticity. Nature communications 7, 10289.

Mellor and Griffith?

Lisman?

Zenisek et al 2003

Naraghi, 1997
"""
@with_kw struct SynapseParams{Tp}
    "polygonal threshold."
    LTP_region::Tp = VPolygon([[6.35, 1.4], [10, 1.4], [6.35, 29.5], [10, 29.5]]) # VPolygon([SVector(6.35,1.4), SVector(10,1.4),SVector(6.35,29.5), SVector(10,29.5)])
    "polygonal threshold."
    LTD_region::Tp = VPolygon([
        [6.35, 1.4],
        [6.35, 23.25],
        [6.35, 29.5],
        [1.85, 11.327205882352938],
        [1.85, 23.25],
        [3.7650354609929075, 1.4],
        [5.650675675675676, 29.5],
    ]) # VPolygon([SVector(6.35,1.4),SVector(6.35,23.25),SVector(6.35,29.5),SVector(1.85,11.327205882352938),SVector(1.85,23.25),SVector(3.7650354609929075,1.4),SVector(5.650675675675676,29.5)])
    "activation rates for thresholds."
    a_D::Float64 = 0.1
    "activation rates for thresholds."
    b_D::Float64 = 0.00002
    "activation rates for thresholds."
    a_P::Float64 = 0.2
    "activation rates for thresholds."
    b_P::Float64 = 0.0001
    "activation rates for thresholds."
    t_D::Float64 = 18000
    "activation rates for thresholds."
    t_P::Float64 = 13000
    "sigmoids controlling the rate of plasticity change."
    K_D::Float64 = 80000.0
    "sigmoids controlling the rate of plasticity change."
    K_P::Float64 = 13000.0
    "plasticity states."
    rest_plstcty::Int64 = 100
    "simulation."
    t_end::Float64 = 100
    "simulation."
    sampling_rate::Float64 = 10
    "biophysical and GHK parameters."
    temp_rates::Float64 = 35.0
    "biophysical and GHK parameters."
    age::Float64 = 60.0
    "biophysical and GHK parameters."
    faraday::Float64 = 96485e-6 * 1e-3
    "biophysical and GHK parameters."
    Ca_ext::Float64 = 2.5e3
    "biophysical and GHK parameters."
    Ca_infty::Float64 = 50e-3
    "biophysical and GHK parameters."
    tau_ca::Float64 = 10.0
    "biophysical and GHK parameters."
    D_Ca::Float64 = 0.3338
    "biophysical and GHK parameters."
    f_Ca::Float64 = 0.1
    "biophysical and GHK parameters."
    perm::Float64 = -0.04583333333333333
    "biophysical and GHK parameters."
    z::Float64 = 2.0
    "biophysical and GHK parameters."
    gas::Float64 = 8.314e-6
    "biophysical and GHK parameters."
    p_release::NTuple{4,Float64} =
        (0.004225803293622208, 1708.4124496514878, 1.3499793762587964, 0.6540248201173222)
    "backpropagation attenuation."
    trec::Float64 = 2000
    "backpropagation attenuation."
    trec_soma::Float64 = 500
    "backpropagation attenuation."
    delta_decay::Float64 = 1.7279e-5
    "backpropagation attenuation."
    p_age_decay_bap::NTuple{3,Float64} =
        (0.13525468256031167, 16.482800452454164, 5.564691354645679)
    "backpropagation attenuation."
    delta_soma::Float64 =
        2.5e-5 *
        (p_age_decay_bap[3] / (1 + exp(p_age_decay_bap[1] * (age - p_age_decay_bap[2]))))
    "backpropagation attenuation."
    delta_aux::Float64 = 2.304e-5
    "backpropagation attenuation."
    injbap::Float64 = 2.0
    "backpropagation attenuation."
    soma_dist::Float64 = 200.0
    "backpropagation attenuation."
    p_dist::NTuple{4,Float64} =
        (0.019719018173341547, 230.3206470553394, 1.4313810030893268, 0.10406540965358434)
    "backpropagation attenuation."
    ϕ_dist::Float64 =
        (p_dist[4] + p_dist[3] / (1 + exp(p_dist[1] * (soma_dist - p_dist[2]))))
    "backpropagation attenuation."
    I_clamp::Float64 = 0.0
    "Na and K."
    gamma_Na::Float64 = 8e2
    "Na and K."
    Erev_Na::Float64 = 50.0
    "Na and K."
    gamma_K::Float64 = 4e1
    "Na and K."
    Erev_K::Float64 = -90.0
    "NMDAr temperature modification."
    p_nmda_frwd::NTuple{4,Float64} =
        (-0.09991802053299291, -37.63132907014948, 1239.0673283348326, -1230.6805720050966)
    "NMDAr temperature modification."
    frwd_T_chng_NMDA::Float64 = (
        p_nmda_frwd[4] +
        p_nmda_frwd[3] / (1 + exp(p_nmda_frwd[1] * (temp_rates - p_nmda_frwd[2])))
    )
    "NMDAr temperature modification."
    p_nmda_bcwd::NTuple{4,Float64} =
        (-0.10605060141396823, 98.99939433046647, 1621.6168608608068, 3.0368551011554143)
    "NMDAr temperature modification."
    bcwd_T_chng_NMDA::Float64 = (
        p_nmda_bcwd[4] +
        p_nmda_bcwd[3] / (1 + exp(p_nmda_bcwd[1] * (temp_rates - p_nmda_bcwd[2])))
    ) # 0.16031*temp_rates - 0.80775
    "NMDAr kinetics (GluN2A type), uM-1ms-1."
    NMDA_N2A_ka::Float64 = frwd_T_chng_NMDA * 34.0 * 1e-3
    "NMDAr kinetics (GluN2A type), uM-1ms-1."
    NMDA_N2A_kb::Float64 = frwd_T_chng_NMDA * 17.0 * 1e-3
    "NMDAr kinetics (GluN2A type), uM-1ms-1."
    NMDA_N2A_kc::Float64 = frwd_T_chng_NMDA * 127.0 * 1e-3
    "NMDAr kinetics (GluN2A type), uM-1ms-1."
    NMDA_N2A_kd::Float64 = frwd_T_chng_NMDA * 580.0 * 1e-3
    "NMDAr kinetics (GluN2A type), uM-1ms-1."
    NMDA_N2A_ke::Float64 = frwd_T_chng_NMDA * 2508.0 * 1e-3
    "NMDAr kinetics (GluN2A type), uM-1ms-1."
    NMDA_N2A_kf::Float64 = frwd_T_chng_NMDA * 3449.0 * 1e-3
    "NMDAr kinetics (GluN2A type), ms-1."
    NMDA_N2A_k_f::Float64 = bcwd_T_chng_NMDA * 662.0 * 1e-3
    "NMDAr kinetics (GluN2A type), ms-1."
    NMDA_N2A_k_e::Float64 = bcwd_T_chng_NMDA * 2167.0 * 1e-3
    "NMDAr kinetics (GluN2A type), ms-1."
    NMDA_N2A_k_d::Float64 = bcwd_T_chng_NMDA * 2610.0 * 1e-3
    "NMDAr kinetics (GluN2A type), ms-1."
    NMDA_N2A_k_c::Float64 = bcwd_T_chng_NMDA * 161.0 * 1e-3
    "NMDAr kinetics (GluN2A type), ms-1."
    NMDA_N2A_k_b::Float64 = bcwd_T_chng_NMDA * 120.0 * 1e-3
    "NMDAr kinetics (GluN2A type), ms-1."
    NMDA_N2A_k_a::Float64 = bcwd_T_chng_NMDA * 60.0 * 1e-3
    "NMDAr kinetics (GluN2B type), uM-1ms-1."
    NMDA_N2B_sa::Float64 = frwd_T_chng_NMDA * 0.25 * 34.0 * 1e-3
    "NMDAr kinetics (GluN2B type), uM-1ms-1."
    NMDA_N2B_sb::Float64 = frwd_T_chng_NMDA * 0.25 * 17.0 * 1e-3
    "NMDAr kinetics (GluN2B type), uM-1ms-1."
    NMDA_N2B_sc::Float64 = frwd_T_chng_NMDA * 0.25 * 127.0 * 1e-3
    "NMDAr kinetics (GluN2B type), uM-1ms-1."
    NMDA_N2B_sd::Float64 = frwd_T_chng_NMDA * 0.25 * 580.0 * 1e-3
    "NMDAr kinetics (GluN2B type), uM-1ms-1."
    NMDA_N2B_se::Float64 = frwd_T_chng_NMDA * 0.25 * 2508.0 * 1e-3
    "NMDAr kinetics (GluN2B type), uM-1ms-1."
    NMDA_N2B_sf::Float64 = frwd_T_chng_NMDA * 0.25 * 3449.0 * 1e-3
    "NMDAr kinetics (GluN2B type), ms-1."
    NMDA_N2B_s_f::Float64 = bcwd_T_chng_NMDA * 0.23 * 662.0 * 1e-3
    "NMDAr kinetics (GluN2B type), ms-1."
    NMDA_N2B_s_e::Float64 = bcwd_T_chng_NMDA * 0.23 * 2167.0 * 1e-3
    "NMDAr kinetics (GluN2B type), ms-1."
    NMDA_N2B_s_d::Float64 = bcwd_T_chng_NMDA * 0.23 * 2610.0 * 1e-3
    "NMDAr kinetics (GluN2B type), ms-1."
    NMDA_N2B_s_c::Float64 = bcwd_T_chng_NMDA * 0.23 * 161.0 * 1e-3
    "NMDAr kinetics (GluN2B type), ms-1."
    NMDA_N2B_s_b::Float64 = bcwd_T_chng_NMDA * 0.23 * 120.0 * 1e-3
    "NMDAr kinetics (GluN2B type), ms-1."
    NMDA_N2B_s_a::Float64 = bcwd_T_chng_NMDA * 0.23 * 60.0 * 1e-3
    "NMDA details."
    p_nmda::NTuple{4,Float64} =
        (0.004477162852447629, 2701.3929349701334, 58.38819453272428, 33.949463268365555)
    "NMDA details."
    gamma_nmda::Float64 =
        (p_nmda[4] + p_nmda[3] / (1 + exp(p_nmda[1] * (Ca_ext - p_nmda[2])))) * 1e-3
    "NMDA details."
    p_age::NTuple{4,Float64} =
        (0.09993657672916968, 25.102347872464193, 0.9642137892004939, 0.5075183905839776)
    "ratio N2B/N2A."
    r_NMDA_age::Float64 =
        rand(Normal(0, 0.05)) + p_age[4] + p_age[3] / (1 + exp(p_age[1] * (age - p_age[2]))) # 0.5+1.6*exp(-age/16.66) + rand(Normal(0,.05))
    "ratio N2B/N2A."
    N_NMDA::Int64 = 15
    "ratio N2B/N2A."
    N_N2B::Int64 = round(N_NMDA * r_NMDA_age / (r_NMDA_age + 1))
    "ratio N2B/N2A, using Sinclair ratio."
    N_N2A::Int64 = round(N_NMDA / (r_NMDA_age + 1))
    "other NMDAr parameters."
    Erev_nmda::Float64 = 0.0
    "other NMDAr parameters."
    Mg::Float64 = 1.0
    "AMPAr temperature modification."
    p_ampa_frwd::NTuple{3,Float64} =
        (-0.4737773089201679, 31.7248285571622, 10.273135485873242)
    "AMPAr temperature modification."
    frwd_T_chng_AMPA::Float64 =
        (p_ampa_frwd[3] / (1 + exp(p_ampa_frwd[1] * (temp_rates - p_ampa_frwd[2])))) # temp_rates*0.78-18.7
    "AMPAr temperature modification."
    p_ampa_bcwd::NTuple{3,Float64} =
        (-0.36705555170278986, 28.976662403966674, 5.134547217640794)
    "AMPAr temperature modification."
    bcwd_T_chng_AMPA::Float64 =
        (p_ampa_bcwd[3] / (1 + exp(p_ampa_bcwd[1] * (temp_rates - p_ampa_bcwd[2])))) # temp_rates*0.37-8.25
    "AMPAr kinetics, uM-1ms-1."
    AMPA_k1::Float64 = frwd_T_chng_AMPA * 1.6 * 1e7 * 1e-6 * 1e-3
    "AMPAr kinetics, ms-1."
    AMPA_k_1::Float64 = bcwd_T_chng_AMPA * 7400 * 1e-3
    "AMPAr kinetics, ms-1."
    AMPA_k_2::Float64 = bcwd_T_chng_AMPA * 0.41 * 1e-3
    "AMPAr kinetics, ms-1."
    AMPA_alpha::Float64 = 2600 * 1e-3
    "AMPAr kinetics, ms-1."
    AMPA_beta::Float64 = 9600 * 1e-3
    "AMPAr kinetics, ms-1."
    AMPA_delta_1::Float64 = 1500 * 1e-3
    "AMPAr kinetics, ms-1."
    AMPA_gamma_1::Float64 = 9.1 * 1e-3
    "AMPAr kinetics, ms-1."
    AMPA_delta_2::Float64 = 170 * 1e-3
    "AMPAr kinetics, ms-1."
    AMPA_gamma_2::Float64 = 42 * 1e-3
    "AMPAr kinetics, ms-1."
    AMPA_delta_0::Float64 = 0.003 * 1e-3
    "AMPAr kinetics, ms-1."
    AMPA_gamma_0::Float64 = 0.83 * 1e-3
    "AMPAr conductances, nS."
    gamma_ampa1::Float64 = 0.5 * 31e-3
    "AMPAr conductances, nS."
    gamma_ampa2::Float64 = 0.5 * 52e-3
    "AMPAr conductances, nS."
    gamma_ampa3::Float64 = 0.5 * 73e-3
    "AMPAr conductances, AMPAr number."
    N_ampa::Int64 = 120
    "AMPAr conductances, AMPAR reversal potential, mV."
    Erev_ampa::Float64 = 0.0
    "GABAr."
    N_GABA::Int64 = 34
    "GABAr."
    p_Cl::NTuple{4,Float64} =
        (0.09151696057098718, 0.6919298240788684, 243.5159017060495, -92.6496083089155)
    "GABAr, Cl reversal potential."
    Erev_Cl::Float64 = (p_Cl[4] + p_Cl[3] / (1 + exp(p_Cl[1] * (age - p_Cl[2]))))
    "GABAr, Cl reversal potential."
    gamma_GABA::Float64 = 35e-3
    "GABAr, Cl reversal potential."
    GABA_r_b1::Float64 = 1e6 * 1e-6 * 1e-3 * 20
    "GABAr, Cl reversal potential."
    GABA_r_u1::Float64 = 1e3 * 4.6e-3
    "GABAr, Cl reversal potential."
    GABA_r_b2::Float64 = 1e6 * 1e-6 * 1e-3 * 10
    "GABAr, Cl reversal potential."
    GABA_r_u2::Float64 = 1e3 * 9.2e-3
    "GABAr, Cl reversal potential."
    GABA_r_ro1::Float64 = 1e3 * 3.3e-3
    "GABAr, Cl reversal potential."
    GABA_r_ro2::Float64 = 1e3 * 10.6e-3
    "GABAr, Cl reversal potential."
    p_GABA::NTuple{4,Float64} =
        (0.19127068198185954, 32.16771140618756, -1.2798050197287802, 1.470692263981145)
    "GABAr, Cl reversal potential."
    GABA_r_c1::Float64 =
        (p_GABA[4] + p_GABA[3] / (1 + exp(p_GABA[1] * (temp_rates - p_GABA[2])))) *
        1e3 *
        9.8e-3
    "GABAr, Cl reversal potential."
    GABA_r_c2::Float64 =
        (p_GABA[4] + p_GABA[3] / (1 + exp(p_GABA[1] * (temp_rates - p_GABA[2])))) * 400e-3
    "passive electrical properties."
    E_leak::Float64 = -70.0
    "passive electrical properties."
    g_leak::Float64 = 4e-6
    "passive electrical properties."
    Cm::Float64 = 0.6e-2
    "passive electrical properties."
    R_a::Float64 = 1e-2
    "morphology, Dendritic properties, dendrite diameter, um."
    D_dend::Float64 = 2.0
    "morphology, Dendritic properties, dendrite length, choosen to tune attenuation, but not modified in BaP adaptation for simplicity sake, um."
    L_dend::Float64 = 1400
    "morphology, Dendritic properties, dendrite surface area, 500 gives dendrite input resistance of 200MOhm, um^2."
    A_dend::Float64 = 2 * pi * (D_dend / 2) * L_dend
    "morphology, Dendritic properties, dendrite volume, um^3."
    Vol_dend::Float64 = pi * ((D_dend / 2)^2) * L_dend
    "morphology, Dendritic properties, dendritic membrane capacitance."
    Cdend::Float64 = Cm * A_dend
    "morphology, Dendritic properties, dendrite cross-sectional area, um^2."
    CS_dend::Float64 = pi * (D_dend / 2) .^ 2
    "morphology, Dendritic properties, nS."
    g_leakdend::Float64 = g_leak * A_dend
    "morphology, Soma properties, soma diameter, um."
    D_soma::Float64 = 30
    "morphology, Soma properties, soma surface area, 500 gives dendrite input resistance of 200MOhm, um^2."
    A_soma::Float64 = pi * D_soma^2
    "morphology, Soma properties, soma membrane capacitance."
    Csoma::Float64 = Cm * A_soma
    "morphology, Soma properties, soma cross-sectional area, um^2."
    CS_soma::Float64 = pi * (D_soma / 2) .^ 2
    "morphology, Soma properties, nS."
    g_leaksoma::Float64 = 15.0
    "morphology, Soma properties, value subject to modifications due to BaP adaptation implementation."
    g_diff::Float64 = D_dend / (4R_a)
    "spine properties, spine head volume [bartol2015], um^3."
    Vol_sp::Float64 = 0.03
    "spine properties, spine head surface area."
    A_sp::Float64 = 4 * pi * ((3 * Vol_sp) / (4 * pi))^(2.0 / 3.0)
    "spine properties, spine membrane capacitance."
    Csp::Float64 = Cm * A_sp
    "spine properties, spine head leak conductance."
    g_leaksp::Float64 = g_leak * A_sp
    "neck properties, spine neck diameter [bartol2015], um."
    D_neck::Float64 = 0.1
    "neck properties, neck length [bartol2015], um."
    L_neck::Float64 = 0.2
    "neck properties, neck cross sectional area, um^2."
    CS_neck::Float64 = pi * (D_neck / 2) .^ 2
    "neck properties."
    g_neck::Float64 = CS_neck / (L_neck * R_a)
    "neck properties."
    tau_diff::Float64 = ((Vol_sp / (2 * D_Ca * D_neck)) + (L_neck^2 / (2 * D_Ca)))
    "synpatic glutamate transient parameters, arbitrary, ms."
    glu_width::Float64 = 1.0 # 0.1 ms for synapse
    "synpatic glutamate transient parameters, arbitrary, mM."
    glu_amp::Float64 = 1e+3
    "synpatic glutamate transient parameters [liu1999]."
    glu_cv::Float64 = 0.5
    "SK channels, number of SK channels."
    N_SK::Int64 = 15
    "SK channels [maylie2004], ns."
    SK_gamma::Float64 = 10e-3
    "SK channels [mellor2016annex], mv."
    SK_Erev::Float64 = -90
    "SK channels [mellor2016annex], uM."
    SK_gating_half::Float64 = 0.33
    "SK channels [mellor2016annex], ms."
    SK_time::Float64 = 6.3
    "SK channels [mellor2016annex], ms."
    SK_hill::Float64 = 6
    "SK channels."
    p_SK_bcwd::NTuple{4,Float64} =
        (0.09391588258147192, 98.85165844770867, -147.61669527876904, 149.37767054612135)
    "SK channels."
    bcwd_SK::Float64 = (
        p_SK_bcwd[4] + p_SK_bcwd[3] / (1 + exp(p_SK_bcwd[1] * (temp_rates - p_SK_bcwd[2])))
    )
    "SK channels."
    p_SK_frwd::NTuple{4,Float64} =
        (-0.334167923607112, 25.590920461511878, 2.2052151559841193, 0.005904170174699533)
    "SK channels."
    frwd_SK::Float64 = (
        p_SK_frwd[4] + p_SK_frwd[3] / (1 + exp(p_SK_frwd[1] * (temp_rates - p_SK_frwd[2])))
    )
    "CaM, CaMKII and CaN concentrations."
    CaM_con::Float64 = 30.0
    "CaM, CaMKII and CaN concentrations, renamed [feng2011], 100um."
    mKCaM_con::Float64 = 70.0
    "CaM, CaMKII and CaN concentrations [lisman?], uM."
    mCaN_con::Float64 = 20.0
    "Chang Pepke model - CaM reactions I."
    kon_1C::Float64 = 5e-3
    "Chang Pepke model - CaM reactions I."
    koff_1C::Float64 = 50e-3
    "Chang Pepke model - CaM reactions I."
    kon_2C::Float64 = 10e-3
    "Chang Pepke model - CaM reactions I."
    koff_2C::Float64 = 10e-3
    "Chang Pepke model - CaM reactions I."
    kon_1N::Float64 = 100e-3
    "Chang Pepke model - CaM reactions I."
    koff_1N::Float64 = 2000e-3
    "Chang Pepke model - CaM reactions I."
    kon_2N::Float64 = 200e-3
    "Chang Pepke model - CaM reactions I."
    koff_2N::Float64 = 500e-3
    "Chang Pepke model - CaM reactions II."
    kf_CaM0::Float64 = 3.8e-6
    "Chang Pepke model - CaM reactions II."
    kb_CaM0::Float64 = 5.5e-3
    "Chang Pepke model - CaM reactions II."
    kf_CaM2C::Float64 = 0.92e-3
    "Chang Pepke model - CaM reactions II."
    kb_CaM2C::Float64 = 6.8e-3
    "Chang Pepke model - CaM reactions II."
    kf_CaM2N::Float64 = 0.12e-3
    "Chang Pepke model - CaM reactions II."
    kb_CaM2N::Float64 = 1.7e-3
    "Chang Pepke model - CaM reactions II."
    kf_CaM4::Float64 = 30e-3
    "Chang Pepke model - CaM reactions II."
    kb_CaM4::Float64 = 1.5e-3
    "Chang Pepke model - CaMKII reactions."
    kon_K1C::Float64 = 44e-3
    "Chang Pepke model - CaMKII reactions."
    koff_K1C::Float64 = 33e-3
    "Chang Pepke model - CaMKII reactions."
    kon_K2C::Float64 = 44e-3
    "Chang Pepke model - CaMKII reactions."
    koff_K2C::Float64 = 0.8e-3
    "Chang Pepke model - CaMKII reactions."
    kon_K1N::Float64 = 76e-3
    "Chang Pepke model - CaMKII reactions."
    koff_K1N::Float64 = 300e-3
    "Chang Pepke model - CaMKII reactions."
    kon_K2N::Float64 = 76e-3
    "Chang Pepke model - CaMKII reactions."
    koff_K2N::Float64 = 20e-3
    "Chang Pepke model - autophosphorilation."
    p_camkii_q10::NTuple{4,Float64} =
        (0.5118207068695309, 45.47503600542303, -161.42634157226917, 162.1718925882677)
    "Chang Pepke model - autophosphorilation."
    q10::Float64 =
        p_camkii_q10[4] +
        p_camkii_q10[3] / (1 + exp(p_camkii_q10[1] * (temp_rates - p_camkii_q10[2]))) # change of temp to fit chang 35C
    "Chang Pepke model - autophosphorilation."
    k1::Float64 = 12.6e-3
    "Chang Pepke model - autophosphorilation."
    k2::Float64 = q10 * 0.33e-3 # q10 * 0.33e-3
    "Chang Pepke model - autophosphorilation."
    k3::Float64 = 4 * q10 * 0.17e-3 # q10 * 0.17e-3
    "Chang Pepke model - autophosphorilation."
    k4::Float64 = 4 * 0.041e-3
    "Chang Pepke model - autophosphorilation."
    k5::Float64 = 4 * q10 * 2 * 0.017e-3 # q10 * 2* 0.017e-3
    "CaM-CaN reactions."
    p_CaN_frwd::NTuple{4,Float64} =
        (-0.29481489145354556, 29.999999999999968, 0.15940019940354327, 0.870299900298228)
    "CaM-CaN reactions, 22C - 4.6e-2 [quintana2005]."
    kcanf::Float64 =
        (
            p_CaN_frwd[4] +
            p_CaN_frwd[3] / (1 + exp(p_CaN_frwd[1] * (temp_rates - p_CaN_frwd[2])))
        ) * 1.75e-2
    "CaM-CaN reactions."
    p_CaN_bcwd::NTuple{4,Float64} =
        (-0.6833299932488973, 26.277500129849113, 0.7114524682690591, 0.29037766196937326)
    "CaM-CaN reactions, 22C - 1.2e-6 [quintana2005]."
    kcanb::Float64 =
        (
            p_CaN_bcwd[4] +
            p_CaN_bcwd[3] / (1 + exp(p_CaN_bcwd[1] * (temp_rates - p_CaN_bcwd[2])))
        ) * 2e-5
    "VGCCs."
    p_frwd_VGCC::NTuple{4,Float64} =
        (1.0485098341579628, 30.66869198447378, -0.3040010721391852, 2.5032059559264357)
    "VGCCs."
    frwd_VGCC::Float64 = (
        p_frwd_VGCC[4] +
        p_frwd_VGCC[3] / (1 + exp(p_frwd_VGCC[1] * (temp_rates - p_frwd_VGCC[2])))
    )
    "VGCCs."
    p_bcwd_VGCC::NTuple{4,Float64} =
        (-0.3302682317933842, 36.279019647221226, 3.2259761593440155, 0.7298285671937866)
    "VGCCs."
    bcwd_VGCC::Float64 = (
        p_bcwd_VGCC[4] +
        p_bcwd_VGCC[3] / (1 + exp(p_bcwd_VGCC[1] * (temp_rates - p_bcwd_VGCC[2])))
    )
    "VGCCs, calcium channel reversal potential, mV."
    Erev_CaT::Float64 = 10.0
    "VGCCs, calcium channel reversal potential, mV."
    Erev_CaR::Float64 = 10.0
    "VGCCs, calcium channel reversal potential, mV."
    Erev_CaL::Float64 = 10.0
    "VGCCs [magee1995], nS."
    gamma_CaT::Float64 = 12e-3
    "VGCCs [magee1995], nS."
    gamma_CaR::Float64 = 17e-3
    "VGCCs [magee1995], nS."
    gamma_CaL::Float64 = 27e-3
    "VGCCs."
    N_caT::Int64 = 3
    "VGCCs."
    N_caR::Int64 = 3
    "VGCCs."
    N_caL::Int64 = 3
    "calcium dye and buffers [zenisek2003,naraghi1997], uMms-1."
    EGTA_kf::Float64 = 2.7e-3
    "calcium dye and buffers, assuming Kd of 0.18uM [naraghi1997] ms-1."
    EGTA_kb::Float64 = 0.18 * EGTA_kf
    "calcium dye and buffers, 0.2 for imagin, 200 for elecrophysiology [tigaret2016] uM."
    EGTA_con::Float64 = 0.0
    "calcium dye and buffers [zenisek2003,naraghi1997], uM-1ms-1."
    BAPTA_kf::Float64 = 0.45
    "calcium dye and buffers, assuming Kd of 0.176uM [naraghi1997], ms-1."
    BAPTA_kb::Float64 = 0.176 * BAPTA_kf
    "calcium dye and buffers, uM."
    BAPTA_con::Float64 = 0.0
    "calcium dye and buffers [bartol2015], uM-1ms-1."
    Imbuf_k_on::Float64 = 0.247
    "calcium dye and buffers [bartol2015], ms-1."
    Imbuf_k_off::Float64 = 0.524
    "calcium dye and buffers."
    K_buff_diss::Float64 = Imbuf_k_off / Imbuf_k_on
    "calcium dye and buffers, 76.7 [bartol2015], uM."
    Imbuf_con::Float64 = 62
    "calcium dye and buffers."
    Imbuf_con_dend::Float64 = Imbuf_con * 4
    "calcium fluorescent dyes, assuming a [Ca] = 1um [bartol2015], ms-1."
    ogb1_kf::Float64 = 0.8
    "calcium fluorescent dyes [bartol2015], ms-1."
    ogb1_kb::Float64 = 0.16
    "calcium fluorescent dyes, assuming a [Ca] = 1um [bartol2015], ms-1."
    fluo4_kf::Float64 = 0.8
    "calcium fluorescent dyes [bartol2015], ms-1."
    fluo4_kb::Float64 = 0.24
    "calcium fluorescent dyes."
    dye::Float64 = 0.0
    "calcium fluorescent dyes [zenisek2003,naraghi1997], uMms-1."
    fluo5f_kf::Float64 = dye * 0.01
    "calcium fluorescent dyes assuming [Kd] = 1.3uM [yasuda2004]."
    fluo5f_kb::Float64 = dye * 26 * fluo5f_kf
    "calcium fluorescent dyes uM [tigaret2016]."
    fluo5f_con::Float64 = dye * 200.0
end
