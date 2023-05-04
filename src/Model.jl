# New synapse model (cian, started 23/03/2018) that replaces NMDA model with
# fully state-based one from Jahr and Stevens, plus three types of VGCCs
# (R-type, T-type and L-type), from Magee and Johhston (1995).

function F_synapse(dxc, xc, xd, p_synapse::SynapseParams, t, events_bap, bap_by_epsp)
    @unpack_SynapseParams p_synapse

    ##### Stochastic channels/receptors
    n1_ampa = xd[14] # ampa subconductance 1
    n2_ampa = xd[15] # ampa subconductance 2
    n3_ampa = xd[16] # ampa subconductance 3
    n1_nmda_A = xd[22] # nmda subconductance 1
    n2_nmda_A = xd[23] # nmda subconductance 2
    n1_nmda_B = xd[44] # nmda subconductance 1
    n2_nmda_B = xd[45] # nmda subconductance 2
    n_car = xd[28] # vgcc-R opened state
    n_cat = xd[32] # vgcc-T opened state
    n_cal = xd[34] + xd[35] # vgcc-L opened states
    n_gaba1 = xd[49] # GABA opened state
    n_gaba2 = xd[50] # GABA opened state

    ##### Continuous variables
    Vsp,
    Vdend,
    Vsoma,
    λ,
    ImbufCa,
    Ca,
    Dye,
    CaM0,
    CaM2C,
    CaM2N,
    CaM4,
    mCaN,
    CaN4,
    mKCaM,
    KCaM0,
    KCaM2N,
    KCaM2C,
    KCaM4,
    PCaM0,
    PCaM2C,
    PCaM2N,
    PCaM4,
    P,
    P2,
    LTD,
    LTP,
    act_D,
    act_P,
    m,
    h,
    n,
    SK,
    λ_age,
    λ_aux = xc

    ##### plasticity prediction regions
    CaMKII = KCaM0 + KCaM2C + KCaM2N + KCaM4 + PCaM0 + PCaM2C + PCaM2N + PCaM4 + P + P2
    CaN = CaN4

    #### activation when it is inside the region
    # this following line allocates 74.779 ns (3 allocations: 144 bytes). The 2 lines count for 25% of the performance
    ∂LTD = SVector(CaN, CaMKII) ∈ LTD_region
    ∂LTP = SVector(CaN, CaMKII) ∈ LTP_region

    ∂act_D = a_D * ∂LTD - b_D * act_D * (1 - ∂LTD)
    ∂act_P = a_P * ∂LTP - b_P * act_P * (1 - ∂LTP)

    ##### Na channel
    m_inf = alpha_m(Vsoma) / (alpha_m(Vsoma) + beta_m(Vsoma))
    m_tau = 1 / (alpha_m(Vsoma) + beta_m(Vsoma))
    ∂m = (m_inf - m) / m_tau
    ∂h = alpha_h(Vsoma) * (1 - h) - beta_h(Vsoma) * h
    I_Na = gamma_Na * (m^3) * h * (Erev_Na - Vsoma)


    ##### K channel
    n_inf = 1 / (1 + alpha_n(Vsoma))
    n_tau = max(50 * beta_n(Vsoma) / (1 + alpha_n(Vsoma)), 2.0)
    ∂n = (n_inf - n) / n_tau
    I_K = gamma_K * n * (Erev_K - Vsoma)

    ##### NMDA
    NMDA = (n1_nmda_A + n2_nmda_A + n1_nmda_B + n2_nmda_B) * B(Vsp, Mg) * gamma_nmda
    Inmda = (Erev_nmda - Vsp) * NMDA # current nmda

    ##### AMPA
    Iampa =
        (Erev_ampa - Vsp) *
        (gamma_ampa1 * n1_ampa + gamma_ampa2 * n2_ampa + gamma_ampa3 * n3_ampa) # current ampa

    ##### GABA
    Igaba = (n_gaba1 + n_gaba2) * (Erev_Cl - Vdend) * gamma_GABA

    ##### Calcium sources (VGCCs currents, and NMDA calcium contribution)
    ΦCa = perm * ghk(Vsp, Ca, Ca_ext, p_synapse) #GHK factor
    Ica_nmda = f_Ca * ΦCa * NMDA
    Icar = gamma_CaR * n_car * ΦCa
    Icat = gamma_CaT * n_cat * ΦCa
    Ical = gamma_CaL * n_cal * ΦCa

    ##### SK channel (not stochastic)
    ∂SK = (SK_chnnl(Ca) * frwd_SK - SK) / (SK_time * bcwd_SK) #SK spine
    Isk = SK_gamma * (SK_Erev - Vsp) * SK * N_SK

    ##### Backpropgation
    # Post input - for experimentally induced BaPs and those induced by EPSPs
    I_BaP =
        inputBaP(t, bap_by_epsp, injbap, I_clamp) + inputBaP(t, events_bap, injbap, I_clamp)
    # Bap decay/attenuation - two component for adaptation in the Bap
    ∂λ = (1 - λ) / trec - delta_decay * (1 / λ_aux) * λ * I_BaP
    ∂λ_aux = (1 - λ_aux) / trec - delta_aux * λ_aux * I_BaP
    gadapt = λ * g_diff * ϕ_dist

    # Bap decay/attenuation - age dependent modification factor
    ∂λ_age = (1 - λ_age) / trec_soma - delta_soma * λ_age * I_BaP

    ##### Voltage
    # Spine
    ∂Vsp =
        (
            Isk +
            Inmda +
            Iampa +
            Icat +
            Icar +
            Ical +
            g_neck * (Vdend - Vsp) +
            g_leak * (E_leak - Vsp)
        ) / (Csp)
    # Dendrite
    ∂Vdend =
        (
            g_neck * (Vsp - Vdend) +
            Igaba +
            g_leakdend * (E_leak - Vdend) +
            gadapt * (Vsoma - Vdend)
        ) / Cdend
    # Soma
    ∂Vsoma =
        (
            (I_BaP + I_Na) * λ_age +
            I_K +
            g_leaksoma * (E_leak - Vsoma) +
            gadapt * (Vdend - Vsoma)
        ) / Csoma


    ##### Buffer and dye (spine only - no neck diffusion)
    ∂ImbufCa = Imbuf_k_on * (Imbuf_con - ImbufCa) * Ca - Imbuf_k_off * ImbufCa
    ∂Dye = 4 * fluo5f_kf * (fluo5f_con - Dye) * Ca - 8 * fluo5f_kb * Dye

    ##### Ca Downstream
    ### CaM-KCaM-rates (coarsed model) from Pepke adapted by
    kf_2C = rates_adapt(kon_1C, kon_2C, koff_1C, kon_2C, Ca)
    kb_2C = rates_adapt(koff_1C, koff_2C, koff_1C, kon_2C, Ca)
    kf_2N = rates_adapt(kon_1N, kon_2N, koff_1N, kon_2N, Ca)
    kb_2N = rates_adapt(koff_1N, koff_2N, koff_1N, kon_2N, Ca)
    kf_K2C = rates_adapt(kon_K1C, kon_K2C, koff_K1C, kon_K2C, Ca)
    kb_K2C = rates_adapt(koff_K1C, koff_K2C, koff_K1C, kon_K2C, Ca)
    kf_K2N = rates_adapt(kon_K1N, kon_K2N, koff_K1N, kon_K2N, Ca)
    kb_K2N = rates_adapt(koff_K1N, koff_K2N, koff_K1N, kon_K2N, Ca)
    F = CaMKII / mKCaM_con

    ∂CaM0 =
        k2 * PCaM0 +
        kb_2C * CaM2C +
        kb_2N * CaM2N +
        kb_CaM0 * KCaM0 +
        -(1 // 2) * kf_2C * (Ca^2) * CaM0 - (1 // 2) * kf_2N * (Ca^2) * CaM0 +
        -kf_CaM0 * CaM0 * mKCaM

    ∂CaM2C =
        kb_2N * CaM4 + kb_CaM2C * KCaM2C + k2 * PCaM2C + +(1 // 2) * kf_2C * (Ca^2) * CaM0 -
        kb_2C * CaM2C - (1 // 2) * kf_2N * (Ca^2) * CaM2C + -kf_CaM2C * CaM2C * mKCaM

    ∂CaM2N =
        kb_2C * CaM4 + kb_CaM2N * KCaM2N + k2 * PCaM2N + +(1 // 2) * kf_2N * (Ca^2) * CaM0 -
        kb_2N * CaM2N - (1 // 2) * kf_2C * (Ca^2) * CaM2N + -kf_CaM2N * CaM2N * mKCaM

    ∂CaM4 =
        k2 * PCaM4 +
        kcanb * CaN4 +
        kb_CaM4 * KCaM4 +
        +(1 // 2) * kf_2N * (Ca^2) * CaM2C +
        (1 // 2) * kf_2C * (Ca^2) * CaM2N - kb_2C * CaM4 + -kb_2N * CaM4 -
        kcanf * CaM4 * mCaN - kf_CaM4 * CaM4 * mKCaM

    ∂mCaN = kcanb * CaN4 - kcanf * CaM4 * mCaN

    ∂CaN4 = kcanf * CaM4 * mCaN - kcanb * CaN4

    ∂mKCaM =
        kb_CaM0 * KCaM0 +
        k3 * P +
        kb_CaM2C * KCaM2C +
        kb_CaM2N * KCaM2N +
        +kb_CaM4 * KCaM4 - kf_CaM0 * CaM0 * mKCaM - kf_CaM2C * CaM2C * mKCaM +
        -kf_CaM2N * CaM2N * mKCaM - kf_CaM4 * CaM4 * mKCaM

    ∂KCaM0 =
        kb_K2C * KCaM2C + kb_K2N * KCaM2N + kf_CaM0 * CaM0 * mKCaM + -kb_CaM0 * KCaM0 -
        (1 // 2) * kf_K2C * (Ca^2) * KCaM0 - F * k1 * KCaM0 +
        -(1 // 2) * kf_K2N * (Ca^2) * KCaM0

    ∂KCaM2N =
        kb_K2C * KCaM4 + kf_CaM2N * CaM2N * mKCaM + +(1 // 2) * kf_K2N * (Ca^2) * KCaM0 -
        kb_CaM2N * KCaM2N - kb_K2N * KCaM2N + -(1 // 2) * kf_K2C * (Ca^2) * KCaM2N -
        F * k1 * KCaM2N

    ∂KCaM2C =
        kb_K2N * KCaM4 + kf_CaM2C * CaM2C * mKCaM + +(1 // 2) * kf_K2C * (Ca^2) * KCaM0 -
        kb_CaM2C * KCaM2C - kb_K2C * KCaM2C + -F * k1 * KCaM2C -
        (1 // 2) * kf_K2N * (Ca^2) * KCaM2C

    ∂KCaM4 =
        kf_CaM4 * CaM4 * mKCaM +
        (1 // 2) * kf_K2C * (Ca^2) * KCaM2N +
        +(1 // 2) * kf_K2N * (Ca^2) * KCaM2C - kb_CaM4 * KCaM4 - kb_K2C * KCaM4 +
        -kb_K2N * KCaM4 - F * k1 * KCaM4

    ∂PCaM0 = F * k1 * KCaM0 - k2 * PCaM0

    ∂PCaM2N = F * k1 * KCaM2N - k2 * PCaM2N

    ∂PCaM2C = F * k1 * KCaM2C - k2 * PCaM2C

    ∂PCaM4 = F * k1 * KCaM4 - k2 * PCaM4

    ∂P = k2 * PCaM0 + k5 * P2 + k2 * PCaM2C + k2 * PCaM2N + k2 * PCaM4 - k3 * P - k4 * P

    ∂P2 = k4 * P - k5 * P2


    ### Postsynaptic Ca
    ∂Ca =
        (Ca_infty - Ca) / tau_ca +
        +(Ica_nmda + Icar + Ical + Icat) / (2 * faraday * A_sp) +
        +(max(Ca_infty, Ca / 3) - Ca) / tau_diff +
        -∂ImbufCa +
        -∂Dye +
        +2kb_2C * CaM2C +
        2kb_2C * CaM4 +
        2kb_2N * CaM2N +
        2kb_2N * CaM4 +
        +2kb_K2C * KCaM2C +
        2kb_K2N * KCaM2N +
        2kb_K2C * KCaM4 +
        2kb_K2N * KCaM4 +
        -kf_2C * (Ca^2) * CaM0 - kf_2N * (Ca^2) * CaM0 - kf_2N * (Ca^2) * CaM2C +
        -kf_2C * (Ca^2) * CaM2N - kf_K2C * (Ca^2) * KCaM0 - kf_K2C * (Ca^2) * KCaM2N +
        -kf_K2N * (Ca^2) * KCaM0 - kf_K2N * (Ca^2) * KCaM2C

    ### dxc update
    dxc[1] = ∂Vsp
    dxc[2] = ∂Vdend
    dxc[3] = ∂Vsoma
    dxc[4] = ∂λ
    dxc[5] = ∂ImbufCa
    dxc[6] = ∂Ca
    dxc[7] = ∂Dye
    dxc[8] = ∂CaM0
    dxc[9] = ∂CaM2C
    dxc[10] = ∂CaM2N
    dxc[11] = ∂CaM4
    dxc[12] = ∂mCaN
    dxc[13] = ∂CaN4
    dxc[14] = ∂mKCaM
    dxc[15] = ∂KCaM0
    dxc[16] = ∂KCaM2N
    dxc[17] = ∂KCaM2C
    dxc[18] = ∂KCaM4
    dxc[19] = ∂PCaM0
    dxc[20] = ∂PCaM2C
    dxc[21] = ∂PCaM2N
    dxc[22] = ∂PCaM4
    dxc[23] = ∂P
    dxc[24] = ∂P2
    dxc[25] = ∂LTD
    dxc[26] = ∂LTP
    dxc[27] = ∂act_D
    dxc[28] = ∂act_P
    dxc[29] = ∂m
    dxc[30] = ∂h
    dxc[31] = ∂n
    dxc[32] = ∂SK
    dxc[33] = ∂λ_age
    dxc[34] = ∂λ_aux

end

function R_synapse(rate, xc, xd, p_synapse::SynapseParams, t, sum_rate, glu = 0)

    @unpack_SynapseParams p_synapse

    ############### Voltage ###################
    Vsp = xc[1]

    ############### Glutamate & GABA ###################
    Glu = glu_amp * glu

    ############### AMPA ###################
    #2line-GO
    rate[1] = 4 * AMPA_k1 * Glu * xd[1]
    rate[2] = 3 * AMPA_k1 * Glu * xd[2]
    rate[3] = 2 * AMPA_k1 * Glu * xd[3]
    rate[4] = 1 * AMPA_k1 * Glu * xd[4]
    #2line-BACK
    rate[5] = 4 * AMPA_k_1 * xd[5]
    rate[6] = 3 * AMPA_k_1 * xd[4]
    rate[7] = 2 * AMPA_k_1 * xd[3]
    rate[8] = 1 * AMPA_k_1 * xd[2]
    #3line-GO
    rate[9] = 3 * AMPA_k1 * Glu * xd[6]
    rate[10] = 3 * AMPA_k1 * Glu * xd[7]
    rate[11] = 2 * AMPA_k1 * Glu * xd[8]
    rate[12] = 1 * AMPA_k1 * Glu * xd[9]
    #3line-BACK
    rate[13] = 3 * AMPA_k_1 * xd[10]
    rate[14] = 2 * AMPA_k_1 * xd[9]
    rate[15] = 1 * AMPA_k_1 * xd[8]
    rate[16] = 1 * AMPA_k_2 * xd[7]
    #4line-GO
    rate[17] = 2 * AMPA_k1 * Glu * xd[11]
    rate[18] = 1 * AMPA_k1 * Glu * xd[12]
    #4line-BACK
    rate[19] = 2 * AMPA_k_1 * xd[13]
    rate[20] = 1 * AMPA_k_1 * xd[12]
    #1column-GO-BACK
    rate[21] = 4 * AMPA_delta_0 * xd[1]
    rate[22] = 1 * AMPA_gamma_0 * xd[6]
    #2column-GO-BACK
    rate[23] = 1 * AMPA_delta_1 * xd[2]
    rate[24] = 1 * AMPA_gamma_1 * xd[7]
    #3column-GO
    rate[25] = 1 * AMPA_alpha * xd[14]
    rate[26] = 2 * AMPA_delta_1 * xd[3]
    rate[27] = 1 * AMPA_delta_2 * xd[8]
    #3column-BACK
    rate[28] = 1 * AMPA_gamma_2 * xd[11]
    rate[29] = 1 * AMPA_gamma_1 * xd[8]
    rate[30] = 2 * AMPA_beta * xd[3]
    #4column-GO
    rate[31] = 1 * AMPA_alpha * xd[15]
    rate[32] = 3 * AMPA_delta_1 * xd[4]
    rate[33] = 2 * AMPA_delta_2 * xd[9]
    #4column-BACK
    rate[34] = 1 * AMPA_gamma_2 * xd[12]
    rate[35] = 1 * AMPA_gamma_1 * xd[9]
    rate[36] = 2 * AMPA_beta * xd[4]
    #5column-GO
    rate[37] = 1 * AMPA_alpha * xd[16]
    rate[38] = 4 * AMPA_delta_1 * xd[5]
    rate[39] = 3 * AMPA_delta_2 * xd[10]
    #5column-BACK
    rate[40] = 1 * AMPA_gamma_2 * xd[13]
    rate[41] = 1 * AMPA_gamma_1 * xd[10]
    rate[42] = 4 * AMPA_beta * xd[5]

    ############### NMDA ###################
    #1line-GO
    rate[43] = NMDA_N2A_ka * xd[17] * Glu
    rate[44] = NMDA_N2A_kb * xd[18] * Glu
    rate[45] = NMDA_N2A_kc * xd[19]
    rate[46] = NMDA_N2A_kd * xd[20]
    rate[47] = NMDA_N2A_ke * xd[21]
    rate[48] = NMDA_N2A_kf * xd[22]
    #1line-BACK
    rate[49] = NMDA_N2A_k_f * xd[23]
    rate[50] = NMDA_N2A_k_e * xd[22]
    rate[51] = NMDA_N2A_k_d * xd[21]
    rate[52] = NMDA_N2A_k_c * xd[20]
    rate[53] = NMDA_N2A_k_b * xd[19]
    rate[54] = NMDA_N2A_k_a * xd[18]

    ################### Sampling ###################
    rate[55] = sampling_rate

    ################### R-type VGCC ###################
    alpha_m_r, beta_m_r = rates_m_r(Vsp)
    alpha_h_r, beta_h_r = rates_h_r(Vsp)

    rate[56] = xd[25] * alpha_m_r * frwd_VGCC
    rate[57] = xd[26] * beta_m_r * bcwd_VGCC
    rate[58] = xd[25] * alpha_h_r * frwd_VGCC
    rate[59] = xd[27] * beta_h_r * bcwd_VGCC
    rate[60] = xd[26] * alpha_h_r * frwd_VGCC
    rate[61] = xd[28] * beta_h_r * bcwd_VGCC
    rate[62] = xd[27] * alpha_m_r * frwd_VGCC
    rate[63] = xd[28] * beta_m_r * bcwd_VGCC

    ################### T-type VGCC  ###################
    alpha_m_t, beta_m_t = rates_m_t(Vsp)
    alpha_h_t, beta_h_t = rates_h_t(Vsp)

    rate[64] = xd[29] * alpha_m_t * frwd_VGCC
    rate[65] = xd[30] * beta_m_t * bcwd_VGCC # this one can have a high rate
    rate[66] = xd[29] * alpha_h_t * frwd_VGCC
    rate[67] = xd[31] * beta_h_t * bcwd_VGCC
    rate[68] = xd[30] * alpha_h_t * frwd_VGCC
    rate[69] = xd[32] * beta_h_t * bcwd_VGCC
    rate[70] = xd[31] * alpha_m_t * frwd_VGCC
    rate[71] = xd[32] * beta_m_t * bcwd_VGCC # this one can have a high rate

    ################### L-type VGCC  ###################
    alpha_l, beta_1_l, beta_2_l = rates_l(Vsp)
    rate[72] = xd[33] * alpha_l * frwd_VGCC
    rate[73] = xd[34] * beta_1_l * bcwd_VGCC
    rate[74] = xd[33] * alpha_l * frwd_VGCC
    rate[75] = xd[35] * beta_2_l * bcwd_VGCC

    ################### LTD/LTP  ###################
    #the 6 lines take 50ns on 200ns, 1/4 of computations are here!!
    D_rate = plasticityRate(xc[27], 2, K_D) / t_D
    P_rate = plasticityRate(xc[28], 2, K_P) / t_P
    rate[76] = xd[36] * D_rate
    rate[77] = xd[37] * P_rate
    rate[78] = xd[36] * P_rate
    rate[79] = xd[38] * D_rate

    ############### NMDA GLUN2B ###################
    #1line-GO
    rate[80] = NMDA_N2B_sa * xd[39] * Glu
    rate[81] = NMDA_N2B_sb * xd[40] * Glu
    rate[82] = NMDA_N2B_sc * xd[41]
    rate[83] = NMDA_N2B_sd * xd[42]
    rate[84] = NMDA_N2B_se * xd[43]
    rate[85] = NMDA_N2B_sf * xd[44]

    #1line-BACK
    rate[86] = NMDA_N2B_s_f * xd[45]
    rate[87] = NMDA_N2B_s_e * xd[44]
    rate[88] = NMDA_N2B_s_d * xd[43]
    rate[89] = NMDA_N2B_s_c * xd[42]
    rate[90] = NMDA_N2B_s_b * xd[41]
    rate[91] = NMDA_N2B_s_a * xd[40]

    ############### GABA ###################
    rate[92] = GABA_r_b1 * xd[46] * Glu #to simplify, we use the same ammount at the same time
    rate[93] = GABA_r_u1 * xd[47]
    rate[94] = GABA_r_b2 * xd[47] * Glu
    rate[95] = GABA_r_u2 * xd[48]
    rate[96] = GABA_r_ro1 * xd[47]
    rate[97] = GABA_r_c1 * xd[49]
    rate[98] = GABA_r_ro2 * xd[48]
    rate[99] = GABA_r_c2 * xd[50]

    bound = 0.0
    if sum_rate == false
        return 0.0, bound
    else
        return sum(rate), bound
    end
end

macro j_jump(i, p_synapse, nu, rate_ex, urate_ex = nothing, rateinterval_ex = nothing)

    assignments = Expr[]

    alpha_beta_regex = r"(alpha|beta)_(m_r|h_r|m_t|h_t|l|1_l|2_l)"
    alpha_beta_matches = Set([m.match for m in eachmatch(alpha_beta_regex, "$rate_ex")])

    if length(alpha_beta_matches) > 0

        for m in ("alpha_1_l", "alpha_2_l", "beta_l")
            if m in alpha_beta_matches
                throw(DomainError(m, "this variable does not exist in the model."))
            end
        end

        push!(assignments, :(Vsp = u[1]))

        if "alpha_m_r" in alpha_beta_matches || "beta_m_r" in alpha_beta_matches
            push!(assignments, :((alpha_m_r, beta_m_r) = rates_m_r(Vsp)))
        end

        if "alpha_h_r" in alpha_beta_matches || "beta_h_r" in alpha_beta_matches
            push!(assignments, :((alpha_h_r, beta_h_r) = rates_h_r(Vsp)))
        end

        if "alpha_m_t" in alpha_beta_matches || "beta_m_t" in alpha_beta_matches
            push!(assignments, :((alpha_m_t, beta_m_t) = rates_m_t(Vsp)))
        end

        if "alpha_h_t" in alpha_beta_matches || "beta_h_t" in alpha_beta_matches
            push!(assignments, :((alpha_h_t, beta_h_t) = rates_h_t(Vsp)))
        end

        if "alpha_l" in alpha_beta_matches ||
           "beta_1_l" in alpha_beta_matches ||
           "beta_2_l" in alpha_beta_matches
            push!(assignments, :((alpha_l, beta_1_l, beta_2_l) = rates_l(Vsp)))
        end

    end

    if occursin("D_rate", "$rate_ex")
        push!(assignments, :(D_rate = plasticityRate(u[27], 2, K_D) / t_D))
    end

    if occursin("P_rate", "$rate_ex")
        push!(assignments, :(P_rate = plasticityRate(u[28], 2, K_D) / t_P))
    end

    ex = Expr[]

    push!(
        ex,
        quote
            @unpack_SynapseParams $(esc(p_synapse))
            @inline @inbounds function rate(u, p, t)
                $(assignments...)
                return $rate_ex
            end
            @inline @inbounds function affect!(integrator)
                for (j, a) in zip(findnz($(esc(nu))[$(esc(i)), :])...)
                    integrator.p.xd[j] += a
                end
            end
        end,
    )

    if !isnothing(urate_ex)
        push!(ex, quote
            max_m_r = rates_m_r(1_000.0)[1]
            max_h_r = rates_h_r(1_000.0)[2]
            max_m_t = rates_m_t(1_000.0)[1]
            max_h_t = rates_h_t(1_000.0)[2]
            max_alpha_l = rates_l(1_000)[1]
            max_beta_1_l, max_beta_2_l = rates_l(-1_000)[2:3]
            @inline @inbounds function urate(u, p, t)
                return $urate_ex
            end
            @inline @inbounds function rateinterval(u, p, t)
                return $rateinterval_ex
            end
        end)
    end

    if isnothing(urate_ex)
        push!(ex, :(ConstantRateJump(rate, affect!)))
    else
        push!(
            ex,
            :(VariableRateJump(rate, affect!; urate = urate, rateinterval = rateinterval)),
        )
    end

    quote
        $(ex...)
    end

end

function J_synapse(p_synapse::SynapseParams, nu)

    # we order the jumps in ther order they appear in the dependency graph
    jumps = JumpSet(;
        constant_jumps = [
            ############### AMPA ###################
            #2line-GO
            @j_jump(1, p_synapse, nu, 4 * AMPA_k1 * p.Glu * p.xd[1]), # 1
            @j_jump(2, p_synapse, nu, 3 * AMPA_k1 * p.Glu * p.xd[2]), # 2
            @j_jump(3, p_synapse, nu, 2 * AMPA_k1 * p.Glu * p.xd[3]), # 3
            @j_jump(4, p_synapse, nu, 1 * AMPA_k1 * p.Glu * p.xd[4]), # 4
            #2line-BACK
            @j_jump(5, p_synapse, nu, 4 * AMPA_k_1 * p.xd[5]), # 5
            @j_jump(6, p_synapse, nu, 3 * AMPA_k_1 * p.xd[4]), # 6
            @j_jump(7, p_synapse, nu, 2 * AMPA_k_1 * p.xd[3]), # 7
            @j_jump(8, p_synapse, nu, 1 * AMPA_k_1 * p.xd[2]), # 8
            #3line-GO
            @j_jump(9, p_synapse, nu, 3 * AMPA_k1 * p.Glu * p.xd[6]), # 9
            @j_jump(10, p_synapse, nu, 3 * AMPA_k1 * p.Glu * p.xd[7]), # 10
            @j_jump(11, p_synapse, nu, 2 * AMPA_k1 * p.Glu * p.xd[8]), # 11
            @j_jump(12, p_synapse, nu, 1 * AMPA_k1 * p.Glu * p.xd[9]), # 12
            #3line-BACK
            @j_jump(13, p_synapse, nu, 3 * AMPA_k_1 * p.xd[10]), # 13
            @j_jump(14, p_synapse, nu, 2 * AMPA_k_1 * p.xd[9]), # 14
            @j_jump(15, p_synapse, nu, 1 * AMPA_k_1 * p.xd[8]), # 15
            @j_jump(16, p_synapse, nu, 1 * AMPA_k_2 * p.xd[7]), # 16
            #4line-GO
            @j_jump(17, p_synapse, nu, 2 * AMPA_k1 * p.Glu * p.xd[11]), # 17
            @j_jump(18, p_synapse, nu, 1 * AMPA_k1 * p.Glu * p.xd[12]), # 18
            #4line-BACK
            @j_jump(19, p_synapse, nu, 2 * AMPA_k_1 * p.xd[13]), # 19
            @j_jump(20, p_synapse, nu, 1 * AMPA_k_1 * p.xd[12]), # 20
            #1column-GO-BACK
            @j_jump(21, p_synapse, nu, 4 * AMPA_delta_0 * p.xd[1]), # 21
            @j_jump(22, p_synapse, nu, 1 * AMPA_gamma_0 * p.xd[6]), # 22
            #2column-GO-BACK
            @j_jump(23, p_synapse, nu, 1 * AMPA_delta_1 * p.xd[2]), # 23
            @j_jump(24, p_synapse, nu, 1 * AMPA_gamma_1 * p.xd[7]), # 24
            #3column-GO
            @j_jump(25, p_synapse, nu, 1 * AMPA_alpha * p.xd[14]), # 25
            @j_jump(26, p_synapse, nu, 2 * AMPA_delta_1 * p.xd[3]), # 26
            @j_jump(27, p_synapse, nu, 1 * AMPA_delta_2 * p.xd[8]), # 27
            #3column-BACK
            @j_jump(28, p_synapse, nu, 1 * AMPA_gamma_2 * p.xd[11]), # 28
            @j_jump(29, p_synapse, nu, 1 * AMPA_gamma_1 * p.xd[8]), # 29
            @j_jump(30, p_synapse, nu, 2 * AMPA_beta * p.xd[3]), # 30
            #4column-GO
            @j_jump(31, p_synapse, nu, 1 * AMPA_alpha * p.xd[15]), # 31
            @j_jump(32, p_synapse, nu, 3 * AMPA_delta_1 * p.xd[4]), # 32
            @j_jump(33, p_synapse, nu, 2 * AMPA_delta_2 * p.xd[9]), # 33
            #4column-BACK
            @j_jump(34, p_synapse, nu, 1 * AMPA_gamma_2 * p.xd[12]), # 34
            @j_jump(35, p_synapse, nu, 1 * AMPA_gamma_1 * p.xd[9]), # 35
            @j_jump(36, p_synapse, nu, 2 * AMPA_beta * p.xd[4]), # 36
            #5column-GO
            @j_jump(37, p_synapse, nu, 1 * AMPA_alpha * p.xd[16]), # 37
            @j_jump(38, p_synapse, nu, 4 * AMPA_delta_1 * p.xd[5]), # 38
            @j_jump(39, p_synapse, nu, 3 * AMPA_delta_2 * p.xd[10]), # 39
            #5column-BACK
            @j_jump(40, p_synapse, nu, 1 * AMPA_gamma_2 * p.xd[13]), # 40
            @j_jump(41, p_synapse, nu, 1 * AMPA_gamma_1 * p.xd[10]), # 41
            @j_jump(42, p_synapse, nu, 4 * AMPA_beta * p.xd[5]), # 42

            ############### NMDA ###################
            #1line-GO
            @j_jump(43, p_synapse, nu, NMDA_N2A_ka * p.xd[17] * p.Glu), # 43
            @j_jump(44, p_synapse, nu, NMDA_N2A_kb * p.xd[18] * p.Glu), # 44
            @j_jump(45, p_synapse, nu, NMDA_N2A_kc * p.xd[19]), # 45
            @j_jump(46, p_synapse, nu, NMDA_N2A_kd * p.xd[20]), # 46
            @j_jump(47, p_synapse, nu, NMDA_N2A_ke * p.xd[21]), # 47
            @j_jump(48, p_synapse, nu, NMDA_N2A_kf * p.xd[22]), # 48
            #1line-BACK
            @j_jump(49, p_synapse, nu, NMDA_N2A_k_f * p.xd[23]), # 49
            @j_jump(50, p_synapse, nu, NMDA_N2A_k_e * p.xd[22]), # 50
            @j_jump(51, p_synapse, nu, NMDA_N2A_k_d * p.xd[21]), # 51
            @j_jump(52, p_synapse, nu, NMDA_N2A_k_c * p.xd[20]), # 52
            @j_jump(53, p_synapse, nu, NMDA_N2A_k_b * p.xd[19]), # 53
            @j_jump(54, p_synapse, nu, NMDA_N2A_k_a * p.xd[18]), # 54

            ############### NMDA GLUN2B ###################
            #1line-GO
            @j_jump(80, p_synapse, nu, NMDA_N2B_sa * p.xd[39] * p.Glu), # 80
            @j_jump(81, p_synapse, nu, NMDA_N2B_sb * p.xd[40] * p.Glu), # 81
            @j_jump(82, p_synapse, nu, NMDA_N2B_sc * p.xd[41]), # 82
            @j_jump(83, p_synapse, nu, NMDA_N2B_sd * p.xd[42]), # 83
            @j_jump(84, p_synapse, nu, NMDA_N2B_se * p.xd[43]), # 84
            @j_jump(85, p_synapse, nu, NMDA_N2B_sf * p.xd[44]), # 85

            #1line-BACK
            @j_jump(86, p_synapse, nu, NMDA_N2B_s_f * p.xd[45]), # 86
            @j_jump(87, p_synapse, nu, NMDA_N2B_s_e * p.xd[44]), # 87
            @j_jump(88, p_synapse, nu, NMDA_N2B_s_d * p.xd[43]), # 88
            @j_jump(89, p_synapse, nu, NMDA_N2B_s_c * p.xd[42]), # 89
            @j_jump(90, p_synapse, nu, NMDA_N2B_s_b * p.xd[41]), # 90
            @j_jump(91, p_synapse, nu, NMDA_N2B_s_a * p.xd[40]), # 91

            ############### GABA ###################
            @j_jump(92, p_synapse, nu, GABA_r_b1 * p.xd[46] * p.Glu), # 92 to simplify, we use the same ammount at the same time)
            @j_jump(93, p_synapse, nu, GABA_r_u1 * p.xd[47]), # 93
            @j_jump(94, p_synapse, nu, GABA_r_b2 * p.xd[47] * p.Glu), # 94
            @j_jump(95, p_synapse, nu, GABA_r_u2 * p.xd[48]), # 95
            @j_jump(96, p_synapse, nu, GABA_r_ro1 * p.xd[47]), # 96
            @j_jump(97, p_synapse, nu, GABA_r_c1 * p.xd[49]), # 97
            @j_jump(98, p_synapse, nu, GABA_r_ro2 * p.xd[48]), # 98
            @j_jump(99, p_synapse, nu, GABA_r_c2 * p.xd[50]), # 99
        ],
        variable_jumps = [
            ################### R-type VGCC ###################
            @j_jump(
                56,
                p_synapse,
                nu,
                p.xd[25] * alpha_m_r * frwd_VGCC,
                p.xd[25] * max_m_r * frwd_VGCC,
                typemax(Float64)
            ), # 56
            @j_jump(
                57,
                p_synapse,
                nu,
                p.xd[26] * beta_m_r * bcwd_VGCC,
                p.xd[26] * max_m_r * frwd_VGCC,
                typemax(Float64)
            ), # 57
            @j_jump(
                58,
                p_synapse,
                nu,
                p.xd[25] * alpha_h_r * frwd_VGCC,
                p.xd[25] * max_h_r * frwd_VGCC,
                typemax(Float64)
            ), # 58
            @j_jump(
                59,
                p_synapse,
                nu,
                p.xd[27] * beta_h_r * bcwd_VGCC,
                p.xd[27] * max_h_r * bcwd_VGCC,
                typemax(Float64)
            ), # 59
            @j_jump(
                60,
                p_synapse,
                nu,
                p.xd[26] * alpha_h_r * frwd_VGCC,
                p.xd[26] * max_h_r * frwd_VGCC,
                typemax(Float64)
            ), # 60
            @j_jump(
                61,
                p_synapse,
                nu,
                p.xd[28] * beta_h_r * bcwd_VGCC,
                p.xd[28] * max_h_r * bcwd_VGCC,
                typemax(Float64)
            ), # 61
            @j_jump(
                62,
                p_synapse,
                nu,
                p.xd[27] * alpha_m_r * frwd_VGCC,
                p.xd[27] * max_m_r * frwd_VGCC,
                typemax(Float64)
            ), # 62
            @j_jump(
                63,
                p_synapse,
                nu,
                p.xd[28] * beta_m_r * bcwd_VGCC,
                p.xd[28] * max_m_r * bcwd_VGCC,
                typemax(Float64)
            ), # 63

            ################### T-type VGCC  ###################
            @j_jump(
                64,
                p_synapse,
                nu,
                p.xd[29] * alpha_m_t * frwd_VGCC,
                p.xd[29] * max_m_t * frwd_VGCC,
                typemax(Float64)
            ), # 64
            @j_jump(
                65,
                p_synapse,
                nu,
                p.xd[30] * beta_m_t * bcwd_VGCC,
                p.xd[30] * max_m_t * bcwd_VGCC,
                typemax(Float64)
            ), # 65 this one can have a high rate
            @j_jump(
                66,
                p_synapse,
                nu,
                p.xd[29] * alpha_h_t * frwd_VGCC,
                p.xd[29] * max_h_t * frwd_VGCC,
                typemax(Float64)
            ), # 66
            @j_jump(
                67,
                p_synapse,
                nu,
                p.xd[31] * beta_h_t * bcwd_VGCC,
                p.xd[31] * max_h_t * bcwd_VGCC,
                typemax(Float64)
            ), # 67
            @j_jump(
                68,
                p_synapse,
                nu,
                p.xd[30] * alpha_h_t * frwd_VGCC,
                p.xd[30] * max_h_t * frwd_VGCC,
                typemax(Float64)
            ), # 68
            @j_jump(
                69,
                p_synapse,
                nu,
                p.xd[32] * beta_h_t * bcwd_VGCC,
                p.xd[32] * max_h_t * bcwd_VGCC,
                typemax(Float64)
            ), # 69
            @j_jump(
                70,
                p_synapse,
                nu,
                p.xd[31] * alpha_m_t * frwd_VGCC,
                p.xd[31] * max_m_t * frwd_VGCC,
                typemax(Float64)
            ), # 70
            @j_jump(
                71,
                p_synapse,
                nu,
                p.xd[32] * beta_m_t * bcwd_VGCC,
                p.xd[32] * max_m_t * bcwd_VGCC,
                typemax(Float64)
            ), # 71, this one can have a high rate

            ################### L-type VGCC  ###################
            @j_jump(
                72,
                p_synapse,
                nu,
                p.xd[33] * alpha_l * frwd_VGCC,
                p.xd[33] * max_alpha_l * frwd_VGCC,
                typemax(Float64)
            ), # 72
            @j_jump(
                73,
                p_synapse,
                nu,
                p.xd[34] * beta_1_l * bcwd_VGCC,
                p.xd[34] * max_beta_1_l * bcwd_VGCC,
                typemax(Float64)
            ), # 73
            @j_jump(
                74,
                p_synapse,
                nu,
                p.xd[33] * alpha_l * frwd_VGCC,
                p.xd[33] * max_alpha_l * frwd_VGCC,
                typemax(Float64)
            ), # 74
            @j_jump(
                75,
                p_synapse,
                nu,
                p.xd[35] * beta_2_l * bcwd_VGCC,
                p.xd[35] * max_beta_2_l * bcwd_VGCC,
                typemax(Float64)
            ), # 75

            ################### LTD/LTP  ###################
            @j_jump(76, p_synapse, nu, p.xd[36] * D_rate, 1, typemax(Float64)), # 76
            @j_jump(77, p_synapse, nu, p.xd[37] * P_rate, 1, typemax(Float64)), # 77
            @j_jump(78, p_synapse, nu, p.xd[36] * P_rate, 1, typemax(Float64)), # 78
            @j_jump(79, p_synapse, nu, p.xd[38] * D_rate, 1, typemax(Float64)), # 79
        ],
    )

    return jumps
end

function buildRxDependencyGraph(nu)
    numrxs, _ = size(nu)
    dep_graph = [Vector{Int}() for n = 1:(numrxs-1)]
    for rx = 1:numrxs
        if rx == 55  # no need to track the Poisson process
            continue
        end
        rx_ix = rx
        if 56 <= rx < 80
            rx_ix += 19
        elseif rx >= 80
            rx_ix -= 25
        end
        for (spec, _) in zip(findnz(nu[rx, :])...)
            # we need to reorder the indices according to the order
            # they apper in the problem
            for (dependent_rx, _) in zip(findnz(nu[:, spec])...)
                # we need to reorder the indices according to the order
                # they apper in the problem
                dependent_rx_ix = dependent_rx
                if 56 <= dependent_rx < 80
                    dependent_rx_ix += 19
                elseif dependent_rx >= 80
                    dependent_rx_ix -= 25
                end
                push!(dep_graph[rx_ix], dependent_rx_ix)
            end
        end
    end
    return dep_graph
end

function SynapseProblem(
    xc,
    xd,
    t1,
    t2,
    events_bap,
    bap_by_epsp,
    glu,
    p_synapse,
    nu,
    algo::T,
    agg = nothing;
    saveat = [],
    save_everystep = isempty(saveat),
    kwargs...,
) where {T<:CHV}
    problem = PDMP.PDMPProblem(
        (dxc, xc, xd, p, t) -> F_synapse(dxc, xc, xd, p, t, events_bap, bap_by_epsp),
        (rate, xc, xd, p, t, sum_rate) -> R_synapse(rate, xc, xd, p, t, sum_rate, glu),
        nu,
        xc,
        xd,
        p_synapse,
        (t1, t2);
        Ncache = 12, # this option is for AD in PreallocationTools
    )
    sol = solve(problem, algo; kwargs...)
    return sol
end

function SynapseProblem(
    xc,
    xd,
    t1,
    t2,
    events_bap,
    bap_by_epsp,
    glu,
    p_synapse,
    nu,
    algo,
    agg;
    jumps = nothing,
    save_positions = (false, true),
    saveat = [],
    save_everystep = isempty(saveat),
    kwargs...,
)
    p = (
        xd0 = copy(xd),
        xd = copy(xd),
        Glu = p_synapse.glu_amp * glu,
        p_synapse = p_synapse,
    )
    oprob = ODEProblem(
        (dxc, xc, p, t) ->
            F_synapse(dxc, xc, p.xd, p.p_synapse, t, events_bap, bap_by_epsp),
        xc,
        (t1, t2),
        p,
    )
    xdsol = SavedValues(typeof(t1), typeof(xd))
    dep_graph = buildRxDependencyGraph(nu)
    callback = _SavingCallback(
        (u, t, integrator) -> copy(integrator.p.xd),
        xdsol;
        save_modified = typeof(save_positions) <: Bool ? save_positions : save_positions[2],
    )
    jprob = JumpProblem(
        oprob,
        agg,
        jumps;
        dep_graph,
        save_positions,
        saveat,
        save_everystep,
        callback,
    )
    sol = (xcsol = solve(jprob, algo; saveat, save_everystep, kwargs...), xdsol = xdsol)
    # @info "Integrator" sol.xcsol.stats
    # with tweaked JumpProcesses.jl
    # total_jumps = (jprob.discrete_jump_aggregation.rejections + jprob.discrete_jump_aggregation.jumps)
    # rejections = jprob.discrete_jump_aggregation.rejections[total_jumps .> 0]
    # total_jumps = total_jumps[total_jumps .> 0]
    # rejection_rate = sum(rejections ./ total_jumps) ./ length(total_jumps)
    # @info "Rejection rate" rejection_rate
    return sol
end
