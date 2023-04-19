macro variable_rate_jump(i, nu, ex)

	assignments = Expr[]    

	alpha_beta_regex = r"(alpha|beta)_(m_r|h_r|m_t|l|1_l|2_l)"
	alpha_beta_matches = Set([m.match for m in eachmatch(alpha_beta_regex, "$ex")])

	if length(alpha_beta_matches) > 0

		for m in ("alpha_1_l", "alpha_2_l", "beta_l")
			if m in alpha_beta_matches
				throw(DomainError(m, "this variable does not exist in the model."))
			end
		end

		push!(assignments, :(Vsp = u[51]))

		if "alpha_m_r" in alpha_beta_matches || "beta_m_r" in alpha_beta_matches
			push!(assignments, :(alpha_m_r, beta_m_r = rates_m_r(Vsp)))
		end

		if "alpha_h_r" in alpha_beta_matches || "beta_h_r" in alpha_beta_matches
			push!(assignments, :(alpha_h_r, beta_h_r = rates_h_r(Vsp)))
		end

		if "alpha_m_t" in alpha_beta_matches || "beta_m_t" in alpha_beta_matches
			push!(assignments, :(alpha_m_t, beta_m_t = rates_m_t(Vsp)))
		end

		if "alpha_l" in alpha_beta_matches || "beta_1_l" in alpha_beta_matches || "beta_2_l" in alpha_beta_matches
			push!(assignments, :(alpha_l, beta_1_l, beta_2_l = rates_l(Vsp)))
		end

	end

	if occursin("D_rate", "$ex")
		push!(assignments, :(D_rate = plasticityRate(u[27], 2, K_D) / t_D))
	end

	if occursin("P_rate", "$ex")
		push!(assignments, :(P_rate = plasticityRate(u[28], 2, K_D) / t_P))
	end

	block_assignments = Expr(:block, assignments...)

	quote
		js, _ = findnz($(esc(nu))[$(esc(i)), :])
		function rate(u, p, t)
			@unpack_SynapseParams p[end]
			$(esc(block_assignments))
			return $ex
		end
		function affect!(u, p, t)
			for j in js
				u[j] += $(esc(nu))[$(esc(i)), j]
			end
		end
		# function urate(u, p, t)
		# end
		# function rateinterval(u, p, t)
		# end
		# return VariableRateJump(rate, affect!; urate=urate, rateinterval=rateinterval)
		VariableRateJump(rate, affect!)
	end

end


function synapse_jumps(p_synapse::SynapseParams, nu, glu)

	@unpack_SynapseParams p_synapse

	############### Glutamate & GABA ###################
	Glu = glu_amp * glu

	p = (
		############### AMPA ###################
		#2line-GO
		4 * AMPA_k1 * Glu, # 1
		3 * AMPA_k1 * Glu, # 2
		2 * AMPA_k1 * Glu, # 3
		1 * AMPA_k1 * Glu, # 4
		#2line-BACK
		4 * AMPA_k_1, # 5
		3 * AMPA_k_1, # 6
		2 * AMPA_k_1, # 7
		1 * AMPA_k_1, # 8
		#3line-GO
		3 * AMPA_k1 * Glu, # 9
		3 * AMPA_k1 * Glu, # 10
		2 * AMPA_k1 * Glu, # 11
		1 * AMPA_k1 * Glu, # 12
		#3line-BACK
		3 * AMPA_k_1, # 13
		2 * AMPA_k_1, # 14
		1 * AMPA_k_1, # 15
		1 * AMPA_k_2, # 16
		#4line-GO
		2 * AMPA_k1 * Glu, # 17
		1 * AMPA_k1 * Glu, # 18
		#4line-BACK
		2 * AMPA_k_1, # 19
		1 * AMPA_k_1, # 20
		#1column-GO-BACK
		4 * AMPA_delta_0, # 21
		1 * AMPA_gamma_0, # 22
		#2column-GO-BACK
		1 * AMPA_delta_1, # 23
		1 * AMPA_gamma_1, # 24
		#3column-GO
		1 * AMPA_alpha, # 25
		2 * AMPA_delta_1, # 26
		1 * AMPA_delta_2, # 27
		#3column-BACK
		1 * AMPA_gamma_2, # 28
		1 * AMPA_gamma_1, # 29
		2 * AMPA_beta, # 30
		#4column-GO
		1 * AMPA_alpha, # 31
		3 * AMPA_delta_1, # 32
		2 * AMPA_delta_2, # 33
		#4column-BACK
		1 * AMPA_gamma_2, # 34
		1 * AMPA_gamma_1, # 35
		2 * AMPA_beta, # 36
		#5column-GO
		1 * AMPA_alpha, # 37
		4 * AMPA_delta_1, # 38
		3 * AMPA_delta_2, # 39
		#5column-BACK
		1 * AMPA_gamma_2, # 40
		1 * AMPA_gamma_1, # 41
		4 * AMPA_beta, # 42

		############### NMDA ###################
		#1line-GO
		NMDA_N2A_ka * Glu, # 43
		NMDA_N2A_kb * Glu, # 44
		NMDA_N2A_kc, # 45
		NMDA_N2A_kd, # 46
		NMDA_N2A_ke, # 47
		NMDA_N2A_kf, # 48

		#1line-BACK
		NMDA_N2A_k_f, # 49
		NMDA_N2A_k_e, # 50
		NMDA_N2A_k_d, # 51
		NMDA_N2A_k_c, # 52
		NMDA_N2A_k_b, # 53
		NMDA_N2A_k_a, # 54

		############### NMDA GLUN2B ###################
		#1line-GO
		NMDA_N2B_sa * Glu, # 80
		NMDA_N2B_sb * Glu, # 81
		NMDA_N2B_sc,  # 82
		NMDA_N2B_sd,  # 93
		NMDA_N2B_se,  # 84
		NMDA_N2B_sf, # 85

		#1line-BACK
		NMDA_N2B_s_f,  # 86
		NMDA_N2B_s_e,  # 87
		NMDA_N2B_s_d,  # 88
		NMDA_N2B_s_c,  # 89
		NMDA_N2B_s_b,  # 90
		NMDA_N2B_s_a,  # 91

		############### GABA ###################
		GABA_r_b1 * Glu, # 92, to simplify, we use the same ammount at the same time
		GABA_r_u1,  # 93
		GABA_r_b2 * Glu, # 94
		GABA_r_u2,  # 95
		GABA_r_ro1, # 96
		GABA_r_c1,  # 97
		GABA_r_ro2, # 98
		GABA_r_c2,  # 99

		############# Original pararameters ############
		p_synapse,

	)

	ma_jumps_idx = vcat(1:54, 80:99)

	reactant_stoich = Vector{Vector{Pair{Int,Int}}}(undef, length(ma_jumps_idx))
	reactants = (nu .< 0)
	for i in ma_jumps_idx
		js = findall(reactants[i, :])
		stoich = [j => -nu[i, j] for j in js]
		stoich_i = i - (i > 54 ? 25 : 0)
		reactant_stoich[stoich_i] = stoich
	end

	nets = (nu .> 0)
	net_stoich = Vector{Vector{Pair{Int,Int}}}(undef, length(ma_jumps_idx))
	for i in ma_jumps_idx
		js, _ = findnz(nets[i, :])
		stoich = [j => nu[i, j] for j in js]
		stoich_i = i - (i > 54 ? 25 : 0)
		net_stoich[stoich_i] = stoich
	end

	param_idxs = 1:length(ma_jumps_idx)

	jumps = [
		MassActionJump(reactant_stoich, net_stoich; scale_rates = false, param_idxs),

		################### R-type VGCC ###################
		@variable_rate_jump(56, nu, u[25] * alpha_m_r * frwd_VGCC), # 56
		@variable_rate_jump(57, nu, u[26] * beta_m_r  * bcwd_VGCC), # 57
		@variable_rate_jump(58, nu, u[25] * alpha_h_r * frwd_VGCC), # 58
		@variable_rate_jump(59, nu, u[27] * beta_h_r  * bcwd_VGCC), # 59
		@variable_rate_jump(60, nu, u[26] * alpha_h_r * frwd_VGCC), # 60
		@variable_rate_jump(61, nu, u[28] * beta_h_r  * bcwd_VGCC), # 61
		@variable_rate_jump(62, nu, u[27] * alpha_m_r * frwd_VGCC), # 62
		@variable_rate_jump(63, nu, u[28] * beta_m_r  * bcwd_VGCC), # 63


		################### T-type VGCC  ###################
		@variable_rate_jump(64, nu, u[29] * alpha_m_t * frwd_VGCC), # 64
		@variable_rate_jump(65, nu, u[30] * beta_m_t  * bcwd_VGCC), # 65 this one can have a high rate
		@variable_rate_jump(66, nu, u[29] * alpha_h_t * frwd_VGCC), # 66
		@variable_rate_jump(67, nu, u[31] * beta_h_t  * bcwd_VGCC), # 67
		@variable_rate_jump(68, nu, u[30] * alpha_h_t * frwd_VGCC), # 68
		@variable_rate_jump(69, nu, u[32] * beta_h_t  * bcwd_VGCC), # 69
		@variable_rate_jump(70, nu, u[31] * alpha_m_t * frwd_VGCC), # 70
		@variable_rate_jump(71, nu, u[32] * beta_m_t  * bcwd_VGCC), # 71, this one can have a high rate


		################### L-type VGCC  ###################
		@variable_rate_jump(72, nu, u[33] * alpha_l  * frwd_VGCC), # 72
		@variable_rate_jump(73, nu, u[34] * beta_1_l * bcwd_VGCC), # 73
		@variable_rate_jump(74, nu, u[33] * alpha_l  * frwd_VGCC), # 74
		@variable_rate_jump(75, nu, u[35] * beta_2_l * bcwd_VGCC), # 75

		################### LTD/LTP  ###################
		@variable_rate_jump(76, nu, u[36] * D_rate), # 76
		@variable_rate_jump(77, nu, u[37] * P_rate), # 77
		@variable_rate_jump(78, nu, u[36] * P_rate), # 78
		@variable_rate_jump(79, nu, u[38] * D_rate), # 79
	]

	return p, jumps
end

function synapse_f!(du, u, p_synapse::SynapseParams, t, events_bap, bap_by_epsp)
	@unpack_SynapseParams p_synapse

	##### Stochastic channels/receptors
	n1_ampa   = u[14] # ampa subconductance 1
	n2_ampa   = u[15] # ampa subconductance 2
	n3_ampa   = u[16] # ampa subconductance 3
	n1_nmda_A = u[22] # nmda subconductance 1
	n2_nmda_A = u[23] # nmda subconductance 2
	n1_nmda_B = u[44] # nmda subconductance 1
	n2_nmda_B = u[45] # nmda subconductance 2
	n_car     = u[28] # vgcc-R opened state
	n_cat     = u[32] # vgcc-T opened state
	n_cal     = u[34] + u[35] # vgcc-L opened states
	n_gaba1   = u[49] # GABA opened state
	n_gaba2   = u[50] # GABA opened state

	##### Continuous variables
	Vsp, Vdend, Vsoma, λ, ImbufCa, Ca, Dye, CaM0, CaM2C, CaM2N, CaM4, mCaN,
	CaN4, mKCaM, KCaM0, KCaM2N, KCaM2C, KCaM4, PCaM0, PCaM2C, PCaM2N, PCaM4, P, P2,
	LTD, LTP, act_D, act_P, m, h, n, SK, λ_age, λ_aux = u[51:end]

	##### plasticity prediction regions
	CaMKII     = KCaM0 + KCaM2C + KCaM2N + KCaM4 + PCaM0 + PCaM2C + PCaM2N + PCaM4 + P + P2
	CaN        = CaN4

	#### activation when it is inside the region
	∂LTD     = SVector(CaN, CaMKII) ∈ LTD_region
	∂LTP     = SVector(CaN, CaMKII) ∈ LTP_region

	∂act_D     = a_D * ∂LTD - b_D * act_D * (1 - ∂LTD)
	∂act_P     = a_P * ∂LTP - b_P * act_P * (1 - ∂LTP)

	##### Na channel
	m_inf = alpha_m(Vsoma) / (alpha_m(Vsoma) + beta_m(Vsoma))
	m_tau = 1 / (alpha_m(Vsoma) + beta_m(Vsoma))
	∂m    = (m_inf - m) / m_tau
	∂h    = alpha_h(Vsoma) * (1 - h) - beta_h(Vsoma) * h
	I_Na  = gamma_Na * (m^3) * h * (Erev_Na - Vsoma)


	##### K channel
	n_inf = 1 / (1 + alpha_n(Vsoma) )
	n_tau = max(50 * beta_n(Vsoma) / (1 + alpha_n(Vsoma)), 2.)
	∂n    = (n_inf - n) / n_tau
	I_K   = gamma_K * n * (Erev_K - Vsoma)

	##### NMDA
	NMDA  = (n1_nmda_A + n2_nmda_A + n1_nmda_B + n2_nmda_B) * B(Vsp, Mg) * gamma_nmda
	Inmda = (Erev_nmda - Vsp) * NMDA # current nmda

	##### AMPA
	Iampa = (Erev_ampa - Vsp) * (gamma_ampa1 * n1_ampa + gamma_ampa2 * n2_ampa + gamma_ampa3 * n3_ampa) # current ampa

	##### GABA
	Igaba = (n_gaba1 + n_gaba2) *  (Erev_Cl - Vdend)  * gamma_GABA

	##### Calcium sources (VGCCs currents, and NMDA calcium contribution)
	ΦCa      = perm * ghk(Vsp, Ca, Ca_ext, p_synapse) #GHK factor
	Ica_nmda = f_Ca * ΦCa * NMDA
	Icar     = gamma_CaR * n_car * ΦCa
	Icat     = gamma_CaT * n_cat * ΦCa
	Ical     = gamma_CaL * n_cal * ΦCa

	##### SK channel (not stochastic)
	∂SK = (SK_chnnl(Ca) * frwd_SK - SK) / (SK_time * bcwd_SK) #SK spine
	Isk = SK_gamma * (SK_Erev - Vsp) * SK * N_SK

	##### Backpropgation
	# Post input - for experimentally induced BaPs and those induced by EPSPs
	I_BaP  = inputBaP(t, bap_by_epsp, injbap, I_clamp) + inputBaP(t, events_bap,  injbap, I_clamp)
	# Bap decay/attenuation - two component for adaptation in the Bap
	∂λ     = (1-λ)/trec - delta_decay * (1/λ_aux) * λ * I_BaP
	∂λ_aux = (1-λ_aux)/trec - delta_aux * λ_aux * I_BaP
	gadapt = λ * g_diff * ϕ_dist

	# Bap decay/attenuation - age dependent modification factor
	∂λ_age = (1-λ_age)/trec_soma - delta_soma * λ_age * I_BaP

	##### Voltage
	# Spine
	∂Vsp   = (Isk + Inmda + Iampa + Icat + Icar + Ical + g_neck * (Vdend - Vsp) + g_leak * (E_leak - Vsp)) / (Csp)
	# Dendrite
	∂Vdend = (g_neck * (Vsp - Vdend) + Igaba + g_leakdend * (E_leak - Vdend) + gadapt * (Vsoma - Vdend)) / Cdend
	# Soma
	∂Vsoma = ((I_BaP + I_Na)*λ_age + I_K + g_leaksoma * (E_leak - Vsoma) +  gadapt * (Vdend - Vsoma)) / Csoma


	##### Buffer and dye (spine only - no neck diffusion)
	∂ImbufCa = Imbuf_k_on * (Imbuf_con - ImbufCa) * Ca - Imbuf_k_off * ImbufCa
	∂Dye     = 4*fluo5f_kf * (fluo5f_con - Dye) * Ca - 8*fluo5f_kb * Dye

	##### Ca Downstream
	### CaM-KCaM-rates (coarsed model) from Pepke adapted by
	kf_2C      = rates_adapt(kon_1C,   kon_2C,   koff_1C,  kon_2C, Ca)
	kb_2C      = rates_adapt(koff_1C,  koff_2C,  koff_1C,  kon_2C, Ca)
	kf_2N      = rates_adapt(kon_1N,   kon_2N,   koff_1N,  kon_2N, Ca)
	kb_2N      = rates_adapt(koff_1N,  koff_2N,  koff_1N,  kon_2N, Ca)
	kf_K2C     = rates_adapt(kon_K1C,  kon_K2C,  koff_K1C, kon_K2C, Ca)
	kb_K2C     = rates_adapt(koff_K1C, koff_K2C, koff_K1C, kon_K2C, Ca)
	kf_K2N     = rates_adapt(kon_K1N,  kon_K2N,  koff_K1N, kon_K2N, Ca)
	kb_K2N     = rates_adapt(koff_K1N, koff_K2N, koff_K1N, kon_K2N, Ca)
	F          = CaMKII / mKCaM_con

	∂CaM0 = k2*PCaM0 + kb_2C*CaM2C + kb_2N*CaM2N + kb_CaM0*KCaM0 +
		-(1//2)*kf_2C*(Ca^2)*CaM0 - (1//2)*kf_2N*(Ca^2)*CaM0 +
		-kf_CaM0*CaM0*mKCaM

	∂CaM2C = kb_2N*CaM4 + kb_CaM2C*KCaM2C + k2*PCaM2C +
		+(1//2)*kf_2C*(Ca^2)*CaM0 - kb_2C*CaM2C - (1//2)*kf_2N*(Ca^2)*CaM2C +
		-kf_CaM2C*CaM2C*mKCaM

	∂CaM2N = kb_2C*CaM4 + kb_CaM2N*KCaM2N + k2*PCaM2N +
		+(1//2)*kf_2N*(Ca^2)*CaM0 - kb_2N*CaM2N - (1//2)*kf_2C*(Ca^2)*CaM2N
		-kf_CaM2N*CaM2N*mKCaM

	∂CaM4 = k2*PCaM4 + kcanb*CaN4 + kb_CaM4*KCaM4 +
		+(1//2)*kf_2N*(Ca^2)*CaM2C + (1//2)*kf_2C*(Ca^2)*CaM2N - kb_2C*CaM4
		-kb_2N*CaM4 - kcanf*CaM4*mCaN - kf_CaM4*CaM4*mKCaM

	∂mCaN = kcanb*CaN4 - kcanf*CaM4*mCaN

	∂CaN4 = kcanf*CaM4*mCaN - kcanb*CaN4

	∂mKCaM = kb_CaM0*KCaM0 + k3*P + kb_CaM2C*KCaM2C + kb_CaM2N*KCaM2N +
		+kb_CaM4*KCaM4 - kf_CaM0*CaM0*mKCaM - kf_CaM2C*CaM2C*mKCaM +
		-kf_CaM2N*CaM2N*mKCaM - kf_CaM4*CaM4*mKCaM

	∂KCaM0 = kb_K2C*KCaM2C + kb_K2N*KCaM2N + kf_CaM0*CaM0*mKCaM +
		-kb_CaM0*KCaM0 - (1//2)*kf_K2C*(Ca^2)*KCaM0 - F*k1*KCaM0 +
		-(1//2)*kf_K2N*(Ca^2)*KCaM0

	∂KCaM2N = kb_K2C*KCaM4 + kf_CaM2N*CaM2N*mKCaM +
		+(1//2)*kf_K2N*(Ca^2)*KCaM0 - kb_CaM2N*KCaM2N - kb_K2N*KCaM2N +
		-(1//2)*kf_K2C*(Ca^2)*KCaM2N - F*k1*KCaM2N

	∂KCaM2C = kb_K2N*KCaM4 + kf_CaM2C*CaM2C*mKCaM +
		+(1//2)*kf_K2C*(Ca^2)*KCaM0 - kb_CaM2C*KCaM2C - kb_K2C*KCaM2C +
		-F*k1*KCaM2C - (1//2)*kf_K2N*(Ca^2)*KCaM2C

	∂KCaM4 = kf_CaM4*CaM4*mKCaM + (1//2)*kf_K2C*(Ca^2)*KCaM2N +
		+(1//2)*kf_K2N*(Ca^2)*KCaM2C - kb_CaM4*KCaM4 - kb_K2C*KCaM4 +
		-kb_K2N*KCaM4 - F*k1*KCaM4

	∂PCaM0 = F*k1*KCaM0 - k2*PCaM0

	∂PCaM2N = F*k1*KCaM2N - k2*PCaM2N

	∂PCaM2C = F*k1*KCaM2C - k2*PCaM2C

	∂PCaM4 = F*k1*KCaM4 - k2*PCaM4

	∂P = k2*PCaM0 + k5*P2 + k2*PCaM2C + k2*PCaM2N + k2*PCaM4 - k3*P - k4*P

	∂P2 = k4*P - k5*P2


	### Postsynaptic Ca
	∂Ca = (Ca_infty - Ca) / tau_ca +
		+(Ica_nmda + Icar + Ical + Icat) / (2 * faraday * A_sp) +
		+(max(Ca_infty, Ca/3) - Ca) / tau_diff +
		-∂ImbufCa +
		-∂Dye +
		+2kb_2C*CaM2C + 2kb_2C*CaM4 + 2kb_2N*CaM2N + 2kb_2N*CaM4 +
		+2kb_K2C*KCaM2C + 2kb_K2N*KCaM2N + 2kb_K2C*KCaM4 + 2kb_K2N*KCaM4 +
		-kf_2C*(Ca^2)*CaM0 - kf_2N*(Ca^2)*CaM0 - kf_2N*(Ca^2)*CaM2C +
		-kf_2C*(Ca^2)*CaM2N - kf_K2C*(Ca^2)*KCaM0 - kf_K2C*(Ca^2)*KCaM2N +
		-kf_K2N*(Ca^2)*KCaM0 - kf_K2N*(Ca^2)*KCaM2C

	### Xdot update
	u[51] =  ∂Vsp
	u[52] =  ∂Vdend
	u[53] =  ∂Vsoma
	u[54] =  ∂λ
	u[55] =  ∂ImbufCa
	u[56] =  ∂Ca
	u[57] =  ∂Dye
	u[58] =  ∂CaM0
	u[59] =  ∂CaM2C
	u[60] = ∂CaM2N
	u[61] = ∂CaM4
	u[62] = ∂mCaN
	u[63] = ∂CaN4
	u[64] = ∂mKCaM
	u[65] = ∂KCaM0
	u[66] = ∂KCaM2N
	u[67] = ∂KCaM2C
	u[68] = ∂KCaM4
	u[69] = ∂PCaM0
	u[70] = ∂PCaM2C
	u[71] = ∂PCaM2N
	u[72] = ∂PCaM4
	u[73] = ∂P
	u[74] = ∂P2
	u[75] = ∂LTD
	u[76] = ∂LTP
	u[77] = ∂act_D
	u[78] = ∂act_P
	u[79] = ∂m
	u[80] = ∂h
	u[81] = ∂n
	u[82] = ∂SK
	u[83] = ∂λ_age
	u[84] = ∂λ_aux

end

function buildDependencyGraph(nu)
	numrxs, _ = size(nu)
	dep_graph = [Vector{Int}() for n in 1:numrxs]
	for rx in 1:numrxs
		for (spec, _) in zip(findnz(nu[rx, :])...)
			for (dependent_rx, _) in zip(findnz(nu[:, spec])...)
				push!(dep_graph[rx], dependent_rx)
			end
		end
	end
	return dep_graph
end

function pdmpsynapse(xc, xd, t1, t2, events_bap, bap_by_epsp, glu, p_synapse, nu, agg::Direct; algo=Tsit5(), save_positions = (false, true), kwargs...)
	u = vcat(xc, xd)
	p, jumps = synapse_jumps(p_synapse, nu, glu)
	oprob = ODEProblem(
		(u, p, t) -> synapse_f!(u, p_synapse, t, events_bap, bap_by_epsp),
		u,
		(t1, t2),
		p
	)
	jprob = JumpProblem(oprob, agg, jumps...; save_positions = save_positions)
	return solve(jprob, algo; kwargs...)
end

function pdmpsynapse(xc, xd, t1, t2, events_bap, bap_by_epsp, glu, p_synapse, nu, agg::CoevolveSynced; algo=Tsit5(), save_positions = (false, true), kwargs...)
	u = vcat(xc, xd)
	jumps, p = synapse_jumps(p_synapse, nu, glu)
	oprob = ODEProblem(
		(u, p, t) -> synapse_f!(u, p_synapse, t, events_bap, bap_by_epsp),
		u,
		(t1, t2),
		p
	)
	dep_graph = buildDependencyGraph(nu)
	jprob = JumpProblem(oprob, agg, jumps...; dep_grap = dep_graph, save_positions = save_positions)
	return solve(jprob, algo; kwargs...)
end
