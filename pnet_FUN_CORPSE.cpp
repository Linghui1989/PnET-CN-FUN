#include "pnet_model.h"

namespace {
	double sum_compartment_c(const pnet_model::FUN_CORPSE_Compartment* comp);
	double sum_compartment_n(const pnet_model::FUN_CORPSE_Compartment* comp);
	double sum_compartment_unprotected_n(const pnet_model::FUN_CORPSE_Compartment* comp);
	double sum_fun_unprotected_n(const pnet_model::share_struct* share);
	double remove_fun_unprotected_n(pnet_model::share_struct* share, double requestedN);
	double compute_current_total_n(const pnet_model::share_struct* share);
	void log_fun_compartment_n_diag(const char* label, int year, int doy,
		double before_total_n, double after_total_n);

	double clamp01(double value)
	{
		if (value < 0.0) return 0.0;
		if (value > 1.0) return 1.0;
		return value;
	}

	double compute_fixation_cost(double soilT_C)
	{
		const double s_FIX = -30.0;
		const double a_FIX = -3.62;
		const double b_FIX = 0.27;
		const double c_FIX = 25.15;
		return s_FIX * exp(a_FIX + b_FIX * soilT_C * (1.0 - 0.5 * soilT_C / c_FIX)) - 2.0 * s_FIX;
	}

	void repartition_pair(double* rhizoValue, double* bulkValue, double rhizoFrac)
	{
		const double total = std::max(0.0, *rhizoValue + *bulkValue);
		const double targetRhizo = total * clamp01(rhizoFrac);
		*rhizoValue = targetRhizo;
		*bulkValue = std::max(0.0, total - targetRhizo);
	}

	void zero_compartment(pnet_model::FUN_CORPSE_Compartment* comp)
	{
		for (int i = 0; i < pnet_model::BG_CLASS_COUNT; ++i)
		{
			comp->unprotC[i] = 0.0;
			comp->unprotN[i] = 0.0;
			comp->protC[i] = 0.0;
			comp->protN[i] = 0.0;
			comp->lastDecompC[i] = 0.0;
			comp->lastDecompN[i] = 0.0;
		}

		comp->livingMicrobeC = 0.0;
		comp->livingMicrobeN = 0.0;
		comp->CO2 = 0.0;
		comp->Rtot = 0.0;
		comp->lastGrossMineralization = 0.0;
		comp->lastGrossImmobilization = 0.0;
		comp->lastProtectionN = 0.0;
		comp->lastDeprotectionN = 0.0;
		comp->lastTurnoverN = 0.0;
	}

	double sum_compartment_c(const pnet_model::FUN_CORPSE_Compartment* comp)
	{
		double total = comp->livingMicrobeC;
		for (int i = 0; i < pnet_model::BG_CLASS_COUNT; ++i)
		{
			total += comp->unprotC[i] + comp->protC[i];
		}
		return total;
	}

	double sum_compartment_n(const pnet_model::FUN_CORPSE_Compartment* comp)
	{
		double total = comp->livingMicrobeN;
		for (int i = 0; i < pnet_model::BG_CLASS_COUNT; ++i)
		{
			total += comp->unprotN[i] + comp->protN[i];
		}
		return total;
	}

	double sum_compartment_unprotected_n(const pnet_model::FUN_CORPSE_Compartment* comp)
	{
		double total = 0.0;
		for (int i = 0; i < pnet_model::BG_CLASS_COUNT; ++i)
		{
			total += std::max(0.0, comp->unprotN[i]);
		}
		return total;
	}

	double sum_fun_unprotected_n(const pnet_model::share_struct* share)
	{
		return sum_compartment_unprotected_n(&share->funBG.organicLayer)
			+ sum_compartment_unprotected_n(&share->funBG.mineralRhizosphere)
			+ sum_compartment_unprotected_n(&share->funBG.mineralBulk);
	}

	double remove_fun_unprotected_n(pnet_model::share_struct* share, double requestedN)
	{
		if (requestedN <= 0.0)
		{
			return 0.0;
		}

		pnet_model::FUN_CORPSE_Compartment* comps[3] = {
			&share->funBG.organicLayer,
			&share->funBG.mineralRhizosphere,
			&share->funBG.mineralBulk
		};

		double totalAvailable = 0.0;
		for (int c = 0; c < 3; ++c)
		{
			totalAvailable += sum_compartment_unprotected_n(comps[c]);
		}

		if (totalAvailable <= 0.0)
		{
			return 0.0;
		}

		double targetRemoval = std::min(requestedN, totalAvailable);
		double removedTotal = 0.0;

		for (int c = 0; c < 3; ++c)
		{
			for (int i = 0; i < pnet_model::BG_CLASS_COUNT; ++i)
			{
				const double poolN = std::max(0.0, comps[c]->unprotN[i]);
				if (poolN <= 0.0)
				{
					continue;
				}

				const double poolShare = poolN / totalAvailable;
				const double removeN = std::min(poolN, targetRemoval * poolShare);
				comps[c]->unprotN[i] -= removeN;
				removedTotal += removeN;
			}
		}

		return removedTotal;
	}

	double compute_current_total_n(const pnet_model::share_struct* share)
	{
		const double folN = std::max(0.0, share->FolMass) * std::max(0.0, share->FolNCon) / 100.0;
		return folN
			+ std::max(0.0, share->WoodMassN)
			+ std::max(0.0, share->RootMassN)
			+ std::max(0.0, share->HON)
			+ std::max(0.0, share->NH4)
			+ std::max(0.0, share->NO3)
			+ std::max(0.0, share->BudN)
			+ std::max(0.0, share->DeadWoodN)
			+ std::max(0.0, share->PlantN);
	}

	void log_fun_compartment_n_diag(const char* label, int year, int doy,
		double before_total_n, double after_total_n)
	{
		const double dN = after_total_n - before_total_n;
		if (fabs(dN) > 1.0e-8)
		{
			printf("FUN NDIAG Compartment %s year %d doy %d dN=%.10f before=%.10f after=%.10f\n",
				label, year, doy, dN, before_total_n, after_total_n);
		}
	}

	void add_split_to_pools(double C, double N, double fastFrac, double slowFrac,
		double* poolC, double* poolN)
	{
		if (C <= 0.0 && N <= 0.0) return;

		const double fastShare = clamp01(fastFrac);
		const double slowShare = clamp01(slowFrac);

		poolC[pnet_model::BG_FAST] += C * fastShare;
		poolC[pnet_model::BG_SLOW] += C * slowShare;
		poolN[pnet_model::BG_FAST] += N * fastShare;
		poolN[pnet_model::BG_SLOW] += N * slowShare;
	}

	double safe_inverse(double value)
	{
		if (value <= 0.0) return 0.0;
		return 1.0 / value;
	}

}

void pnet_model::ReadFUNCORPSEParams(share_struct* share, char* inputname)
{
	FUN_CORPSE_Params* params = &share->funParams;
	char note[400];
	char paramPath[500];
	char fallbackPath[500];
	FILE* fileParams = NULL;

	params->vmaxref[BG_FAST] = 1500.0;
	params->vmaxref[BG_SLOW] = 28.0;
	params->vmaxref[BG_DEADMIC] = 1000.0;
	params->Ea[BG_FAST] = 37000.0;
	params->Ea[BG_SLOW] = 54000.0;
	params->Ea[BG_DEADMIC] = 50000.0;
	params->kC[BG_FAST] = 0.01;
	params->kC[BG_SLOW] = 0.01;
	params->kC[BG_DEADMIC] = 0.01;
	params->eup[BG_FAST] = 0.6;
	params->eup[BG_SLOW] = 0.1;
	params->eup[BG_DEADMIC] = 0.6;
	params->nup[BG_FAST] = 0.3;
	params->nup[BG_SLOW] = 0.3;
	params->nup[BG_DEADMIC] = 0.3;

	params->max_immobilization_rate = 1.0;
	params->gas_diffusion_exp = 2.5;
	params->minMicrobeC = 5e-4;
	params->Tmic = 0.25;
	params->et = 0.6;
	params->CN_microb_default = 10.0;
	params->tProtected = 100.0;
	params->protection_rate[BG_FAST] = 0.6;
	params->protection_rate[BG_SLOW] = 0.0;
	params->protection_rate[BG_DEADMIC] = 4.0;
	params->frac_turnover_min = 0.2;
	params->frac_turnover_slow = 0.3;

	params->leaf_fast_frac = 0.5;
	params->leaf_slow_frac = 0.5;
	params->root_fast_frac = 0.5;
	params->root_slow_frac = 0.5;
	params->deadwood_fast_frac = 0.1;
	params->deadwood_slow_frac = 0.9;
	params->initial_protected_fraction = 0.5;
	params->initial_organic_fraction = 0.1887;
	params->initial_rhizo_fraction = 0.0;
	params->initial_active_n_scalar_organic = 1.0;
	params->initial_active_n_scalar_mineral = 1.0;
	params->initial_bulk_slow_n_scalar = 1.0;
	params->initial_bulk_deadmic_n_scalar = 1.0;
	params->litter_transfer_to_soil = 6.0e-5;
	params->rhizo_fraction = 0.2;
	params->root_inputs_only_rhiz = 0;

	params->f_exudate_Nacq = 1.0;
	params->f_myc_Nacq = 1.0;
	params->C_Nacq_max_frac_GPP = 0.05;
	params->K_Ndef = 0.05;
	params->exudate_CN = 30.0;
	params->Vmax_myc_N = 0.05;
	params->K_myc_C = 0.5;
	params->K_myc_N = 0.05;
	params->eta_myc_to_plant = 1.0;
	params->myc_cost_n = 0.075 / 365.0 / 2.0;
	params->myc_cost_c = 425000.0 / 365.0;
	params->nonmyc_cost_n = 0.6 / 365.0;
	params->nonmyc_cost_c = 20000.0 / 365.0;
	params->storage_cost = 2.0e-3;
	params->litter_fungi_frac = 0.0;
	params->decomp_scalar_global = 1.0;
	params->decomp_scalar_organic = 1.0;
	params->decomp_scalar_rhizo = 1.0;
	params->decomp_scalar_bulk = 1.0;
	params->gross_nmin_cap_frac_yr = 0.12;
	params->nitrification_cap_frac_day = 0.05;
	params->ndrain_cap_frac_day = 0.03;

	params->CN_microb_min = 8.0;
	params->CN_microb_max = 20.0;
	params->K_Navail_CN = 0.05;

	params->use_dynamic_microbe_cn = 1;
	params->enable_root_exudation = 1;
	params->enable_mycorrhizal_uptake = 1;

	strcpy(paramPath, inputname);
	for (int i = (int)strlen(paramPath) - 1; i >= 0; --i)
	{
		if (paramPath[i] == '\\' || paramPath[i] == '/')
		{
			paramPath[i + 1] = '\0';
			break;
		}
	}
	strcat(paramPath, "fun_corpse_params.txt");

	fileParams = fopen(paramPath, "r");
	if (fileParams == NULL)
	{
		sprintf(fallbackPath, "%s", "fun_corpse_params.txt");
		fileParams = fopen(fallbackPath, "r");
	}

	if (fileParams == NULL)
	{
		printf("FUN-CORPSE parameter file not found, using defaults.\n");
		return;
	}

	const char* paramSource = (fileParams != NULL) ? paramPath : fallbackPath;

	for (int i = 0; i < 2; ++i) fgets(note, sizeof(note), fileParams);

	fscanf(fileParams, "%lf", &params->vmaxref[BG_FAST]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->vmaxref[BG_SLOW]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->vmaxref[BG_DEADMIC]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->Ea[BG_FAST]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->Ea[BG_SLOW]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->Ea[BG_DEADMIC]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->kC[BG_FAST]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->kC[BG_SLOW]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->kC[BG_DEADMIC]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->eup[BG_FAST]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->eup[BG_SLOW]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->eup[BG_DEADMIC]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->nup[BG_FAST]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->nup[BG_SLOW]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->nup[BG_DEADMIC]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->max_immobilization_rate); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->gas_diffusion_exp); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->minMicrobeC); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->Tmic); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->et); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->CN_microb_default); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->tProtected); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->protection_rate[BG_FAST]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->protection_rate[BG_SLOW]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->protection_rate[BG_DEADMIC]); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->frac_turnover_min); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->frac_turnover_slow); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->leaf_fast_frac); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->leaf_slow_frac); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->root_fast_frac); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->root_slow_frac); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->deadwood_fast_frac); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->deadwood_slow_frac); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->initial_protected_fraction); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->initial_organic_fraction); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->initial_rhizo_fraction); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->initial_active_n_scalar_organic); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->initial_active_n_scalar_mineral); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->initial_bulk_slow_n_scalar); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->initial_bulk_deadmic_n_scalar); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->litter_transfer_to_soil); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->rhizo_fraction); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%d", &params->root_inputs_only_rhiz); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->f_exudate_Nacq); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->f_myc_Nacq); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->C_Nacq_max_frac_GPP); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->K_Ndef); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->exudate_CN); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->Vmax_myc_N); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->K_myc_C); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->K_myc_N); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->eta_myc_to_plant); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->CN_microb_min); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->CN_microb_max); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%lf", &params->K_Navail_CN); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%d", &params->use_dynamic_microbe_cn); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%d", &params->enable_root_exudation); fgets(note, sizeof(note), fileParams);
	fscanf(fileParams, "%d", &params->enable_mycorrhizal_uptake); fgets(note, sizeof(note), fileParams);

	if (1 == fscanf(fileParams, "%lf", &params->myc_cost_n)) fgets(note, sizeof(note), fileParams);
	if (1 == fscanf(fileParams, "%lf", &params->myc_cost_c)) fgets(note, sizeof(note), fileParams);
	if (1 == fscanf(fileParams, "%lf", &params->nonmyc_cost_n)) fgets(note, sizeof(note), fileParams);
	if (1 == fscanf(fileParams, "%lf", &params->nonmyc_cost_c)) fgets(note, sizeof(note), fileParams);
	if (1 == fscanf(fileParams, "%lf", &params->storage_cost)) fgets(note, sizeof(note), fileParams);
	if (1 == fscanf(fileParams, "%lf", &params->litter_fungi_frac)) fgets(note, sizeof(note), fileParams);
	if (1 == fscanf(fileParams, "%lf", &params->decomp_scalar_global)) fgets(note, sizeof(note), fileParams);
	if (1 == fscanf(fileParams, "%lf", &params->decomp_scalar_organic)) fgets(note, sizeof(note), fileParams);
	if (1 == fscanf(fileParams, "%lf", &params->decomp_scalar_rhizo)) fgets(note, sizeof(note), fileParams);
	if (1 == fscanf(fileParams, "%lf", &params->decomp_scalar_bulk)) fgets(note, sizeof(note), fileParams);
	if (1 == fscanf(fileParams, "%lf", &params->gross_nmin_cap_frac_yr)) fgets(note, sizeof(note), fileParams);
	if (1 == fscanf(fileParams, "%lf", &params->nitrification_cap_frac_day)) fgets(note, sizeof(note), fileParams);
	if (1 == fscanf(fileParams, "%lf", &params->ndrain_cap_frac_day)) fgets(note, sizeof(note), fileParams);

	fclose(fileParams);
	printf("FUN params loaded from %s\n", paramSource);
	printf("FUN init params: activeN_org=%.6f activeN_min=%.6f bulkSlowN=%.6f bulkDeadmicN=%.6f decomp_global=%.6f grossNCap=%.6f nitrifCap=%.6f ndrainCap=%.6f\n",
		params->initial_active_n_scalar_organic,
		params->initial_active_n_scalar_mineral,
		params->initial_bulk_slow_n_scalar,
		params->initial_bulk_deadmic_n_scalar,
		params->decomp_scalar_global,
		params->gross_nmin_cap_frac_yr,
		params->nitrification_cap_frac_day,
		params->ndrain_cap_frac_day);
}

void pnet_model::ReadFUNRhizoSeasonality(share_struct* share, char* inputname)
{
	char note[400];
	char seasonPath[500];
	char fallbackPath[500];
	FILE* fileSeason = NULL;
	const double fallbackFrac = clamp01(share->funParams.rhizo_fraction);

	for (int doy = 1; doy <= 366; ++doy)
	{
		share->FUNRhizoFractionByDOY[doy] = fallbackFrac;
	}
	share->FUNRhizoFractionCurrent = fallbackFrac;
	share->FUNUseRhizoSeasonality = 0;

#ifdef FUN_DISABLE_RHIZO_SEASONALITY
	printf("FUN rhizosphere seasonality disabled at compile time, using fixed rhizo_fraction %.4f.\n", fallbackFrac);
	return;
#endif

	strcpy(seasonPath, inputname);
	for (int i = (int)strlen(seasonPath) - 1; i >= 0; --i)
	{
		if (seasonPath[i] == '\\' || seasonPath[i] == '/')
		{
			seasonPath[i + 1] = '\0';
			break;
		}
	}
	strcat(seasonPath, "fun_rhizo_daily.txt");

	fileSeason = fopen(seasonPath, "r");
	if (fileSeason == NULL)
	{
		sprintf(fallbackPath, "%s", "fun_rhizo_daily.txt");
		fileSeason = fopen(fallbackPath, "r");
	}

	if (fileSeason == NULL)
	{
		printf("FUN rhizosphere seasonality file not found, using fixed rhizo_fraction %.4f.\n", fallbackFrac);
		return;
	}

	for (int i = 0; i < 2; ++i) fgets(note, sizeof(note), fileSeason);

	int doy = 0;
	double rhizoFrac = fallbackFrac;
	int loaded = 0;
	while (2 == fscanf(fileSeason, "%d %lf", &doy, &rhizoFrac))
	{
		if (doy >= 1 && doy <= 366)
		{
			share->FUNRhizoFractionByDOY[doy] = clamp01(rhizoFrac);
			loaded++;
		}
		fgets(note, sizeof(note), fileSeason);
	}
	fclose(fileSeason);

	if (loaded >= 365)
	{
		if (loaded == 365)
		{
			share->FUNRhizoFractionByDOY[366] = share->FUNRhizoFractionByDOY[365];
		}
		share->FUNUseRhizoSeasonality = 1;
		share->FUNRhizoFractionCurrent = share->FUNRhizoFractionByDOY[1];
		printf("Loaded FUN daily rhizosphere seasonality from file (%d DOY values).\n", loaded);
	}
	else
	{
		printf("FUN rhizosphere seasonality file incomplete (%d DOY values), using fixed rhizo_fraction %.4f.\n", loaded, fallbackFrac);
		for (int day = 1; day <= 366; ++day)
		{
			share->FUNRhizoFractionByDOY[day] = fallbackFrac;
		}
	}
}

void pnet_model::ResetFUN_CORPSEFluxDiagnostics(share_struct* share)
{
	share->funBG.rootExudateC = 0.0;
	share->funBG.mycorrhizalC = 0.0;
	share->funBG.mycorrhizalNToPlant = 0.0;
	share->funBG.grossNMobilization = 0.0;
	share->funBG.grossNImmobilization = 0.0;
	share->funBG.netNMineralization = 0.0;
	share->funBG.microbialTurnoverN = 0.0;
	share->funBG.deprotectionN = 0.0;
	share->funBG.protectionN = 0.0;
}

void pnet_model::EstimateFUNAllocationDemand(veg_struct* veg, clim_struct* clim, int rstep, share_struct* share, int CN_Mode)
{
	share->FUNPotentialGrowthC = 0.0;
	share->FUNPotentialNdemand = 0.0;
	share->FUNPotentialRootProdC = 0.0;
	share->FUNPotentialWoodProdC = 0.0;

	if (CN_Mode != 1)
	{
		return;
	}

	double budNNeed = 0.0;
	double rootNNeed = 0.0;
	double woodNNeed = 0.0;
	double storageNNeed = 0.0;
	double budCPotential = 0.0;
	double TMult = (exp(0.1 * (share->Tave - 7.1)) * 0.68);
	double rootCAdd = veg->RootAllocA * (share->dayspan / 365.0) + veg->RootAllocB * share->FolProdCMo;
	double plantCAfterPsn = share->PlantC + share->NetPsnMo - share->FolGRespMo;
	double rootCPool = share->RootC;

	if (share->FolGDD == 1 && veg->FolReten != 1)
	{
		budCPotential = std::max(0.0, share->BudC);
		budNNeed = (share->BudC / veg->CFracBiomass) * share->FolNCon / 100.0;
	}

	if (plantCAfterPsn >= rootCAdd)
	{
		rootCPool += rootCAdd;
	}

	double rootAllocFrac = (share->dayspan / 365.0) * TMult;
	if (rootAllocFrac > 1.0) rootAllocFrac = 1.0;
	if (rootAllocFrac < 0.0) rootAllocFrac = 0.0;
	double rootAllocCMo = rootAllocFrac * rootCPool;
	double rootProdCMoPotential = rootAllocCMo / (1.0 + veg->RootMRespFrac + veg->GRespFrac);

	double woodProdCMoPotential = 0.0;
	if (share->GDDTot >= veg->GDDWoodStart)
	{
		double GDDWoodEff = (share->GDDTot - veg->GDDWoodStart) / (veg->GDDWoodEnd - veg->GDDWoodStart);
		if (GDDWoodEff > 1.0) GDDWoodEff = 1.0;
		if (GDDWoodEff < 0.0) GDDWoodEff = 0.0;
		double delGDDWoodEff = GDDWoodEff - share->OldGDDWoodEff;
		woodProdCMoPotential = share->WoodC * delGDDWoodEff;
	}

	rootNNeed = rootProdCMoPotential / veg->CFracBiomass * veg->RLPctN * share->NRatio;
	woodNNeed = woodProdCMoPotential / veg->CFracBiomass * veg->WLPctN * share->NRatio;
	// Rebuild PlantN storage gradually across the year so FUN uptake can remain
	// active after the initial flush instead of shutting off once growth demand drops.
	storageNNeed = std::max(0.0, veg->MaxNStore - share->PlantN) * (share->dayspan / 365.0);

	share->FUNPotentialNdemand = std::max(0.0, budNNeed + rootNNeed + woodNNeed + storageNNeed);
	share->FUNPotentialGrowthC = std::max(0.0, plantCAfterPsn);
	share->FUNPotentialRootProdC = std::max(0.0, rootProdCMoPotential);
	share->FUNPotentialWoodProdC = std::max(0.0, woodProdCMoPotential);
}

void pnet_model::InitializeFUN_CORPSEBelowground(site_struct* site, share_struct* share)
{
	FUN_CORPSE_Compartment* org = &share->funBG.organicLayer;
	FUN_CORPSE_Compartment* rhizo = &share->funBG.mineralRhizosphere;
	FUN_CORPSE_Compartment* bulk = &share->funBG.mineralBulk;
	FUN_CORPSE_Params* params = &share->funParams;

	zero_compartment(org);
	zero_compartment(rhizo);
	zero_compartment(bulk);
	ResetFUN_CORPSEFluxDiagnostics(share);
	share->FUNPotentialGrowthC = 0.0;
	share->FUNPotentialNdemand = 0.0;
	share->FUNPotentialRootProdC = 0.0;
	share->FUNPotentialWoodProdC = 0.0;

	share->funBG.microbialCN = params->CN_microb_default;

	const double organicFrac = clamp01(params->initial_organic_fraction);
	const double mineralRhizoFrac = clamp01(params->initial_rhizo_fraction);
	const double totalC0 = std::max(0.0, site->HOM);

	// MATLAB single-plot initialization uses explicit pool fractions rather than
	// putting all initial SOM into one slow pool. Use those relative fractions here.
	const double matlabLitterC[BG_CLASS_COUNT] = { 5.6, 1520.0, 3.1 };
	const double matlabProtectedC[BG_CLASS_COUNT] = { 560.0, 0.0, 2400.0 };
	const double matlabLitterN[BG_CLASS_COUNT] = { 0.1, 90.8, 0.35 };
	const double matlabProtectedN[BG_CLASS_COUNT] = { 27.0, 0.0, 130.0 };
	const double matlabMicrobeC = 14.2;
	const double matlabMicrobeN = matlabMicrobeC / std::max(1.0, params->CN_microb_default);
	const double matlabTotalC = matlabLitterC[BG_FAST] + matlabLitterC[BG_SLOW] + matlabLitterC[BG_DEADMIC]
		+ matlabProtectedC[BG_FAST] + matlabProtectedC[BG_SLOW] + matlabProtectedC[BG_DEADMIC]
		+ matlabMicrobeC;
	const double matlabTotalN = matlabLitterN[BG_FAST] + matlabLitterN[BG_SLOW] + matlabLitterN[BG_DEADMIC]
		+ matlabProtectedN[BG_FAST] + matlabProtectedN[BG_SLOW] + matlabProtectedN[BG_DEADMIC]
		+ matlabMicrobeN;
	const double matlabCN = matlabTotalC / std::max(1.0e-12, matlabTotalN);

	const double organicC0 = totalC0 * organicFrac;
	const double mineralC0 = totalC0 - organicC0;
	const double rhizoC0 = mineralC0 * mineralRhizoFrac;
	const double bulkC0 = mineralC0 - rhizoC0;

	const double totalN0 = (site->HON > 0.0)
		? std::max(0.0, site->HON)
		: totalC0 / std::max(1.0e-12, matlabCN);

	const double cFracUnprotected[BG_CLASS_COUNT] = {
		matlabLitterC[BG_FAST] / matlabTotalC,
		matlabLitterC[BG_SLOW] / matlabTotalC,
		matlabLitterC[BG_DEADMIC] / matlabTotalC
	};
	const double cFracProtected[BG_CLASS_COUNT] = {
		matlabProtectedC[BG_FAST] / matlabTotalC,
		matlabProtectedC[BG_SLOW] / matlabTotalC,
		matlabProtectedC[BG_DEADMIC] / matlabTotalC
	};
	const double nFracUnprotected[BG_CLASS_COUNT] = {
		matlabLitterN[BG_FAST] / matlabTotalN,
		matlabLitterN[BG_SLOW] / matlabTotalN,
		matlabLitterN[BG_DEADMIC] / matlabTotalN
	};
	const double nFracProtected[BG_CLASS_COUNT] = {
		matlabProtectedN[BG_FAST] / matlabTotalN,
		matlabProtectedN[BG_SLOW] / matlabTotalN,
		matlabProtectedN[BG_DEADMIC] / matlabTotalN
	};
	const double cFracMicrobe = matlabMicrobeC / matlabTotalC;
	const double nFracMicrobe = matlabMicrobeN / matlabTotalN;

	const double organicN0 = totalN0 * organicFrac;
	const double mineralN0 = std::max(0.0, totalN0 - organicN0);
	const double rhizoN0_scaled = mineralN0 * mineralRhizoFrac;
	const double bulkN0_scaled = std::max(0.0, mineralN0 - rhizoN0_scaled);

	auto assign_compartment = [&](FUN_CORPSE_Compartment* comp, double totalC, double totalN)
	{
		for (int i = 0; i < BG_CLASS_COUNT; ++i)
		{
			comp->unprotC[i] = totalC * cFracUnprotected[i];
			comp->protC[i] = totalC * cFracProtected[i];
			comp->unprotN[i] = totalN * nFracUnprotected[i];
			comp->protN[i] = totalN * nFracProtected[i];
		}
		comp->livingMicrobeC = totalC * cFracMicrobe;
		comp->livingMicrobeN = totalN * nFracMicrobe;
	};

	assign_compartment(org, organicC0, organicN0);
	assign_compartment(rhizo, rhizoC0, rhizoN0_scaled);
	assign_compartment(bulk, bulkC0, bulkN0_scaled);

	// Reassign the initial bulk slow pool between protected and unprotected
	// states so startup mineralization can be tuned without changing total bulk
	// slow C or N.
	const double initialProtectedFrac = clamp01(params->initial_protected_fraction);
	const double bulkSlowC0 = std::max(0.0, bulk->unprotC[BG_SLOW] + bulk->protC[BG_SLOW]);
	const double bulkSlowN0 = std::max(0.0, bulk->unprotN[BG_SLOW] + bulk->protN[BG_SLOW]);
	bulk->protC[BG_SLOW] = bulkSlowC0 * initialProtectedFrac;
	bulk->unprotC[BG_SLOW] = std::max(0.0, bulkSlowC0 - bulk->protC[BG_SLOW]);
	bulk->protN[BG_SLOW] = bulkSlowN0 * initialProtectedFrac;
	bulk->unprotN[BG_SLOW] = std::max(0.0, bulkSlowN0 - bulk->protN[BG_SLOW]);

	auto scale_active_pool_n = [&](FUN_CORPSE_Compartment* comp, double scalar)
	{
		const double nScale = std::max(0.0, scalar);
		comp->unprotN[BG_FAST] *= nScale;
		comp->unprotN[BG_DEADMIC] *= nScale;
	};

	scale_active_pool_n(org, params->initial_active_n_scalar_organic);
	scale_active_pool_n(rhizo, params->initial_active_n_scalar_mineral);
	scale_active_pool_n(bulk, params->initial_active_n_scalar_mineral);

	const double bulkSlowNScale = std::max(0.0, params->initial_bulk_slow_n_scalar);
	const double bulkDeadmicNScale = std::max(0.0, params->initial_bulk_deadmic_n_scalar);
	bulk->unprotN[BG_SLOW] *= bulkSlowNScale;
	bulk->protN[BG_SLOW] *= bulkSlowNScale;
	bulk->unprotN[BG_DEADMIC] *= bulkDeadmicNScale;
	bulk->protN[BG_DEADMIC] *= bulkDeadmicNScale;

	share->HOM = totalC0;
	share->HON = totalN0;

	SyncLegacySoilPoolsFromFUN(share);
}

double pnet_model::ComputeDynamicMicrobialCN(share_struct* share)
{
	FUN_CORPSE_Params* params = &share->funParams;
	if (params->use_dynamic_microbe_cn == 0)
	{
		return params->CN_microb_default;
	}

	const double Navail = std::max(0.0, share->NH4 + share->NO3);
	return params->CN_microb_min
		+ (params->CN_microb_max - params->CN_microb_min)
		* params->K_Navail_CN / (params->K_Navail_CN + Navail + 1.0e-12);
}

double pnet_model::ComputePlantNDeficit(veg_struct* veg, share_struct* share)
{
	return std::max(0.0, veg->MaxNStore - share->PlantN);
}

double pnet_model::ComputeAdditionalNacqCAllocation(veg_struct* veg, share_struct* share)
{
	double N_deficit = ComputePlantNDeficit(veg, share);
	if (N_deficit <= 0.0) return 0.0;

	double C_Nacq_max = share->funParams.C_Nacq_max_frac_GPP * std::max(0.0, share->GrsPsnMo);
	if (C_Nacq_max <= 0.0) return 0.0;

	double C_Nacq = C_Nacq_max * N_deficit / (share->funParams.K_Ndef + N_deficit);
	C_Nacq = std::min(C_Nacq, std::max(0.0, share->PlantC));
	return std::max(0.0, C_Nacq);
}

double pnet_model::GetFUNRhizoFractionForDOY(share_struct* share, int doy)
{
#ifdef FUN_DISABLE_RHIZO_SEASONALITY
	return clamp01(share->funParams.rhizo_fraction);
#else
	if (doy < 1) doy = 1;
	if (doy > 366) doy = 366;

	if (share->FUNUseRhizoSeasonality != 0)
	{
		return clamp01(share->FUNRhizoFractionByDOY[doy]);
	}

	return clamp01(share->funParams.rhizo_fraction);
#endif
}

void pnet_model::ApplyFUNRhizoSeasonality(share_struct* share, int doy)
{
#ifdef FUN_DISABLE_RHIZO_SEASONALITY
	share->FUNRhizoFractionCurrent = clamp01(share->funParams.rhizo_fraction);
	return;
#else
	FUN_CORPSE_Compartment* rhizo = &share->funBG.mineralRhizosphere;
	FUN_CORPSE_Compartment* bulk = &share->funBG.mineralBulk;
	const double currentRhizoFrac = GetFUNRhizoFractionForDOY(share, doy);

	if (doy >= 365)
	{
		for (int i = 0; i < BG_CLASS_COUNT; ++i)
		{
			double totalUnprotC = rhizo->unprotC[i] + bulk->unprotC[i];
			double totalUnprotN = rhizo->unprotN[i] + bulk->unprotN[i];
			double totalProtC = rhizo->protC[i] + bulk->protC[i];
			double totalProtN = rhizo->protN[i] + bulk->protN[i];

			rhizo->unprotC[i] = totalUnprotC * currentRhizoFrac;
			rhizo->unprotN[i] = totalUnprotN * currentRhizoFrac;
			rhizo->protC[i] = totalProtC * currentRhizoFrac;
			rhizo->protN[i] = totalProtN * currentRhizoFrac;

			bulk->unprotC[i] = totalUnprotC * (1.0 - currentRhizoFrac);
			bulk->unprotN[i] = totalUnprotN * (1.0 - currentRhizoFrac);
			bulk->protC[i] = totalProtC * (1.0 - currentRhizoFrac);
			bulk->protN[i] = totalProtN * (1.0 - currentRhizoFrac);
		}

		double totalMicrobeC = rhizo->livingMicrobeC + bulk->livingMicrobeC;
		double totalMicrobeN = rhizo->livingMicrobeN + bulk->livingMicrobeN;
		rhizo->livingMicrobeC = totalMicrobeC * currentRhizoFrac;
		rhizo->livingMicrobeN = totalMicrobeN * currentRhizoFrac;
		bulk->livingMicrobeC = totalMicrobeC * (1.0 - currentRhizoFrac);
		bulk->livingMicrobeN = totalMicrobeN * (1.0 - currentRhizoFrac);
		share->FUNRhizoFractionCurrent = currentRhizoFrac;
		return;
	}

	const double nextRhizoFrac = GetFUNRhizoFractionForDOY(share, doy + 1);
	const double rhizoDiff = nextRhizoFrac - currentRhizoFrac;

	if (rhizoDiff > 0.0)
	{
		for (int i = 0; i < BG_CLASS_COUNT; ++i)
		{
			double movedUnprotC = bulk->unprotC[i] * rhizoDiff;
			double movedUnprotN = bulk->unprotN[i] * rhizoDiff;
			double movedProtC = bulk->protC[i] * rhizoDiff;
			double movedProtN = bulk->protN[i] * rhizoDiff;

			bulk->unprotC[i] -= movedUnprotC;
			bulk->unprotN[i] -= movedUnprotN;
			bulk->protC[i] -= movedProtC;
			bulk->protN[i] -= movedProtN;

			rhizo->unprotC[i] += movedUnprotC;
			rhizo->unprotN[i] += movedUnprotN;
			rhizo->protC[i] += movedProtC;
			rhizo->protN[i] += movedProtN;
		}

		double movedMicrobeC = bulk->livingMicrobeC * rhizoDiff;
		double movedMicrobeN = bulk->livingMicrobeN * rhizoDiff;
		bulk->livingMicrobeC -= movedMicrobeC;
		bulk->livingMicrobeN -= movedMicrobeN;
		rhizo->livingMicrobeC += movedMicrobeC;
		rhizo->livingMicrobeN += movedMicrobeN;
	}
	else if (rhizoDiff < 0.0)
	{
		const double moveFrac = -rhizoDiff;
		for (int i = 0; i < BG_CLASS_COUNT; ++i)
		{
			double movedUnprotC = rhizo->unprotC[i] * moveFrac;
			double movedUnprotN = rhizo->unprotN[i] * moveFrac;
			double movedProtC = rhizo->protC[i] * moveFrac;
			double movedProtN = rhizo->protN[i] * moveFrac;

			rhizo->unprotC[i] -= movedUnprotC;
			rhizo->unprotN[i] -= movedUnprotN;
			rhizo->protC[i] -= movedProtC;
			rhizo->protN[i] -= movedProtN;

			bulk->unprotC[i] += movedUnprotC;
			bulk->unprotN[i] += movedUnprotN;
			bulk->protC[i] += movedProtC;
			bulk->protN[i] += movedProtN;
		}

		double movedMicrobeC = rhizo->livingMicrobeC * moveFrac;
		double movedMicrobeN = rhizo->livingMicrobeN * moveFrac;
		rhizo->livingMicrobeC -= movedMicrobeC;
		rhizo->livingMicrobeN -= movedMicrobeN;
		bulk->livingMicrobeC += movedMicrobeC;
		bulk->livingMicrobeN += movedMicrobeN;
	}

	share->FUNRhizoFractionCurrent = currentRhizoFrac;
#endif
}

double pnet_model::ComputeFUNPassiveNUptake(site_struct* site, veg_struct* veg, share_struct* share, double NdemandTotal)
{
	(void)NdemandTotal;

	double soilN = std::max(0.0, share->NH4 + share->NO3);
	const double plantCapacity = std::max(0.0, veg->MaxNStore - share->PlantN);
	if (soilN <= 0.0 || plantCapacity <= 0.0)
	{
		return 0.0;
	}

	if (share->RootMass < 0.9 * share->OldRoot && share->OldRoot > 0.0)
	{
		share->UptakeEff = share->RootMass / share->OldRoot;
	}
	else
	{
		share->UptakeEff = 1.0;
		share->OldRoot = 0.0;
	}

	double TMult = (exp(0.1 * (share->Tave - 7.1)) * 0.68);
	double RootNSinkStr = share->RootNSinkEff * TMult;
	RootNSinkStr = RootNSinkStr * (share->dayspan / 30.4) * 1.5;
	if (site->SnowPack > 0)
	{
		RootNSinkStr = 0.0;
	}
	else if (share->GDDTot > 50 && share->FolMass == 0)
	{
		RootNSinkStr *= 0.3;
	}

	double passiveN = soilN * RootNSinkStr * share->UptakeEff;
	const double soilCap = 0.8 * soilN;
	passiveN = std::min(passiveN, soilCap);
	passiveN = std::min(passiveN, plantCapacity);

	if (passiveN <= 0.0)
	{
		return 0.0;
	}

	double uptakeFrac = passiveN / std::max(soilN, 1.0e-12);
	double NH4Up = std::min(share->NH4, share->NH4 * uptakeFrac);
	double NO3Up = std::min(share->NO3, passiveN - NH4Up);

	share->NH4 -= NH4Up;
	share->NO3 -= NO3Up;
	share->PlantN += passiveN;
	share->PlantNUptakeYr += passiveN;
	share->PlantNUptakeLast = passiveN;

	return passiveN;
}

void pnet_model::ApplyFUNCostBasedNUptake(site_struct* site, veg_struct* veg, share_struct* share, int doy, double Cavailable, double NdemandTotal)
{
	FUN_CORPSE_Params* params = &share->funParams;
	const double dt_years = share->dayspan / 365.0;
	double passiveN = ComputeFUNPassiveNUptake(site, veg, share, NdemandTotal);
	double soilN = std::max(0.0, share->NH4 + share->NO3);
	double soilOrganicN = std::max(0.0, sum_fun_unprotected_n(share));
	// Keep the legacy PlantN pool as storage for AllocateMo rather than consuming
	// it here first; otherwise FUN acquisition is masked and root exudation stays off.
	double freeN = passiveN;
	double Ndeficit = std::max(0.0, NdemandTotal - freeN);
	double rootBiomass = std::max(1.0e-6, share->RootMass);
	double activeCost = 0.0;
	double nonmycCost = 0.0;
	double fixCost = std::numeric_limits<double>::infinity();
	double activeN = 0.0;
	double nonmycN = 0.0;
	double fixN = 0.0;
	double plantCN = 0.0;
	const double totalLivingMicrobeN =
		std::max(0.0, share->funBG.organicLayer.livingMicrobeN)
		+ std::max(0.0, share->funBG.mineralRhizosphere.livingMicrobeN)
		+ std::max(0.0, share->funBG.mineralBulk.livingMicrobeN);
	const double microbialFixCapacityN =
		totalLivingMicrobeN * dt_years / std::max(params->Tmic, 1.0e-12);

	Cavailable = std::max(0.0, Cavailable);
	if (Cavailable <= 0.0)
	{
		return;
	}
	if (Ndeficit <= 0.0)
	{
		share->PlantNUptakeLast = passiveN;
		return;
	}

	if (params->enable_root_exudation != 0 && params->f_exudate_Nacq > 0.0 && soilN > 0.0)
	{
		nonmycCost = params->nonmyc_cost_n / std::max(soilN, 1.0e-12)
			+ params->nonmyc_cost_c / rootBiomass;
	}

	if (params->enable_mycorrhizal_uptake != 0 && params->f_myc_Nacq > 0.0 && soilOrganicN > 0.0)
	{
		activeCost = params->myc_cost_n / std::max(soilOrganicN, 1.0e-12)
			+ params->myc_cost_c / rootBiomass;
	}

	fixCost = compute_fixation_cost(share->Tave);
	if (!(fixCost > 0.0))
	{
		fixCost = std::numeric_limits<double>::infinity();
	}

	double recCostAcq = 0.0;
	if (activeCost > 0.0)
	{
		recCostAcq += 1.0 / activeCost;
	}
	if (nonmycCost > 0.0)
	{
		recCostAcq += 1.0 / nonmycCost;
	}
	if (std::isfinite(fixCost))
	{
		recCostAcq += 1.0 / fixCost;
	}
	if (recCostAcq <= 0.0)
	{
		return;
	}

	double costAcq = 1.0 / recCostAcq;
	plantCN = Cavailable / std::max(NdemandTotal, 1.0e-12);
	plantCN = std::max(1.0, std::min(500.0, plantCN));

	double Cacq = (Cavailable - freeN * plantCN) / (1.0 + plantCN / std::max(costAcq, 1.0e-12));
	Cacq = std::max(0.0, Cacq);
	if (Cacq <= 0.0)
	{
		return;
	}

	if (activeCost > 0.0)
	{
		activeN = Cacq / activeCost;
		const double mycSupplyCapacityN =
			params->Vmax_myc_N * std::max(1.0, (double)share->dayspan)
			* (rootBiomass / std::max(params->K_myc_C + rootBiomass, 1.0e-12))
			* (soilOrganicN / std::max(params->K_myc_N + soilOrganicN, 1.0e-12));
		activeN = std::min(activeN, std::max(0.0, mycSupplyCapacityN));
		activeN = std::min(activeN, soilOrganicN);
	}
	if (nonmycCost > 0.0)
	{
		nonmycN = Cacq / nonmycCost;
	}
	if (std::isfinite(fixCost))
	{
		// Temporary diagnostic switch: disable biological N fixation so the
		// remaining ecosystem N imbalance can be assessed without this external input.
		fixN = 0.0;
	}

	double totalMineralRemoval = nonmycN;
	if (totalMineralRemoval > soilN && totalMineralRemoval > 0.0)
	{
		double scale = soilN / totalMineralRemoval;
		nonmycN *= scale;
		totalMineralRemoval = soilN;
	}

	double totalPotentialN = activeN + nonmycN + fixN;
	double Nacq = std::min(Ndeficit, totalPotentialN);
	if (Nacq <= 0.0 || totalPotentialN <= 0.0)
	{
		return;
	}
	if (Nacq < totalPotentialN)
	{
		double scale = Nacq / totalPotentialN;
		activeN *= scale;
		nonmycN *= scale;
		fixN *= scale;
		totalMineralRemoval = nonmycN;
	}

	activeN = remove_fun_unprotected_n(share, activeN);
	totalMineralRemoval = nonmycN;

	double Cspent = std::min(Cavailable, (activeN + nonmycN + fixN) * costAcq);
	double mycGrossN = activeN;
	double mycNToPlant = mycGrossN * clamp01(params->eta_myc_to_plant);
	double fungalRetainedN = std::max(0.0, mycGrossN - mycNToPlant);
	double nh4Frac = share->NH4 / std::max(soilN, 1.0e-12);
	double nh4Removal = std::min(share->NH4, totalMineralRemoval * nh4Frac);
	double no3Removal = std::min(share->NO3, totalMineralRemoval - nh4Removal);
	double rhizoCflux = (nonmycN + fixN) * costAcq;
	double fungalProdC = mycGrossN * costAcq;
	double rhizFungiFrac = GetFUNRhizoFractionForDOY(share, doy);
	double litterFungiFrac = clamp01(params->litter_fungi_frac);

	share->NH4 -= nh4Removal;
	share->NO3 -= no3Removal;
	share->PlantN += nonmycN + mycNToPlant + fixN;
	share->PlantNUptakeYr += nonmycN + mycNToPlant + fixN;
	share->PlantNUptakeLast = passiveN + nonmycN + mycNToPlant + fixN;
	if (share->PlantN > veg->MaxNStore)
	{
		share->FUNPlantNOverflowYr += (share->PlantN - veg->MaxNStore);
	}
	share->FUNFixNYr += fixN;
	share->PlantC -= Cspent;
	if (share->PlantC < 0.0) share->PlantC = 0.0;

	share->funBG.rootExudateC += rhizoCflux;
	share->funBG.mycorrhizalC += fungalProdC;
	share->funBG.mycorrhizalNToPlant += mycNToPlant;
	share->RootExudateCYr += rhizoCflux;
	share->MycorrhizalCYr += fungalProdC;
	share->MycorrhizalNToPlantYr += mycNToPlant;

	share->funBG.mineralRhizosphere.unprotC[BG_FAST] += share->funBG.rootExudateC;
	share->funBG.mineralRhizosphere.unprotC[BG_FAST] += fungalProdC * rhizFungiFrac * (1.0 - litterFungiFrac);
	share->funBG.mineralBulk.unprotC[BG_FAST] += fungalProdC * (1.0 - rhizFungiFrac) * (1.0 - litterFungiFrac);
	share->funBG.organicLayer.unprotC[BG_FAST] += fungalProdC * litterFungiFrac;

	if (fungalRetainedN > 0.0)
	{
		share->funBG.mineralRhizosphere.unprotN[BG_FAST] += fungalRetainedN * rhizFungiFrac;
		share->funBG.mineralBulk.unprotN[BG_FAST] += fungalRetainedN * (1.0 - rhizFungiFrac);
	}

	SyncLegacySoilPoolsFromFUN(share);
}

void pnet_model::AllocateMoFUN(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share, int CN_Mode)
{
	const double demandThisStep = std::max(0.0, share->FUNPotentialNdemand);
	const double potentialRootProdC = std::max(0.0, share->FUNPotentialRootProdC);
	const double potentialWoodProdC = std::max(0.0, share->FUNPotentialWoodProdC);
	const double totalN_before = compute_current_total_n(share);

	AllocateMo(veg, share, rstep, CN_Mode);
	const double totalN_after = compute_current_total_n(share);
	if (fabs(totalN_after - totalN_before) > 1.0e-8)
	{
		printf("FUN NDIAG AllocateMoFUN year %d doy %d dN=%.10f before=%.10f after=%.10f\n",
			clim->year[rstep], clim->doy[rstep], totalN_after - totalN_before, totalN_before, totalN_after);
	}

	share->FUNNdemandYr += demandThisStep;
	share->FUNNdemandGapYr += std::max(0.0, demandThisStep - std::max(0.0, share->PlantNUptakeLast));
	share->FUNPotentialRootProdCYr += potentialRootProdC;
	share->FUNPotentialWoodProdCYr += potentialWoodProdC;
	share->FUNNLimitRootProdCYr += std::max(0.0, potentialRootProdC - std::max(0.0, share->RootProdCMo));
	share->FUNNLimitWoodProdCYr += std::max(0.0, potentialWoodProdC - std::max(0.0, share->WoodProdCMo));

	share->FUNPotentialGrowthC = 0.0;
	share->FUNPotentialNdemand = 0.0;
	share->FUNPotentialRootProdC = 0.0;
	share->FUNPotentialWoodProdC = 0.0;
}

void pnet_model::TransferOrganicLayerToMineralSoil(share_struct* share)
{
	FUN_CORPSE_Compartment* org = &share->funBG.organicLayer;
	FUN_CORPSE_Compartment* rhizo = &share->funBG.mineralRhizosphere;
	FUN_CORPSE_Compartment* bulk = &share->funBG.mineralBulk;
	const double transferFrac = clamp01(share->funParams.litter_transfer_to_soil);
	const double rhizoFrac = share->FUNRhizoFractionCurrent;

	if (transferFrac <= 0.0) return;

	for (int i = 0; i < BG_CLASS_COUNT; ++i)
	{
		double movedUnprotC = org->unprotC[i] * transferFrac;
		double movedUnprotN = org->unprotN[i] * transferFrac;
		double movedProtC = org->protC[i] * transferFrac;
		double movedProtN = org->protN[i] * transferFrac;

		org->unprotC[i] -= movedUnprotC;
		org->unprotN[i] -= movedUnprotN;
		org->protC[i] -= movedProtC;
		org->protN[i] -= movedProtN;

		rhizo->unprotC[i] += movedUnprotC * rhizoFrac;
		rhizo->unprotN[i] += movedUnprotN * rhizoFrac;
		rhizo->protC[i] += movedProtC * rhizoFrac;
		rhizo->protN[i] += movedProtN * rhizoFrac;

		bulk->unprotC[i] += movedUnprotC * (1.0 - rhizoFrac);
		bulk->unprotN[i] += movedUnprotN * (1.0 - rhizoFrac);
		bulk->protC[i] += movedProtC * (1.0 - rhizoFrac);
		bulk->protN[i] += movedProtN * (1.0 - rhizoFrac);
	}
}

void pnet_model::SyncLegacySoilPoolsFromFUN(share_struct* share)
{
	share->HOM = sum_compartment_c(&share->funBG.organicLayer)
		+ sum_compartment_c(&share->funBG.mineralRhizosphere)
		+ sum_compartment_c(&share->funBG.mineralBulk);
	share->HON = sum_compartment_n(&share->funBG.organicLayer)
		+ sum_compartment_n(&share->funBG.mineralRhizosphere)
		+ sum_compartment_n(&share->funBG.mineralBulk);
}

void pnet_model::UpdateCORPSECompartment(const char* label, FUN_CORPSE_Compartment* comp, double soilT_K,
	double theta, double* NH4, double* NO3, double dt_years, FUN_CORPSE_Params* params,
	double microbialCN, double decompScalar)
{
	const double gasConst = 8.314;
	const double refTemp = 293.15;
	const double thetaEff = clamp01(theta);
	const double unprotN_before = comp->unprotN[BG_FAST] + comp->unprotN[BG_SLOW] + comp->unprotN[BG_DEADMIC];
	const double protN_before = comp->protN[BG_FAST] + comp->protN[BG_SLOW] + comp->protN[BG_DEADMIC];
	const double microbeN_before = comp->livingMicrobeN;
	const double mineralN_before = std::max(0.0, *NH4) + std::max(0.0, *NO3);
	double potentialResp[BG_CLASS_COUNT] = { 0.0, 0.0, 0.0 };
	double decomposedN[BG_CLASS_COUNT] = { 0.0, 0.0, 0.0 };
	double carbon_supply = 0.0;
	double nitrogen_supply = 0.0;
	double totalUnprotC = 0.0;

	for (int i = 0; i < BG_CLASS_COUNT; ++i)
	{
		comp->lastDecompC[i] = 0.0;
		comp->lastDecompN[i] = 0.0;
		comp->lastProtectionN = 0.0;
		comp->lastDeprotectionN = 0.0;
		comp->lastTurnoverN = 0.0;
		totalUnprotC += comp->unprotC[i];
	}

	const double targetMicrobeC = std::max(0.01, params->minMicrobeC * totalUnprotC);
	if (comp->livingMicrobeC < targetMicrobeC && totalUnprotC > 0.0)
	{
		double recruitC = targetMicrobeC - std::max(0.0, comp->livingMicrobeC);
		double recruitedN = 0.0;
		for (int i = 0; i < BG_CLASS_COUNT; ++i)
		{
			if (comp->unprotC[i] <= 0.0) continue;

			const double poolShare = comp->unprotC[i] / std::max(totalUnprotC, 1.0e-12);
			const double movedC = std::min(comp->unprotC[i], recruitC * poolShare);
			double movedN = 0.0;
			if (comp->unprotC[i] > 0.0 && comp->unprotN[i] > 0.0)
			{
				movedN = std::min(comp->unprotN[i], movedC * comp->unprotN[i] / std::max(comp->unprotC[i], 1.0e-12));
			}

			comp->unprotC[i] -= movedC;
			comp->unprotN[i] -= movedN;
			comp->livingMicrobeC += movedC;
			recruitedN += movedN;
		}

		comp->livingMicrobeN += recruitedN;
	}

	for (int i = 0; i < BG_CLASS_COUNT; ++i)
	{
		if (comp->unprotC[i] <= 0.0) continue;

		double tempScalar = exp((-params->Ea[i] / gasConst) * (1.0 / soilT_K - 1.0 / refTemp));
		double waterScalar = pow(thetaEff, 3.0) * pow(std::max(0.0, 1.0 - thetaEff), params->gas_diffusion_exp);
		double heterotrophC = std::max(0.0, comp->livingMicrobeC);
		double substrateTerm = 0.0;
		double denom = comp->unprotC[i] * params->kC[i] + heterotrophC;
		if (denom > 0.0)
		{
			substrateTerm = comp->unprotC[i] * heterotrophC / denom;
		}
		potentialResp[i] = params->vmaxref[i] * decompScalar * tempScalar * waterScalar * substrateTerm;
		potentialResp[i] = std::min(potentialResp[i], comp->unprotC[i] / std::max(dt_years, 1.0e-12));

		if (comp->unprotN[i] > 0.0)
		{
			decomposedN[i] = potentialResp[i] * comp->unprotN[i] / std::max(comp->unprotC[i], 1.0e-12);
		}

		carbon_supply += potentialResp[i] * params->eup[i];
		nitrogen_supply += decomposedN[i] * params->nup[i];
	}

	double microbeTurnover = std::max(0.0,
		(comp->livingMicrobeC - params->minMicrobeC * std::max(1.0, totalUnprotC)) / std::max(params->Tmic, 1.0e-12));
	double maintenance_resp = microbeTurnover * (1.0 - params->et);
	double availableGrowthC = std::max(0.0, carbon_supply - maintenance_resp);
	double mineralN_total = std::max(0.0, *NH4 + *NO3);
	double IMM_N_max = mineralN_total * params->max_immobilization_rate * 365.0;
	IMM_N_max = std::min(IMM_N_max, mineralN_total / std::max(dt_years, 1.0e-12));

	double N_required_for_growth = availableGrowthC / std::max(1.0, microbialCN);
	double grossImmobilization = 0.0;
	double grossMineralizationFromStoich = 0.0;
	double microbialGrowthC = 0.0;
	double overflowResp = 0.0;
	double CN_imbalance_term = 0.0;

	if (availableGrowthC > (nitrogen_supply + IMM_N_max) * microbialCN)
	{
		grossImmobilization = IMM_N_max;
		microbialGrowthC = (nitrogen_supply + IMM_N_max) * microbialCN;
		overflowResp = availableGrowthC - microbialGrowthC;
		CN_imbalance_term = -IMM_N_max;
	}
	else if (availableGrowthC >= nitrogen_supply * microbialCN)
	{
		grossImmobilization = std::max(0.0, N_required_for_growth - nitrogen_supply);
		microbialGrowthC = availableGrowthC;
		CN_imbalance_term = -grossImmobilization;
	}
	else
	{
		grossMineralizationFromStoich = std::max(0.0, nitrogen_supply - N_required_for_growth);
		microbialGrowthC = availableGrowthC;
		CN_imbalance_term = grossMineralizationFromStoich;
	}

	double deadmic_C_produced = dt_years * microbeTurnover * params->et;
	const double microbeNPerC = std::max(0.0, comp->livingMicrobeN) / std::max(comp->livingMicrobeC, 1.0e-12);
	double total_turnover_N = dt_years * microbeTurnover * microbeNPerC;
	double structural_turnover_N = total_turnover_N * params->et;
	double respiratory_turnover_N = std::max(0.0, total_turnover_N - structural_turnover_N);
	double turnover_N_min = structural_turnover_N * params->frac_turnover_min;
	double turnover_N_slow = structural_turnover_N * params->frac_turnover_slow;
	double deadmic_N_produced = std::max(0.0, structural_turnover_N - turnover_N_min - turnover_N_slow);

	for (int i = 0; i < BG_CLASS_COUNT; ++i)
	{
		double dC = dt_years * potentialResp[i];
		double dN = dt_years * decomposedN[i];
		comp->unprotC[i] = std::max(0.0, comp->unprotC[i] - dC);
		comp->unprotN[i] = std::max(0.0, comp->unprotN[i] - dN);
		comp->lastDecompC[i] = dC;
		comp->lastDecompN[i] = dN;
	}

	comp->unprotC[BG_DEADMIC] += deadmic_C_produced * (1.0 - params->frac_turnover_slow);
	comp->unprotC[BG_SLOW] += deadmic_C_produced * params->frac_turnover_slow;
	comp->unprotN[BG_DEADMIC] += deadmic_N_produced;
	comp->unprotN[BG_SLOW] += turnover_N_slow;

	double grossMineralization = 0.0;
	for (int i = 0; i < BG_CLASS_COUNT; ++i)
	{
		grossMineralization += (1.0 - params->nup[i]) * decomposedN[i];
	}
	grossMineralization += std::max(0.0, CN_imbalance_term)
		+ (turnover_N_min + respiratory_turnover_N) / std::max(dt_years, 1.0e-12);
	const double totalOrganicN_before = std::max(0.0, unprotN_before + protN_before + microbeN_before);
	const double grossNMinCap = std::max(0.0, totalOrganicN_before * params->gross_nmin_cap_frac_yr);
	if (grossNMinCap > 0.0 && grossMineralization > grossNMinCap)
	{
		const double excessGrossMineralization = grossMineralization - grossNMinCap;
		grossMineralization = grossNMinCap;
		comp->unprotN[BG_SLOW] += excessGrossMineralization * dt_years;
	}
	grossImmobilization = std::max(0.0, -CN_imbalance_term);

	*NH4 += grossMineralization * dt_years;
	double immobilizedN = grossImmobilization * dt_years;
	double fromNH4 = std::min(*NH4, immobilizedN);
	*NH4 -= fromNH4;
	double remaining = immobilizedN - fromNH4;
	double fromNO3 = std::min(*NO3, remaining);
	*NO3 -= fromNO3;
	double actualImmobilizedN = fromNH4 + fromNO3;
	double immobilizationShortfallN = std::max(0.0, immobilizedN - actualImmobilizedN);

	comp->livingMicrobeC = std::max(0.0, comp->livingMicrobeC + dt_years * (microbialGrowthC - microbeTurnover));
	const double microbialGrowthN_actual =
		dt_years * nitrogen_supply
		+ actualImmobilizedN
		- dt_years * std::max(0.0, grossMineralizationFromStoich);
	comp->livingMicrobeN = std::max(0.0, comp->livingMicrobeN + microbialGrowthN_actual - total_turnover_N);

	for (int i = 0; i < BG_CLASS_COUNT; ++i)
	{
		double protectedCturnover = comp->protC[i] / std::max(params->tProtected, 1.0e-12);
		double protectedNturnover = comp->protN[i] / std::max(params->tProtected, 1.0e-12);
		double movedC = dt_years * protectedCturnover;
		double movedN = dt_years * protectedNturnover;

		comp->protC[i] = std::max(0.0, comp->protC[i] - movedC);
		comp->protN[i] = std::max(0.0, comp->protN[i] - movedN);
		comp->unprotC[i] += movedC;
		comp->unprotN[i] += movedN;
		comp->lastDeprotectionN += movedN;

		double newProtectedC = params->protection_rate[i] * comp->unprotC[i] * dt_years;
		newProtectedC = std::min(newProtectedC, comp->unprotC[i]);
		double newProtectedN = 0.0;
		if (comp->unprotC[i] > 0.0)
		{
			newProtectedN = newProtectedC * comp->unprotN[i] / comp->unprotC[i];
		}

		comp->unprotC[i] -= newProtectedC;
		comp->unprotN[i] -= newProtectedN;
		comp->protC[i] += newProtectedC;
		comp->protN[i] += newProtectedN;
		comp->lastProtectionN += newProtectedN;
	}

	comp->CO2 = 0.0;
	for (int i = 0; i < BG_CLASS_COUNT; ++i)
	{
		comp->CO2 += dt_years * potentialResp[i] * std::max(0.0, 1.0 - params->eup[i]);
	}
	comp->CO2 += dt_years * (overflowResp + maintenance_resp);
	comp->Rtot = comp->CO2;
	comp->lastGrossMineralization = grossMineralization * dt_years;
	comp->lastGrossImmobilization = grossImmobilization * dt_years;
	comp->lastTurnoverN = total_turnover_N;
	const double unprotN_after = comp->unprotN[BG_FAST] + comp->unprotN[BG_SLOW] + comp->unprotN[BG_DEADMIC];
	const double protN_after = comp->protN[BG_FAST] + comp->protN[BG_SLOW] + comp->protN[BG_DEADMIC];
	const double microbeN_after = comp->livingMicrobeN;
	const double mineralN_after = std::max(0.0, *NH4) + std::max(0.0, *NO3);
	const double total_before = unprotN_before + protN_before + microbeN_before + mineralN_before;
	const double total_after = unprotN_after + protN_after + microbeN_after + mineralN_after;
	const double residualN = total_after - total_before;
	if (fabs(residualN) > 1.0e-8)
	{
		printf("FUN NDIAG Residual %s dN=%.10f dunprot=%.10f dprot=%.10f dmicrobe=%.10f dmineral=%.10f grossMin=%.10f grossImmob=%.10f turnover=%.10f minTurn=%.10f slowTurn=%.10f deadmic=%.10f\n",
			label,
			residualN,
			unprotN_after - unprotN_before,
			protN_after - protN_before,
			microbeN_after - microbeN_before,
			mineralN_after - mineralN_before,
			comp->lastGrossMineralization,
			comp->lastGrossImmobilization,
			total_turnover_N,
			turnover_N_min,
			turnover_N_slow,
			deadmic_N_produced);
	}
	if (immobilizationShortfallN > 1.0e-8)
	{
		printf("FUN NDIAG ImmobShortfall shortfall=%.10f target=%.10f actual=%.10f mineral_after_gmin=%.10f microbialCN=%.10f\n",
			immobilizationShortfallN, immobilizedN, actualImmobilizedN, std::max(0.0, *NH4 + *NO3), microbialCN);
	}
}

void pnet_model::UpdateFUN_CORPSEBelowground(site_struct* site, veg_struct* veg, clim_struct* clim,
	int rstep, share_struct* share, double TMult, double WMult)
{
	const double dt_years = share->dayspan / 365.0;
	const double dt_days = std::max(1.0, (double)share->dayspan);
	const double litterC_per_day = veg->CFracBiomass / dt_days;
	const double litterN_per_day = 1.0 / dt_days;
	const double theta = clamp01(share->Water / std::max(1.0, site->WHC));
	const double soilT_K = share->Tave + 273.15;
	FUN_CORPSE_Params* params = &share->funParams;

	ResetFUN_CORPSEFluxDiagnostics(share);
	share->funBG.microbialCN = ComputeDynamicMicrobialCN(share);
	share->BGMicrobialCNLast = share->funBG.microbialCN;

	add_split_to_pools(share->FolLitM * litterC_per_day, share->FolLitN * litterN_per_day, params->leaf_fast_frac, params->leaf_slow_frac,
		share->funBG.organicLayer.unprotC, share->funBG.organicLayer.unprotN);
	add_split_to_pools(share->WoodLitM * litterC_per_day, share->WoodLitN * litterN_per_day, params->deadwood_fast_frac, params->deadwood_slow_frac,
		share->funBG.organicLayer.unprotC, share->funBG.organicLayer.unprotN);
	const double rootLitC = share->RootLitM * litterC_per_day;
	const double rootRhizoFrac = (params->root_inputs_only_rhiz != 0) ? 1.0 : share->FUNRhizoFractionCurrent;
	add_split_to_pools(rootLitC * rootRhizoFrac, share->RootLitN * litterN_per_day * rootRhizoFrac, params->root_fast_frac, params->root_slow_frac,
		share->funBG.mineralRhizosphere.unprotC, share->funBG.mineralRhizosphere.unprotN);
	add_split_to_pools(rootLitC * (1.0 - rootRhizoFrac), share->RootLitN * litterN_per_day * (1.0 - rootRhizoFrac), params->root_fast_frac, params->root_slow_frac,
		share->funBG.mineralBulk.unprotC, share->funBG.mineralBulk.unprotN);

	double organic_n_before = sum_compartment_n(&share->funBG.organicLayer) + std::max(0.0, share->NH4) + std::max(0.0, share->NO3);
	UpdateCORPSECompartment("organic", &share->funBG.organicLayer, soilT_K, theta, &share->NH4, &share->NO3, dt_years, params, share->funBG.microbialCN,
		params->decomp_scalar_global * params->decomp_scalar_organic);
	double organic_n_after = sum_compartment_n(&share->funBG.organicLayer) + std::max(0.0, share->NH4) + std::max(0.0, share->NO3);
	log_fun_compartment_n_diag("organic", clim->year[rstep], clim->doy[rstep], organic_n_before, organic_n_after);
	ApplyFUNRhizoSeasonality(share, clim->doy[rstep]);
	TransferOrganicLayerToMineralSoil(share);
	double rhizo_n_before = sum_compartment_n(&share->funBG.mineralRhizosphere) + std::max(0.0, share->NH4) + std::max(0.0, share->NO3);
	UpdateCORPSECompartment("rhizo", &share->funBG.mineralRhizosphere, soilT_K, theta, &share->NH4, &share->NO3, dt_years, params, share->funBG.microbialCN,
		params->decomp_scalar_global * params->decomp_scalar_rhizo);
	double rhizo_n_after = sum_compartment_n(&share->funBG.mineralRhizosphere) + std::max(0.0, share->NH4) + std::max(0.0, share->NO3);
	log_fun_compartment_n_diag("rhizo", clim->year[rstep], clim->doy[rstep], rhizo_n_before, rhizo_n_after);
	double bulk_n_before = sum_compartment_n(&share->funBG.mineralBulk) + std::max(0.0, share->NH4) + std::max(0.0, share->NO3);
	UpdateCORPSECompartment("bulk", &share->funBG.mineralBulk, soilT_K, theta, &share->NH4, &share->NO3, dt_years, params, share->funBG.microbialCN,
		params->decomp_scalar_global * params->decomp_scalar_bulk);
	double bulk_n_after = sum_compartment_n(&share->funBG.mineralBulk) + std::max(0.0, share->NH4) + std::max(0.0, share->NO3);
	log_fun_compartment_n_diag("bulk", clim->year[rstep], clim->doy[rstep], bulk_n_before, bulk_n_after);

	FUN_CORPSE_Compartment* comps[3] = {
		&share->funBG.organicLayer,
		&share->funBG.mineralRhizosphere,
		&share->funBG.mineralBulk
	};

	double soilResp = 0.0;
	for (int i = 0; i < 3; ++i)
	{
		soilResp += comps[i]->CO2;
		share->funBG.grossNMobilization += comps[i]->lastGrossMineralization;
		share->funBG.grossNImmobilization += comps[i]->lastGrossImmobilization;
		share->funBG.protectionN += comps[i]->lastProtectionN;
		share->funBG.deprotectionN += comps[i]->lastDeprotectionN;
		share->funBG.microbialTurnoverN += comps[i]->lastTurnoverN;
	}

	share->BGOrgGrossNMineralizationYr += share->funBG.organicLayer.lastGrossMineralization;
	share->BGRhizGrossNMineralizationYr += share->funBG.mineralRhizosphere.lastGrossMineralization;
	share->BGBulkGrossNMineralizationYr += share->funBG.mineralBulk.lastGrossMineralization;
	share->BGOrgGrossNImmobilizationYr += share->funBG.organicLayer.lastGrossImmobilization;
	share->BGRhizGrossNImmobilizationYr += share->funBG.mineralRhizosphere.lastGrossImmobilization;
	share->BGBulkGrossNImmobilizationYr += share->funBG.mineralBulk.lastGrossImmobilization;
	share->BGOrgNetNMineralizationYr += share->funBG.organicLayer.lastGrossMineralization - share->funBG.organicLayer.lastGrossImmobilization;
	share->BGRhizNetNMineralizationYr += share->funBG.mineralRhizosphere.lastGrossMineralization - share->funBG.mineralRhizosphere.lastGrossImmobilization;
	share->BGBulkNetNMineralizationYr += share->funBG.mineralBulk.lastGrossMineralization - share->funBG.mineralBulk.lastGrossImmobilization;

	share->funBG.netNMineralization = share->funBG.grossNMobilization - share->funBG.grossNImmobilization;
	share->SoilDecResp = soilResp;
	share->SoilDecRespYr += soilResp;
	share->GrossNMinYr += share->funBG.grossNMobilization;
	share->GrossNImmobYr += share->funBG.grossNImmobilization;
	share->NetNMinYr += share->funBG.netNMineralization;
	share->NetCBal -= soilResp;

	share->BGGrossNMobilizationYr += share->funBG.grossNMobilization;
	share->BGGrossNImmobilizationYr += share->funBG.grossNImmobilization;
	share->BGNetNMineralizationYr += share->funBG.netNMineralization;
	share->BGMicrobialTurnoverNYr += share->funBG.microbialTurnoverN;
	share->BGProtectionNYr += share->funBG.protectionN;
	share->BGDeprotectionNYr += share->funBG.deprotectionN;

	SyncLegacySoilPoolsFromFUN(share);
}

void pnet_model::PnETActiveRootUptake(site_struct* site, veg_struct* veg, share_struct* share, double TMult)
{
	double RootNSinkStr = share->RootNSinkEff * TMult;
	RootNSinkStr = RootNSinkStr * (share->dayspan / 30.4) * 1.5;

	if (site->SnowPack > 0)
	{
		RootNSinkStr = 0.0;
	}
	else if (share->GDDTot > 50 && share->FolMass == 0)
	{
		RootNSinkStr *= 0.3;
	}

	if (share->RootMass < 0.9 * share->OldRoot && share->OldRoot > 0.0)
	{
		share->UptakeEff = share->RootMass / share->OldRoot;
	}
	else
	{
		share->UptakeEff = 1.0;
		share->OldRoot = 0.0;
	}

	double PlantNUptake = (share->NH4 + share->NO3) * RootNSinkStr * share->UptakeEff;
	if ((PlantNUptake + share->PlantN) > veg->MaxNStore)
	{
		PlantNUptake = veg->MaxNStore - share->PlantN;
		if ((share->NO3 + share->NH4) > 0.0)
		{
			RootNSinkStr = PlantNUptake / (share->NO3 + share->NH4);
		}
	}

	if (PlantNUptake < 0.0)
	{
		PlantNUptake = 0.0;
		RootNSinkStr = 0.0;
	}

	share->PlantN += PlantNUptake;
	share->PlantNUptakeYr += PlantNUptake;
	share->PlantNUptakeLast = PlantNUptake;

	double NH4Up = share->NH4 * RootNSinkStr;
	double NO3Up = share->NO3 * RootNSinkStr;
	share->NH4 -= NH4Up;
	share->NO3 -= NO3Up;
}

void pnet_model::DecompFUN_CORPSE(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share)
{
	double TMult, WMult, NetNitr, fNRatioNit;
	double FN2O, FNO, FN2, FNOd, FN2Od, NitLostGas, Kn, Fdeni;
	double R_NOx_N2O, Wnit, den_T, Q10, den_W, Wcrit, R_N2_N2O, den_CO2, Ccrit, Cmax;
	int iday;
	bool appliedFUNUptake = false;

	FN2O = 0.0;
	FNO = 0.0;
	FN2 = 0.0;
	FN2Od = 0.0;
	Kn = 0.01;

	for (iday = 1; iday <= share->dayspan; iday++)
	{
		if (share->SnowStatus == 1)
		{
			share->SnowNO3 += clim->NO3dep[rstep] / share->dayspan;
			share->SnowNH4 += clim->NH4dep[rstep] / share->dayspan;
		}
		else
		{
			share->NO3 += clim->NO3dep[rstep] / share->dayspan;
			share->NH4 += clim->NH4dep[rstep] / share->dayspan;
		}

		share->NdepTot += clim->NO3dep[rstep] / share->dayspan + clim->NH4dep[rstep] / share->dayspan;

		TMult = (exp(0.1 * (share->Tave - 7.1)) * 0.68);
		WMult = share->MeanSoilMoistEff;
		if (site->WaterStress == 0) WMult = 1.0;

		const double totalN_before_updatefun = compute_current_total_n(share);
		UpdateFUN_CORPSEBelowground(site, veg, clim, rstep, share, TMult, WMult);
		const double totalN_after_updatefun = compute_current_total_n(share);
		if (fabs(totalN_after_updatefun - totalN_before_updatefun) > 1.0e-8)
		{
			printf("FUN NDIAG UpdateFUN year %d doy %d dN=%.10f before=%.10f after=%.10f\n",
				clim->year[rstep], clim->doy[rstep], totalN_after_updatefun - totalN_before_updatefun,
				totalN_before_updatefun, totalN_after_updatefun);
		}
		if (!appliedFUNUptake && share->FUNPotentialGrowthC > 0.0 && share->FUNPotentialNdemand > 0.0)
		{
			const double totalN_before_fun = compute_current_total_n(share);
			ApplyFUNCostBasedNUptake(site, veg, share, clim->doy[rstep], share->FUNPotentialGrowthC, share->FUNPotentialNdemand);
			const double totalN_after_fun = compute_current_total_n(share);
			if (fabs(totalN_after_fun - totalN_before_fun) > 1.0e-8)
			{
				printf("FUN NDIAG ApplyFUN year %d doy %d dN=%.10f before=%.10f after=%.10f\n",
					clim->year[rstep], clim->doy[rstep], totalN_after_fun - totalN_before_fun, totalN_before_fun, totalN_after_fun);
			}
			appliedFUNUptake = true;
		}

		double nit_season_factor;
		int doy = clim->doy[rstep];
		if (doy <= 90 || doy >= 335)
			nit_season_factor = 0.462;
		else if (doy <= 181)
			nit_season_factor = 0.856;
		else if (doy <= 273)
			nit_season_factor = 1.000;
		else
			nit_season_factor = 0.819;

		fNRatioNit = share->dayspan / 30.4;
		const double NH4BeforeNitr = std::max(0.0, share->NH4);
		NetNitr = share->NH4 * share->NRatioNit * fNRatioNit * nit_season_factor;
		const double nitrificationCap = NH4BeforeNitr * std::max(0.0, share->funParams.nitrification_cap_frac_day) * share->dayspan;
		if (nitrificationCap > 0.0 && NetNitr > nitrificationCap)
		{
			NetNitr = nitrificationCap;
		}

		share->WFPS = share->Water / site->WHC * 0.75;
		if (share->WFPS > 1.0) share->WFPS = 1.0;

		Wnit = 1.0;
		NitLostGas = Kn * NetNitr * Wnit;
		NitLostGas = 0.0;

		R_NOx_N2O = pow(10.0, -3.79 * share->WFPS + 2.73);
		FN2O = NitLostGas / (1.0 + R_NOx_N2O);
		FNO = NitLostGas - FN2O;

		share->NO3 += NetNitr;
		share->NH4 -= NetNitr + NitLostGas;

		Q10 = 2.0;
		share->Tsoil = share->Tave;
		if (0.0 < share->Tsoil) den_T = pow(Q10, 0.1 * (share->Tsoil - 20.0));
		else den_T = 0.01;

		Wcrit = 0.75;
		den_W = 1.0 - 1.0 / (1.0 + exp(12.5 * (share->WFPS - Wcrit)));

		Ccrit = 1.8;
		Cmax = 3.0;
		den_CO2 = 1.0 - 1.0 / (1.0 + exp(6.8 * (share->SoilDecResp - Ccrit) / Cmax));
		den_CO2 = 1.0;

		site->Kden = 0.03;
		if (den_W > den_CO2) den_W = den_CO2;
		Fdeni = share->NO3 * site->Kden * den_T * den_W;

		double RnA, RnB;
		share->RnMax = 10.0;
		share->RnX1 = 0.7;
		share->RnY1 = 1.5;
		share->RnX2 = 1.0;
		share->RnY2 = share->RnMax;

		RnB = log(share->RnY2 / share->RnY1) / (share->RnX2 - share->RnX1);
		RnA = share->RnY1 * exp(-RnB * share->RnX1);
		R_N2_N2O = RnA * exp(RnB * share->WFPS);

		FN2Od = Fdeni / (1.0 + R_NOx_N2O + R_N2_N2O);
		FN2Od = 0.0;
		FNOd = R_NOx_N2O * FN2Od;
		FN2 = FN2Od * R_N2_N2O;

		share->FluxN2OYrDe += FN2Od;
		share->FluxNOYrDe += FNOd;
		share->FluxN2OYr += FN2O + FN2Od;
		share->FluxNOYr += FNO + FNOd;
		share->FluxN2Yr += FN2;

		share->NO3 -= FN2Od + FN2 + FNOd;
		share->NetNitrYr += NetNitr;
	}
}

void pnet_model::pnet_fun_corpse()
{
	char inputname[350];
	char climdayfile[350];
	int CN_Mode = 1;
	int rstep, ndays;
	int ystep = 1;
	int doyEnd;
	int YrStartSpin, YrEndSpin;

	veg_struct veg;
	site_struct site;
	clim_struct clim;
	out_struct out;
	share_struct share;

	clim_struct climmon, climday, climhr, climdayspin;

	FILE* fileoutD;
	sprintf(inputname, "%s%s%s", exePath, PathInput, "input.txt");
	sprintf(climdayfile, "%s%s%s", exePath, PathInput, "climateday.clim");

	clim.timestep = 0;
	climday.timestep = 1;
	climdayspin.timestep = 1;

	ReadInput(&site, &veg, &share, inputname);
	ReadFUNCORPSEParams(&share, inputname);
	ReadFUNRhizoSeasonality(&share, inputname);
	ReadClimDay(&climday, climdayfile);
	initvars(&site, &veg, &share);

	YrStartSpin = 1001;
	YrEndSpin = 2023;

	share.yrspin = YrEndSpin - YrStartSpin + 1;
	climday.length = share.yrspin * 365;

	memset_out(share.yrspin, &out);

	share.ifO3EffectOnPSN = 1;
	share.kO3EffScalar = 1.12;
	share.kO3Eff = share.kO3Eff * share.kO3EffScalar;
	share.WUEO3Eff = 0;

	Disturbance(&site, &veg, &share, 0);
	Disturbance(&site, &veg, &share, 1);

	printf("   ===============================================\n\n");
	printf("             PnET-FUN-CORPSE model starts daily run\n");

	ndays = climday.length;
	share.ISA = 0.00;

	for (rstep = 1; rstep <= ndays; rstep++)
	{
		AtmEnviron(&site, &climday, rstep, &share);
		Phenology(&veg, &climday, rstep, &share, 1);
		Photosyn(&site, &veg, &climday, rstep, &share);
		Waterbal(&site, &veg, &climday, rstep, &share);
		EstimateFUNAllocationDemand(&veg, &climday, rstep, &share, CN_Mode);
		DecompFUN_CORPSE(&site, &veg, &climday, rstep, &share);
		AllocateMoFUN(&site, &veg, &climday, rstep, &share, CN_Mode);
		{
			const double totalN_before_cntrans = compute_current_total_n(&share);
			Phenology(&veg, &climday, rstep, &share, 2);
			const double totalN_after_phen2 = compute_current_total_n(&share);
			if (fabs(totalN_after_phen2 - totalN_before_cntrans) > 1.0e-8)
			{
				printf("FUN NDIAG Phenology2 year %d doy %d dN=%.10f before=%.10f after=%.10f\n",
					climday.year[rstep], climday.doy[rstep], totalN_after_phen2 - totalN_before_cntrans,
					totalN_before_cntrans, totalN_after_phen2);
			}
			const double totalN_before_cn = totalN_after_phen2;
			CNTrans(&site, &veg, &climday, rstep, &share);
			const double totalN_after_cn = compute_current_total_n(&share);
			if (fabs(totalN_after_cn - totalN_before_cn) > 1.0e-8)
			{
				printf("FUN NDIAG CNTrans year %d doy %d dN=%.10f before=%.10f after=%.10f\n",
					climday.year[rstep], climday.doy[rstep], totalN_after_cn - totalN_before_cn,
					totalN_before_cn, totalN_after_cn);
			}
		}
		Leach(&climday, rstep, &share);

		if (climday.year[rstep] >= site.nYearStart && climday.year[rstep] <= site.nYearEnd)
			WriteoutDay(&site, &veg, &climday, &share, rstep, fileoutD);

		if (climday.doy[rstep] >= 365)
		{
			if (is_leapyear(climday.year[rstep]))
			{
				doyEnd = 366;
				doyEnd = 365;
			}
			else doyEnd = 365;

			if (climday.doy[rstep] == doyEnd)
			{
				int annualRstep = (rstep < ndays) ? (rstep + 1) : rstep;
				const double totalN_before_allocateyr = compute_current_total_n(&share);
				AllocateYr(&site, &veg, &climday, annualRstep, &share, CN_Mode);
				const double totalN_after_allocateyr = compute_current_total_n(&share);
				if (fabs(totalN_after_allocateyr - totalN_before_allocateyr) > 1.0e-8)
				{
					printf("FUN NDIAG AllocateYr year %d doy %d dN=%.10f before=%.10f after=%.10f\n",
						climday.year[rstep], climday.doy[rstep], totalN_after_allocateyr - totalN_before_allocateyr,
						totalN_before_allocateyr, totalN_after_allocateyr);
				}
				storeoutput(&veg, &share, &out, annualRstep, &ystep, 1);
				YearInit(&share);
			}
		}

		printf(" PnET-FUN-CORPSE executing year %d  doy  %d  \n", climday.year[rstep], climday.doy[rstep]);
	}

	WriteoutYr(&climday, &out);
	memfree_all(&site, &clim, &share, &out);
	Disturbance(&site, &veg, &share, 100);

	printf("   ===============================================\n");
}
