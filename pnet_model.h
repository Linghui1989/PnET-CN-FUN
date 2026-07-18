// pnet_model.h: interface for the pnet_model class.
//
//////////////////////////////////////////////////////////////////////


#include <stdio.h>
#include <stdlib.h>
//#include <cstddef.h>
//#include <stddef.h>
#include <algorithm>    // std::min max
//#include <io.h>
#include <float.h>
#include <string.h>

//#include <string.h>
#include <math.h>
#include <numeric>


#include "pnet_constants.h"
#include <cmath>



class pnet_model
{
	public:
	pnet_model();
	virtual ~pnet_model();

	enum { BG_FAST = 0, BG_SLOW = 1, BG_DEADMIC = 2, BG_CLASS_COUNT = 3 }; // FUN-CORPSE substrate classes: fast, slow, and dead microbial material


	// paths to files
	char sep[3];			// file separator. Wins:\\; Linux: /
	char exePath[300];		// program location
	char PathInput[100];		//input directory path
	char PathOutSite[30];	//output directory path for site mode
	char PathOutRegion[30];	//output directory path for regional mode
	char PathInter[30];		//intermediate directory path to store temporary files
	char PathLib[30];		// directory path to library folder. not used currently
	char PathRegion[30];	// path to regional input files.



	// model type

	int modelmode;  // 1: site modeling or 2:regional modeling
	int modeltype;  // model type, 0: pnet-day, 1: pnet-ii, 2: pnet-cn, 3: pnet-daily, 4: pnet-fun-corpse



	typedef struct
	{

		int	timestep;  // timestep of the climate data: 0: monthly; 1: daily; 2: hourly
		int	length; //length of climate record

		int* year; 						// year
		int* month;   					// month of the doy
		int* day;  						// date of the doy
		double* hour;						// hour of the doy  // only for hourly version
		int* ifday;						// if day time or night time for hours
		int* doy;						// day of year

		double* tmax;						// oC
		double* tmin;						// oC
		double*	tmean;						// oC

		double* par;						//  umol/m2/s 
		double* prec;						// cm
		double*	vpd;						// kPa

		double* O3;							// ppm-h, dose of hourly concentrations above a throshold  of 40 ppb ozone dose values summed over the growing season (may through oct)??
		//The D40 dose is accumulated over the entire growing season and is calculated as the sum of all
		//daytime hourly values > 40 ppb after subtracting 40 from from May through October.
		double* CO2;						// ppm
		double* NH4dep;						// g N/m2
		double* NO3dep;						// g N/m2
		int *ifleap;						// leap year, 1: yes; 0: no
	
	} clim_struct;

	//a structure to hold site data
	typedef struct
	{
		double	WHC;						// cm
		double	Lat;						// decimal degree without minute and second
		double	WaterStress;				// 0 or 1: 1 for no stress on photosynthesis
		double	SnowPack;					// cm equivalent water
		int		CO2gsEffect;				// 1: Yes; 0: NO. CO2 effect on gs to modify WUE and O3 effect.
		int		O3EffectOnWUE;				// 1: Yes; 0: NO, if there is an effect of O3 on WUE
		int		Nabs;				// 1: Yes; 0: NO, N deposition has absulate unit, i.e, gN m-2

		int		nYearStart;  		// simulation first year
		int		nYearEnd;   		 // end year of simulation
		double	SWater;  			// soil water
		double	SoilMoistFact; 		// soil moisture factor to soil respiration
		double	FastFlowFrac;		//Fraction of water inputs lost directly to drainage through macro pore or runoff
		double  WaterReleaseFact;	//Soil water release ability for evapotranspiration
		double	DWater;				// water stress index
		double  HOM;			//soil organic matter,	g/m2
		double	HON;			//soil organic N, gNm2
		double	NH4;			//soil ammonium N, gNm2
		double	Kho;			// soil decomposition constant, yr-1

		double	SoilRespA;		//the slope in the linear relation between soil resp CO2 and mean monthly temp, only in PnET-II
		double	SoilRespB;		//the intercept in the linear relation (g/m^2) between soil resp and mean monthly temp,only in PnET-II
		double	NImmobA;		// N immoblization rate parameter
		double	NImmobB;		// linear coefficients for fraction of mineralized N reimmobilized as a function of SOM C:N

		char ClimateFileName[300];		// climate file path, only used in GUI version

		// disturbance
		int		agstart;					// agriculture start year
		int		agstop;						// agriculture end year
		double	agrem;						// removal fraction in agriculture


		int	distyrs;					    // number of disturbances, added by J.X.
		double* distintensity;				// disturbance motality, 0-1
		int* distyear;						// disturbance year
		double* distremove;					// removal fraction for aboveground biomass
		double* distsoilloss;				// removal fraction for soil biomass				
		double* folregen;					// foliar regeneration rate next year, g/m-2
		int* distdoy;						// disturbance day of the year.
		double* distRootloss;				// disturbance of root masss change

		int FertYrStart;					// N fertilization start year
		int FertDOYStart;					// N fertilization start day of year
		int FertYrEnd;						// N fertilization end year
		int FertDOYEnd;						// N fertilization end day of year
		double FertNH4;						// N fertilization  NH4 rate gN/m2
		double FertNO3;						// N fertilization  NO3 rate gN/m2
		double FertUrea;					// N fertilization  urea rate gN/m2

		double Kden;
		double Albedo;						// dynamic snow albedo, in %

		int agnum;
		//	double
		int yrurstart;					// urbanization start year
		int yrurstart2;					// urbanization stage 2

	} site_struct;


	//a structure to hold veg data (parameters)
	typedef struct
	{
		double	AmaxA;				// the intercept in a linear relation max net Photosynthesis as a function of N.
		double	AmaxB;				//the slope in the function. Max rate units are umol CO2 m-2 leaf s-1
		double	AmaxFrac;			//Daily Amax as a fraction of the early morning instantaneous rate
		double	BaseFolRespFrac;	//respiration as a fraction of max photosynthesis
		double	CFracBiomass;		//Carbon as fraction of foliage mass
		double	DVPD1;				//Coefficients for converting VPD to DVPD (kPa-1	
		double	DVPD2;				//DVPD being fraction loss of photosynthesis
		double	FolMassMax;			//g/m2
		double	FolMassMin;			////g/m2
		double	GDDFolStart;		//growing degree days at which foliar production begins
		double	GRespFrac;			//Growth respiration, fraction of allocation
		double	GDDFolEnd;			//growing degree days at which foliar production ends
		double	HalfSat;			//Half saturation light level (umol m-2 s-1)
		double	PsnTMin;			//Minimum temperature for photosynthesis (C)
		double	PsnTOpt;			//Optimum temperature for photosynthesis (C)
		double	RespQ10;			//Q10 value for foliar respiration
		double	SLWdel;				//change in SLW with increasing foliar mass above (g/m-2 g-1)
		double	SLWmax;				//specific leaf weight (top of canopy) g/m-2
		double	SenescStart;		//first day of year when leaf drop could potentially start
		double	k;					//canopy light attenuation constant exp -kT
		double	FolNCon;			//2.2% for estimate from spectral data HF, HW, pine, spruce
		double	FolRelGrowMax;		//
		double	RootAllocA;			//Relationship between foliar and root allocation. Intercept
		double	RootAllocB;			//Relationship between foliar and root allocation. Slope
		double	WoodMRespA;			//Wood Maintenance respiration as a fraction of gross photowynthesis
		double	RootMRespFrac;		//Ratio of fine root maintenance respiration to biomass production
		double	PlantCReserveFrac;	//Fraction of PlantC held in reserve after allocation to BudC
		double	MinWoodFolRatio;	//Min ratio of carbon allocation to wood and foliage
		double	WUEConst;			//Constant in equation for water use efficiency as a function of VPD, g CO2/kg water
		double	PrecIntFrac;		//Fraction of precipitation intercepted and evaporated
		double	FastFlowFrac;		//Fraction of water inputs lost directly to drainage
		double	f;					//Soil water release parameter for evapotranspiration
		double	SoilRespA;			//the slope in the linear relation between soil resp CO2 and mean monthly temp
		double	SoilRespB;			//the intercept in the linear relation (g/m^2) between soil resp and mean monthly temp
		double	SoilMoistFact;		//effect of soil moisture
		double	WoodTurnover;		//fractional mortality of live wood per year, live wood to dead wood
		double	RootTurnoverA;		//constant term of 2-deg polynomial describing fine root turnover
		double	RootTurnoverB;		//as a function of annual net N mineralization. Linear coefficient
		double	RootTurnoverC;		//quadratic coefficient
		double	FolReten;			//Foliar retention time (yr)
		double	Kho;				//Decomposition constant for SOM pool yr-1 (SOM*(1-exp(-Kho*Tmult))
		double	NImmobA;			//linear coefficients for fraction of mineralized N
		double	NImmobB;			//reimmobilized as a function of SOM C:N
		double	MaxNStore;			//max N content in PlantN pool g N m-2
		double	GDDWoodStart;		//Growing degree days at which wood production begins
		double	GDDWoodEnd;			//Growing degree days at which wood production ends
		double	FolNConRange;		//max fractional increase in N concentration
		double	FolNRetrans;		//fraction of foliage N retransfer to plant N, remainder in litter 
		double	FLPctN;				//min % N concentration in foliar litter
		double	RLPctN;				//min N % cincentration in root litter
		double	WLPctN;				//min N % cincentration in wood litter
		double	WoodLitLossRate;	// fraction of dead wood to litter  yr-1
		double	WoodLitCLoss;		// fraction of litter decayed as CO2 

		int		TreeCode;			// not used currently.
		int		Age;				// not used currently.
		double	FolRespFrac;	    //Daily respiration as a fraction of monthly temperature results


	}  veg_struct;

	// FUN-CORPSE state for one compartment: organic layer, mineral rhizosphere, or mineral bulk soil
	typedef struct
	{
		double unprotC[BG_CLASS_COUNT];		// unprotected carbon pools for fast, slow, and dead microbial material, g C m-2
		double unprotN[BG_CLASS_COUNT];		// unprotected nitrogen pools for fast, slow, and dead microbial material, g N m-2
		double protC[BG_CLASS_COUNT];		// protected carbon pools for fast, slow, and dead microbial material, g C m-2
		double protN[BG_CLASS_COUNT];		// protected nitrogen pools for fast, slow, and dead microbial material, g N m-2
		double livingMicrobeC;			// living microbial biomass carbon in this compartment, g C m-2
		double livingMicrobeN;			// living microbial biomass nitrogen in this compartment, g N m-2
		double CO2;					// CO2 released from this compartment during the current step, g C m-2
		double Rtot;					// total respiration from this compartment during the current step, g C m-2

		double lastDecompC[BG_CLASS_COUNT];	// diagnostic carbon decomposed from each substrate class during the current step, g C m-2
		double lastDecompN[BG_CLASS_COUNT];	// diagnostic nitrogen decomposed from each substrate class during the current step, g N m-2
		double lastGrossMineralization;		// diagnostic gross N mineralization in the current step, g N m-2
		double lastGrossImmobilization;		// diagnostic gross N immobilization in the current step, g N m-2
		double lastProtectionN;			// diagnostic N moved into protected pools during the current step, g N m-2
		double lastDeprotectionN;		// diagnostic N released from protected pools during the current step, g N m-2
		double lastTurnoverN;			// diagnostic microbial turnover N produced in the current step, g N m-2

	} FUN_CORPSE_Compartment;

	// FUN-CORPSE parameter set used by all belowground compartments
	typedef struct
	{
		double vmaxref[BG_CLASS_COUNT];		// reference maximum decomposition rate for each substrate class
		double Ea[BG_CLASS_COUNT];		// activation energy for each substrate class, J mol-1
		double kC[BG_CLASS_COUNT];		// half-saturation constant for substrate carbon, g C m-2
		double eup[BG_CLASS_COUNT];		// microbial carbon uptake efficiency for each substrate class
		double nup[BG_CLASS_COUNT];		// microbial nitrogen uptake efficiency for each substrate class

		double max_immobilization_rate;		// maximum daily fraction of mineral N that microbes can immobilize
		double gas_diffusion_exp;		// exponent controlling moisture limitation on substrate access
		double minMicrobeC;			// minimum microbial biomass as a fraction of unprotected C
		double Tmic;				// microbial turnover time, years
		double et;				// fraction of turnover routed to litter rather than immediate CO2
		double CN_microb_default;		// fixed microbial C:N ratio when dynamic microbial C:N is disabled
		double tProtected;			// turnover time of protected SOM pools, years
		double protection_rate[BG_CLASS_COUNT];	// rate of transfer from unprotected to protected pools for each substrate class
		double frac_turnover_min;		// fraction of microbial turnover N released directly to mineral N
		double frac_turnover_slow;		// fraction of microbial turnover routed to slow substrate pools

		double leaf_fast_frac;			// fraction of foliar litter entering the fast substrate pool
		double leaf_slow_frac;			// fraction of foliar litter entering the slow substrate pool
		double root_fast_frac;			// fraction of root litter entering the fast substrate pool
		double root_slow_frac;			// fraction of root litter entering the slow substrate pool
		double deadwood_fast_frac;		// fraction of deadwood-derived litter entering the fast substrate pool
		double deadwood_slow_frac;		// fraction of deadwood-derived litter entering the slow substrate pool
		double initial_protected_fraction;	// initial fraction of legacy HOM/HON assigned to protected bulk-soil pools
		double initial_organic_fraction;	// initial fraction of legacy HOM/HON assigned to the organic horizon
		double initial_rhizo_fraction;		// initial fraction of legacy mineral HOM/HON assigned to the rhizosphere compartment
		double initial_active_n_scalar_organic;	// scalar applied to initial organic-horizon active-pool N (fast plus dead microbial unprotected N) to tune early N mineralization
		double initial_active_n_scalar_mineral;	// scalar applied to initial mineral-soil active-pool N (fast plus dead microbial unprotected N in rhizo and bulk) to tune early N mineralization
		double initial_bulk_slow_n_scalar;	// scalar applied to initial bulk-soil slow-pool N (protected plus unprotected) to tune early mineral-horizon N mineralization
		double initial_bulk_deadmic_n_scalar;	// scalar applied to initial bulk-soil dead-microbial-pool N (protected plus unprotected) to tune early mineral-horizon N mineralization
		double litter_transfer_to_soil;		// daily fraction of organic-layer material transferred to mineral soil
		double rhizo_fraction;			// fallback fixed rhizosphere fraction used when no daily seasonal forcing file is available
		int root_inputs_only_rhiz;		// flag for routing root-derived inputs only to rhizosphere

		double f_exudate_Nacq;			// fraction of additional N acquisition C allocated to root exudation
		double f_myc_Nacq;			// fraction of additional N acquisition C allocated to mycorrhiza
		double C_Nacq_max_frac_GPP;		// maximum additional belowground N acquisition C as a fraction of daily GPP
		double K_Ndef;				// half-saturation constant for plant N deficit response, g N m-2
		double exudate_CN;			// assumed exudate C:N ratio when exudate N is represented
		double Vmax_myc_N;			// maximum mycorrhizal organic N acquisition rate, g N m-2 d-1
		double K_myc_C;				// half-saturation constant for mycorrhizal C allocation, g C m-2 d-1
		double K_myc_N;				// half-saturation constant for organic N availability, g N m-2
		double eta_myc_to_plant;		// fraction of mycorrhizally acquired N transferred to the plant
		double myc_cost_n;			// coefficient on mineral N availability in the FUN mycorrhizal acquisition cost function
		double myc_cost_c;			// coefficient on root biomass in the FUN mycorrhizal acquisition cost function
		double nonmyc_cost_n;			// coefficient on mineral N availability in the non-mycorrhizal acquisition cost function
		double nonmyc_cost_c;			// coefficient on root biomass in the non-mycorrhizal acquisition cost function
		double storage_cost;			// baseline storage/retranslocation acquisition cost used in the simplified FUN optimizer
		double litter_fungi_frac;		// fraction of mycorrhizal fungal production routed to the litter layer fast pool
		double decomp_scalar_global;		// global scalar multiplier applied equally to decomposition rates in all horizons and substrate pools
		double decomp_scalar_organic;		// scalar multiplier on decomposition rates in the organic horizon
		double decomp_scalar_rhizo;		// scalar multiplier on decomposition rates in the rhizosphere horizon
		double decomp_scalar_bulk;		// scalar multiplier on decomposition rates in the bulk-soil horizon
		double gross_nmin_cap_frac_yr;		// FUN-only upper bound on gross N mineralization as a fraction of compartment organic N per year
		double nitrification_cap_frac_day;	// FUN-only upper bound on nitrification as a fraction of the NH4 pool per day
		double ndrain_cap_frac_day;		// FUN-only upper bound on nitrate leaching as a fraction of the NO3 pool per day

		double CN_microb_min;			// lower bound of dynamic microbial C:N ratio
		double CN_microb_max;			// upper bound of dynamic microbial C:N ratio
		double K_Navail_CN;			// half-saturation constant linking available N to microbial C:N, g N m-2

		int use_dynamic_microbe_cn;		// flag for dynamic versus fixed microbial C:N
		int enable_root_exudation;		// flag for root exudation pathway
		int enable_mycorrhizal_uptake;		// flag for mycorrhizal N uptake pathway

	} FUN_CORPSE_Params;

	// FUN-CORPSE belowground state aggregated across all soil compartments
	typedef struct
	{
		FUN_CORPSE_Compartment organicLayer;		// FUN-CORPSE state for the surface organic layer
		FUN_CORPSE_Compartment mineralRhizosphere;	// FUN-CORPSE state for mineral rhizosphere soil
		FUN_CORPSE_Compartment mineralBulk;		// FUN-CORPSE state for mineral bulk soil

		double rootExudateC;			// root exudate carbon added during the current step, g C m-2
		double mycorrhizalC;			// carbon allocated to mycorrhiza during the current step, g C m-2
		double mycorrhizalNToPlant;		// nitrogen transferred from mycorrhiza to plant during the current step, g N m-2
		double microbialCN;			// active microbial C:N ratio used during the current step

		double grossNMobilization;		// total gross N mobilization across all compartments during the current step, g N m-2
		double grossNImmobilization;		// total gross N immobilization across all compartments during the current step, g N m-2
		double netNMineralization;		// total net N mineralization across all compartments during the current step, g N m-2
		double microbialTurnoverN;		// total microbial turnover N across all compartments during the current step, g N m-2
		double deprotectionN;			// total N released from protected pools during the current step, g N m-2
		double protectionN;			// total N moved into protected pools during the current step, g N m-2

	} FUN_CORPSE_State;

	//a structure to hold share data
	typedef struct
	{
		double Tave;					// oC
		double Tday;					// oC
		double Tnight;					// oC
		double DayLength;				// second
		double NightLength;				// second
		double VPD;						// kpa
		double dayspan; 				// numbers of days of each timestep
		double GDDTot;					// total growth degree days
		double OldGDDFolEff;			// to calculate foliar growth
		double FolLitM;					// foliar litter mass, g/m2
		double PosCBalMass;				// possible C mass at balance point
		double PosCBalMassTot;			// total potential mass
		double PosCBalMassIx;			// total potential mass days
		double LAI;						// leaf area index
		double DVPD;					// vpd effect on photosynthesis
		double RH;						// Relative humidity
		double Tw;						// wet-bulb, stull
		double DayResp;					// foliar respiration at daytime
		double NightResp;				// foliar respiration at nighttime
		double CanopyNetPsn;			// canopy net photosynthesis, gpp-foliar respiration
		double CanopyGrossPsn;			// gpp,//g C/m2 ground/day
		double Dwatertot;				// total water stress at growing season
		double DwaterIx;				// total water stressed days
		double GrsPsnMo;				// monthly gpp
		double NetPsnMo;				// monthly net psn
		double FolGRespMo;				// foliar growth respiration
		double WoodMRespYr;				// yearly wood maintenance respiration
		double CanopyGrossPsnActMo;		//monthly gpp modified by water stress and other stress
		double FolProdCYr;				// foliar npp, g C m-2
		double FolProdCMo;				// foliar npp each time step, g C m-2
		double FolGRespYr;				// foliar yearly growth respiration, g C m-2
		double RootProdCYr;				// root npp, g C m-2
		double RootMRespYr;				// root yearly maintenance respiration, g C m-2
		double RootGRespYr;				// root yearly growth respiration, g C m-2
		double SoilRespMo;				// soil respiration, g C m-2
		double SoilRespYr;				// soil yearly respiration, g C m-2
		double OldGDDWoodEff;			// to calculate wood growth
		double WoodProdCYr;				// wood npp, g C m-2
		double WoodGRespYr;				// wood yearly growth respiration, g C m-2
		double TotPsn;					// total net psn
		double MeanSoilMoistEff;		// soil moisture effect on som decay
		double Drainage;				// water drainage, as runoff + leaching
		double TotDrain;				// total drainage
		double TotEvap;					// total evaporation
		double TotTrans;				// total transpiration
		double TotPrec;					// total precipitation
		double TotWater;				// total soil water
		double TotGrossPsn;				// total gpp
		double NPPFolYr;				// foliar npp, g dry matter m-2
		double NPPWoodYr;				// wood npp, g dry matter m-2
		double NPPRootYr;				// root npp, g dry matter m-2
		double ET;						// ET
		double NEP;						// net ecosystem production, -NEE
		double NetCBal;					// net C balance, NEE
		double SoilDecResp;				// soil decomposition
		double BudN;					// bud N or foliar total N
		double SoilDecRespYr;			// yearly soil decomposition
		double WoodDecRespYr;			// yearly deadwood decay loss as CO2
		double DelAmax;					// Amax adjustor for CO2 effect
		double DWUE;					// WUE adjustor for CO2 effect
		double CanopyDO3Pot;			// potential O3 effect on photosynthesis for the whole canopy
		double DroughtO3Frac;			// drought effect on O3 effect
		double TotO3Dose;				// total O3 dose
		double RootMassN;				//root N
		double TotalLitterMYr;			// total yearly litter mass
		double TotalLitterNYr;			// total yearly litter N
		double GrossNImmobYr;			// total yearly gross N immoblized
		double GrossNMinYr;				// total yearly gross N mineralization
		double PlantNUptakeYr;			// yearly plant uptake N
		double NetNitrYr;				// yearly net Nitrification rate
		double NetNMinYr;				// yearly net mineralization rate
		double FracDrain;				// proportion of drainage to total water
		double NDrainYr;				// yearly N drained out
		double NDrain;					// drained N, gN m-2
		double NDrainMgL;				// drained N concentration in water, mgN l-1
		double WoodDecResp;				// wood decay
		double TotalLitterM;			// total litter mass
		double TotalLitterN;			// total litter N
		double FolN;					// total foliar N
		double FolC;					// total foliar C
		double TotalN;					// total N
		double TotalM;					// total mass
		double NO3;						// NO3 content
		double NH4;						// NH4 content
		double FolNCon;				// foliar N concentration
		double FolNConOld; 				// to store FolN for output.
		double NdepTot;				// total N deposition
		double NdepTot_prev;		// total N deposition in previous year


		//Shared variables with initial conditions
		double FolMass;   //double FolMass=veg.FolMassMin;   In PnET-Day only,g/m2
		double BudC;    //double BudC=(veg.FolMassMax-double FolMass)*veg.CFracBiomass;  In PnET-Day only
		double Water;		// soil water
		double DeadWoodM;	// dead wood mass
		double WoodC;		//  wood C pool for wood growth
		double PlantC;		// plant C pool to store non structure C
		double RootC;		// Root C pool for root dynamic growth
		double LightEffMin;	// minimum light effect for next year foliar growth
		double NRatio;			// N stress index
		double PlantN;			// plant N pool
		double WoodMass;		// wood mass
		double RootMass;		// root mass
		double HOM;				// soil som
		double HON;				// soil son
		double RootNSinkEff;	// root N uptake capability
		double WUEO3Eff;		// O3 effect on WUE
		double WoodMassN;		// live wood total N
		double DeadWoodN;		// dead wood total N
		double NRatioNit;		// Nitrification constant determined by Nratio
		double NetNMinLastYr; 	// previous year net N mineralizatio rate
		double DWater;	  		// water stress for plant growth

		double LightEffCBal;	// light effect at foliar light compensation point.
		double LightEffCBalTot;	// total light effect at foliar light compensation point at growing season
		double LightEffCBalIx;	// number of days for LightEffCBal > 0.

		double O3Effect[51];	// O3 effect for each canopy layer

		double avgPCBM;			// average light effect
		double AvgDWater;		// average water stress

		double TaveYr;			// annual average air T, degree
		double PARYr;			// annual average PAR, umol m-2 s-1 at daytime

		double WoodProdCMo;
		double RootProdCMo;

		double Tsoil;
		double WFPS;

		double RnY2;
		double RnX2;
		double RnY1;
		double RnX1;
		double RnMax;

		double FluxN2OYrDe;
		double FluxNOYrDe;
		double FluxN2OYr;
		double FluxNOYr;
		double FluxN2Yr;

		double TotCanopyGrossPsnMax;


		int O3EffectOnPSN;
		int ifO3EffectOnPSN;
		double kO3Eff;		// empirically derived ozone response coefficient, 2.6 * power(10,-6)/ ppb-O3 for hardwood forests
		double CanopyDO3;	//  O3 effect on photosynthesis for the whole canopy
		double CanopyDO3Avg;	// Avg O3 effect on photosynthesis for the whole canopy for the whole growing season
		double CanopyDO3Tot;			// Total O3 effect on photosynthesis for the whole canopy for the whole growing season
		double CanopyDO3TotDay;			// Total days of O3 effect on photosynthesis for the whole canopy for the whole growing season
		double kO3EffScalar;		// empirically derived ozone response coefficient, 2.6 * power(10,-6)/ ppb-O3 for hardwood forests

		double EnvMaxFol;  //Linghui Meng 05072020 
		double SppMaxFol;	//Linghui Meng 05072020
		double FolGDD;
		int yrspin; // number of year for spin up
		int distid;


		double ISA0;					// The original ISA
		double ISA1;					// ISA in 2020
		double ISA2;					// ISA change from 2020-2050	
		double ISA;					//Linghui add ISA to the PnET model 20220811
		int In;						// the status of urbanization, 0 is off, 1 is urbanzationing, 2 is urbanzationed

		//parameter for new Phenology teets et al., 2023 //Linghui 07252024
		double CDDTot;
		double Tsum; //new temperature summary for phenology
		double Sc;
		double param;

		//add parameter for soil freezing
		double freeze;
		double freezeDoy;

		double ColdContent;
		double SnowLiq; //snow liquid 
		double SnowNH4;
		double SnowNO3;

		double SnowStatus; //check if there has snowpack

		//add parameter for root loss
		double OldRoot;
		double UptakeEff;

		FUN_CORPSE_State funBG;		// FUN-CORPSE belowground dynamic state shared across functions
		FUN_CORPSE_Params funParams;	// FUN-CORPSE parameter set shared across functions
		double FUNRhizoFractionByDOY[367];	// daily rhizosphere fraction forcing indexed by day of year, 1-366
		double FUNRhizoFractionCurrent;	// rhizosphere fraction applied in the current simulation step
		int FUNUseRhizoSeasonality;		// flag for enabling daily rhizosphere-fraction forcing from file

		double FolLitN;			// foliar litter nitrogen available to FUN-CORPSE mapping in the current step, g N m-2
		double RootLitM;			// root litter mass routed to FUN-CORPSE in the current step, g dry matter m-2
		double RootLitN;			// root litter nitrogen routed to FUN-CORPSE in the current step, g N m-2
		double WoodLitM;			// deadwood-derived litter mass routed to FUN-CORPSE in the current step, g dry matter m-2
		double WoodLitN;			// deadwood-derived litter nitrogen routed to FUN-CORPSE in the current step, g N m-2
		double PlantNUptakeLast;		// PnET active inorganic N uptake during the most recent step, g N m-2
		double FUNPotentialGrowthC;		// FUN-only pending growth carbon used to solve the current day's external N acquisition before final allocation, g C m-2
		double FUNPotentialNdemand;		// FUN-only pending plant N demand paired with FUNPotentialGrowthC for the current day, g N m-2
		double FUNPotentialRootProdC;		// FUN-only potential root production carbon before N limitation for the current step, g C m-2
		double FUNPotentialWoodProdC;		// FUN-only potential wood production carbon before N limitation for the current step, g C m-2
		double RootExudateCYr;		// yearly cumulative root exudate carbon input, g C m-2
		double MycorrhizalCYr;		// yearly cumulative mycorrhizal carbon allocation, g C m-2
		double MycorrhizalNToPlantYr;	// yearly cumulative mycorrhizal N transfer to plant, g N m-2
		double FUNFixNYr;			// yearly cumulative biological N fixation added through the FUN acquisition pathway, g N m-2
		double FUNNdemandYr;			// yearly cumulative FUN-side plant nitrogen demand diagnostic, g N m-2
		double FUNNdemandGapYr;		// yearly cumulative unmet FUN-side plant nitrogen demand before storage use, g N m-2
		double FUNPotentialRootProdCYr;	// yearly cumulative potential root production carbon before N limitation, g C m-2
		double FUNPotentialWoodProdCYr;	// yearly cumulative potential wood production carbon before N limitation, g C m-2
		double FUNNLimitRootProdCYr;	// yearly cumulative root production carbon lost to N limitation, g C m-2
		double FUNNLimitWoodProdCYr;	// yearly cumulative wood production carbon lost to N limitation, g C m-2
		double FUNPlantNOverflowYr;	// yearly cumulative amount by which FUN-side uptake pushes PlantN above MaxNStore before later allocation or return flow, g N m-2
		double FUNOverflowToNH4Yr;	// yearly cumulative PlantN overflow returned to NH4 in FUN runs, g N m-2
		double FUNOverflowToNO3Yr;	// yearly cumulative PlantN overflow returned to NO3 in FUN runs, g N m-2
		double BGGrossNMobilizationYr;	// yearly cumulative belowground gross N mobilization, g N m-2
		double BGGrossNImmobilizationYr;	// yearly cumulative belowground gross N immobilization, g N m-2
		double BGNetNMineralizationYr;	// yearly cumulative belowground net N mineralization, g N m-2
		double BGOrgGrossNMineralizationYr;	// yearly cumulative organic-horizon gross N mineralization, g N m-2
		double BGRhizGrossNMineralizationYr;	// yearly cumulative rhizosphere gross N mineralization, g N m-2
		double BGBulkGrossNMineralizationYr;	// yearly cumulative bulk-soil gross N mineralization, g N m-2
		double BGOrgGrossNImmobilizationYr;	// yearly cumulative organic-horizon gross N immobilization, g N m-2
		double BGRhizGrossNImmobilizationYr;	// yearly cumulative rhizosphere gross N immobilization, g N m-2
		double BGBulkGrossNImmobilizationYr;	// yearly cumulative bulk-soil gross N immobilization, g N m-2
		double BGOrgNetNMineralizationYr;	// yearly cumulative organic-horizon net N mineralization, g N m-2
		double BGRhizNetNMineralizationYr;	// yearly cumulative rhizosphere net N mineralization, g N m-2
		double BGBulkNetNMineralizationYr;	// yearly cumulative bulk-soil net N mineralization, g N m-2
		double BGMicrobialTurnoverNYr;	// yearly cumulative belowground microbial turnover N, g N m-2
		double BGProtectionNYr;		// yearly cumulative N transferred into protected pools, g N m-2
		double BGDeprotectionNYr;		// yearly cumulative N released from protected pools, g N m-2
		double BGMicrobialCNLast;		// most recent microbial C:N ratio used by FUN-CORPSE
		double AnnualStartTotalN;		// total ecosystem nitrogen at the start of the current year, g N m-2
		double AnnualStartSoilOrgN;		// soil organic nitrogen mirror (HON/FUN soil pools) at the start of the current year, g N m-2
		double AnnualStartMineralN;		// mineral nitrogen (NH4 + NO3) at the start of the current year, g N m-2
		double AnnualStartPlantStoreN;	// plant storage nitrogen (PlantN + BudN) at the start of the current year, g N m-2
		double AnnualStartVegStructN;	// structural vegetation nitrogen (foliar + wood + root) at the start of the current year, g N m-2
		double AnnualStartDeadWoodN;	// dead wood nitrogen at the start of the current year, g N m-2

	} share_struct;

	//a structure to hold output data
	typedef struct
	{
		int length;			// length of years

		double* grosspsn;			// monthly gross psn
		double* netpsn;				// monthly net psn
		double* netcbal;
		double* vpd;
		double* folmass;
		double* plantnMo;
		double* nppfol;
		double* nppwood;
		double* npproot;
		double* nep;
		double* gpp;
		double* waterstress;
		double* trans;
		double* soilwater;
		double* psn;				// yearly net psn
		double* drain;
		double* prec;
		double* evap;
		double* et;
		double* plantc;
		double* budc;
		double* woodc;
		double* rootc;
		double* plantnYr;
		double* budn;
		double* ndrain;
		double* netnmin;
		double* grossnmin;
		double* nplantuptake;
		double* grossnimob;
		double* bgnetnmin;		// yearly belowground net N mineralization diagnostic from FUN-CORPSE
		double* orggrossnmin;		// yearly organic-horizon gross N mineralization from FUN-CORPSE
		double* rhizgrossnmin;		// yearly rhizosphere gross N mineralization from FUN-CORPSE
		double* bulkgrossnmin;		// yearly bulk-soil gross N mineralization from FUN-CORPSE
		double* mineralgrossnmin;	// yearly mineral-soil gross N mineralization from FUN-CORPSE
		double* orgnetnmin_layer;	// yearly organic-horizon net N mineralization from FUN-CORPSE
		double* rhiznetnmin_layer;	// yearly rhizosphere net N mineralization from FUN-CORPSE
		double* bulknetnmin_layer;	// yearly bulk-soil net N mineralization from FUN-CORPSE
		double* mineralnetnmin_layer;	// yearly mineral-soil net N mineralization from FUN-CORPSE
		double* rootexudatec;		// yearly root exudate carbon diagnostic from FUN-CORPSE
		double* mycnplant;		// yearly mycorrhizal N transfer to plant diagnostic from FUN-CORPSE
		double* funfixn;		// yearly biological N fixation diagnostic from the FUN acquisition pathway, g N m-2
		double* littern;
		double* netnitrif;
		double* nratio;
		double* foln;
		double* litm;
		double* litn;
		double* rmresp;
		double* rgresp;
		double* decresp;
		double* totalsoilresp;		// yearly total soil respiration using a unified definition: root maintenance plus soil decomposition
		double* decwresp; //linghui 0729
		double* orgfastc;
		double* orgfastn;
		double* orgslowc;
		double* orgslown;
		double* orgdeadc;
		double* orgdeadn;
		double* rhizfastc;
		double* rhizfastn;
		double* rhizslowc;
		double* rhizslown;
		double* rhizdeadc;
		double* rhizdeadn;
		double* bulkfastc;
		double* bulkfastn;
		double* bulkslowc;
		double* bulkslown;
		double* bulkdeadc;
		double* bulkdeadn;
		double* orgmicn;		// yearly microbial nitrogen in the organic horizon FUN-CORPSE compartment
		double* rhizmicn;		// yearly microbial nitrogen in the rhizosphere FUN-CORPSE compartment
		double* bulkmicn;		// yearly microbial nitrogen in the bulk soil FUN-CORPSE compartment
		double* totalmicn;		// yearly total microbial nitrogen across all FUN-CORPSE compartments
		double* orgprotn;		// yearly protected nitrogen in the organic horizon across all substrate classes, g N m-2
		double* rhizprotn;		// yearly protected nitrogen in the rhizosphere across all substrate classes, g N m-2
		double* bulkprotn;		// yearly protected nitrogen in the bulk soil across all substrate classes, g N m-2
		double* totalprotn;		// yearly total protected nitrogen across all FUN-CORPSE compartments, g N m-2
		double* protectionnyr;		// yearly cumulative N transferred into protected pools, g N m-2
		double* deprotectionnyr;	// yearly cumulative N released from protected pools, g N m-2
		double* funndemand;		// yearly cumulative FUN-side plant nitrogen demand diagnostic, g N m-2
		double* funndemandgap;		// yearly cumulative unmet FUN-side plant nitrogen demand before storage use, g N m-2
		double* funpotrootc;		// yearly cumulative potential root production carbon before N limitation, g C m-2
		double* funpotwoodc;		// yearly cumulative potential wood production carbon before N limitation, g C m-2
		double* funnlimitrootc;		// yearly cumulative root production carbon lost to N limitation, g C m-2
		double* funnlimitwoodc;		// yearly cumulative wood production carbon lost to N limitation, g C m-2
		double* funplantnoverflow;	// yearly cumulative temporary PlantN overflow created immediately after FUN uptake, g N m-2
		double* funoverflowtonh4;	// yearly cumulative PlantN overflow returned to NH4, g N m-2
		double* funoverflowtono3;	// yearly cumulative PlantN overflow returned to NO3, g N m-2
		double* nstart_total;		// total ecosystem nitrogen at the start of each year, g N m-2
		double* nend_total;		// total ecosystem nitrogen at the end of each year, g N m-2
		double* dn_total;		// annual change in total ecosystem nitrogen, g N m-2 yr-1
		double* nstart_soilorg;		// soil organic nitrogen at the start of each year, g N m-2
		double* nend_soilorg;		// soil organic nitrogen at the end of each year, g N m-2
		double* dn_soilorg;		// annual change in soil organic nitrogen, g N m-2 yr-1
		double* nstart_mineral;		// mineral nitrogen (NH4 + NO3) at the start of each year, g N m-2
		double* nend_mineral;		// mineral nitrogen (NH4 + NO3) at the end of each year, g N m-2
		double* dn_mineral;		// annual change in mineral nitrogen, g N m-2 yr-1
		double* nstart_plantstore;	// plant storage nitrogen (PlantN + BudN) at the start of each year, g N m-2
		double* nend_plantstore;	// plant storage nitrogen (PlantN + BudN) at the end of each year, g N m-2
		double* dn_plantstore;		// annual change in plant storage nitrogen, g N m-2 yr-1
		double* nstart_vegstruct;	// structural vegetation nitrogen at the start of each year, g N m-2
		double* nend_vegstruct;		// structural vegetation nitrogen at the end of each year, g N m-2
		double* dn_vegstruct;		// annual change in structural vegetation nitrogen, g N m-2 yr-1
		double* nstart_deadwood;	// dead wood nitrogen at the start of each year, g N m-2
		double* nend_deadwood;		// dead wood nitrogen at the end of each year, g N m-2
		double* dn_deadwood;		// annual change in dead wood nitrogen, g N m-2 yr-1
		double* ngas;			// annual gaseous nitrogen loss (N2O + NO + N2), g N m-2 yr-1
		double* nbalance_resid;		// annual total-N residual: dTotalN - (Ndep - NDrain - Ngas), g N m-2 yr-1
		double* folm;
		double* deadwoodm;
		double* woodm;
		double* rootm;
		double* hom;
		double* hon;
		double* ndep;

		double* fn2o;
		double* fno;
		double* fn2;
		double* fn2ode;
		double* fnode;
		double* v1;			// spare variables for temporary use // =Tavg
		double* v2;			// par
		double* v3;			// par
//		double* RootNSinkEff;



	} out_struct;


	// structure for meteorological variables

	struct meteorology {

		double ustar;                   // friction velocity, m s-1
		double ustarnew;                // updated friction velocity with new H, m s-1
		double rhova_g;                 // absolute humidity, g m-3
		double rhova_kg;                // absolute humidity, kg m-3
		double sensible_heat_flux;      // sensible heat flux, W M-2
		double H_old;                   // old sensible heat flux, W m-2
		double air_density;             // air density, kg m-3
		double T_Kelvin;                // absolute air temperature, K
		double dispersion[sze3][sze];   // Lagrangian dispersion matrix, s m-1
		double zl;                      // z over L, height normalized by the Obukhov length
		double press_kpa;               // station pressure, kPa
		double press_bars;              // station pressure, bars
		double press_Pa;                // pressure, Pa
		double pstat273;                // gas constant computations
		double air_density_mole;        // air density, mole m-3
		double relative_humidity;       // relative humidity, ea/es(T)
		double vpd;                     // vapor pressure deficit

	} met;


	/****************************
	**                         **
	**  STRUCTURE DEFINITIONS  **
	**                         **
	****************************/
	typedef struct
	{
		//	file init;             /* initialization file */
		//	file inmet;            /* input meteorological data file */
		//	file outmet;           /* output meteorological data file */
		char init[128];
		char inmet[128];
		char outmet[128];
		int nhead;             /* number of header lines in input met file */
		int ndays;             /* number of days of data in input file */
		int indewpt;           /* input dewpoint temperature flag (0=NO, 1=YES) */
		int outhum;            /* output humidity flag            (0=VPD, 1=VP) */
		int inyear;            /* input year flag                 (0=NO, 1=YES) */
							   //	char systime[100];     /* system time at start of simulation */
							   //  char header[200];      /* header string, written to all output files */
							   //   char outprefix[200];   /* output filename prefix */
	} control_struct;

	typedef struct
	{
		double base_elev;      /* base elevation, meters */
		double base_isoh;      /* base annual precip isohyet, cm */
		double site_lat;       /* site latitude, dec. degrees (- for south) */
		double site_elev;      /* site elevation, meters */
		double site_slp;       /* site slope, degrees */
		double site_asp;       /* site aspect, degrees */
		double site_isoh;      /* site annual precip isohyet, cm */
		double site_ehoriz;    /* site east horizon, degrees */
		double site_whoriz;    /* site west horizon, degrees */
		double tmax_lr;        /* maximum temperature lapse rate, deg C/1000m */
		double tmin_lr;        /* minimum temperature lapse rate, deg C/1000m */
	} parameter_struct;

	typedef struct
	{
		int *year;             /* array of year values */
		int *mon;             /* array of month values */
		int *day;             /* array of day values */


		int *yday;             /* array of yearday values */
		double *tmax;          /* array of base maximum temperature values */
		double *tmin;          /* array of base minimum temperature values */
		double *prcp;          /* array of base daily precipitation values */
		double *tdew;          /* array of base dewpoint temperature values */
		double *s_tmax;        /* array of site tmax values */
		double *s_tmin;        /* array of site tmin values */
		double *s_tday;        /* array of site daylight temperature values */
		double *s_tmean;        /* array of site daily average temperature values */
		double *s_prcp;        /* array of site prcp values */
		double *s_prcpw;        /* array of site prcp values in liquid phase  */
		double *s_hum;         /* array of site humidity values (VPD or VP, Pa) */
		double *s_srad;        /* array of site shortwave radiation values */
		double *s_dayl;        /* array of site daylength values */
		double *s_swe;         /* array of site snow water equivalent values (cm) */
		double *s_smw;         /* array of site snow melted water values (cm) */

	} data_struct;


	void PathDefine(int wins);			// to define the directory path
	void pnet_run(int wins);			// control specific pnet run


	void pnet_site();					// site mode pnet
//	void pnet_region();					// regional pnet
//	void ReadRegionFile();				//
//	void GetRegionPath(char*);


	// memory set and free functions
	void memset_out(int nyrs, out_struct* out);
	void memset_climate(clim_struct* clim);
	void memfree_all(site_struct* site, clim_struct* clim, share_struct* share, out_struct* out);
	void pnet_memfree(site_struct* site, clim_struct* clim, share_struct* share, out_struct* out);
	void memfree_climate(clim_struct* clim);
	void memfree_out(out_struct* out);
	void memfree_site(site_struct* site);
	//	void memfree_share(share_struct* share);



	// specific for PnET-Day 
	void pnet_day();

	// specific for PnET-II 
	void pnet_ii();
	void SoilResp(veg_struct* veg, share_struct* share, int rstep);


	// specific for PnET-CN
	void pnet_cn();
	int getdays(int year, int doy); // calculate number of days for each month
	void AllocateMo(veg_struct* veg, share_struct* share, int rstep, int CN_Mode);
	void AllocateMoFUN(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share, int CN_Mode);	// FUN-only final allocation wrapper: run standard AllocateMo after FUN uptake has updated PlantN for the current day
	void AllocateYr(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share, int CN_Mode);
	void AtmEnviron(site_struct* site, clim_struct* clim, int rstep, share_struct* share);
	void CNTrans(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share);
	void Decomp(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share);
	void Leach(clim_struct* clim, int restep, share_struct* share);
	void Phenology(veg_struct* veg, clim_struct* clim, int rstep, share_struct* share, int GrowthPhase);
	void Photosyn(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share);
	void Waterbal(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share);

	void initvars(site_struct* site, veg_struct* veg, share_struct* share);
	void init_out(int nyrs, out_struct* out);
	void YearInit(share_struct* share);

	void ReadInput(site_struct* site, veg_struct* veg, share_struct* share, char* inputname);
	void ReadClim(clim_struct* clim, char* climname);
	void ReadFUNCORPSEParams(share_struct* share, char* inputname);
	void ReadFUNRhizoSeasonality(share_struct* share, char* inputname);	// read daily rhizosphere fraction forcing for FUN-CORPSE from the input folder

	void storeoutput(veg_struct* veg, share_struct* share, out_struct* out, int rstep, int* ystep, int NewYear);
	void WriteoutMo(site_struct* site, veg_struct* veg, clim_struct* clim, share_struct* share, int rstep, FILE*& fileoutM);
	void WriteoutYr(clim_struct* clim, out_struct* out);
	void WriteoutYr_netcdf(clim_struct* clim, out_struct* out, char* filename);
//	void WriteoutRgn(int ngrid, out_struct* out, char* filename);
//	void StoreoutRgn(int grid, int nyrs, out_struct* out, out_struct* outrgn);




	// specific for PnET-Daily
	void pnet_daily();
	void pnet_fun_corpse();
	void ReadClimDay(clim_struct* clim, char* climname);
	void WriteoutDay(site_struct* site, veg_struct* veg, clim_struct* clim, share_struct* share, int rstep, FILE*& fileout);
	int getmonth(int year, int doy);  // get the month of doy
	int is_leapyear(int year); // 1: leap year; 0: non leap year

	void Photosyn_hour(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share);
	void ReadClimHour(clim_struct* clim, char* climname);

	// pre and post processes
	void clim_day_mon(site_struct* site, share_struct* share, clim_struct* climday, clim_struct* climmon); // daily climate to monthly
	void clim_day_fill(site_struct* site, share_struct* share, clim_struct* climday, clim_struct* climmon); // daily climate to monthly
	void write_clim(clim_struct* clim, int i);




	int Disturbance(site_struct* site, veg_struct* veg, share_struct* share, int scenario);
	void epscor_region_merrimack_future();
	void pnet_region_test(int argc, char** argv);
	void pnet_region(int argc, char** argv);
//	double CO2_Daily(int year, int doy, int ctl);


	/********************************
	**                             **
	**    FUNCTION PROTOTYPES      **
	**                             **
	********************************/


	void MTCLIM(clim_struct* clim, control_struct *ctrl, parameter_struct *p, data_struct *data);


	int read_init(control_struct *ctrl, parameter_struct *p);
	int data_alloc(control_struct *ctrl, data_struct *data);
	int read_inmet(control_struct *ctrl, data_struct *data);

	int calc_tair(control_struct *ctrl, parameter_struct *p,
		data_struct *data);
	int calc_prcp(control_struct *ctrl, parameter_struct *p,
		data_struct *data);
	int snowpack(control_struct *ctrl, parameter_struct *p,
		data_struct *data);
	int calc_srad_humidity(control_struct *ctrl, parameter_struct *p,
		data_struct *data);
	int calc_srad_humidity_iterative(control_struct *ctrl, parameter_struct *p, data_struct *data);

	//int calc_srad_humidity_iterative( control_struct *ctrl,	parameter_struct *p, data_struct *data,
	//	double *dtr,double *sm_dtr, double *parray,double *window,double *t_fmax,double *tdew,double *save_pet );

	int write_out(control_struct *ctrl, data_struct *data, clim_struct* clim);
	int data_free(control_struct *ctrl, data_struct *data);
	double calc_pet(double rad, double ta, double pa, double dayl);
	double atm_pres(double elev);
	int pulled_boxcar(double *input, double *output, int n, int w, int w_flag);




	



	double Lat;						// latitude, degree
	double Lon;						// longitude, degree
	double CanopyN;					// canopy N ,%
	double Elev;						// elevation,m
	double Slope;					// slope, degree
	double Aspect;					// aspect, degree (N ,0, E 90, S 180)
	double FracDeci;						//Fraction of Deciduous stands,0-1
	double FracEver;						//Fraction of Evergreen stands,0-1, 1-FracDeci
	double FolNDeci;					//Foliar N of Deciduous stands,%
	double FolNEver;					//Foliar N of Evergreen stands,%

	int grid;




	int epscor_region();
	int epscor_region_HB_obs();
	int epscor_region_merrimack();  // contemporary run and validation
	int  epscor_region_post();

	void Photosyn_daily(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share);



	void preprocess_climate(site_struct* site, share_struct* share, clim_struct* climday, clim_struct* climmon);

	// FUN-CORPSE initialization and update routines
	void InitializeFUN_CORPSEBelowground(site_struct* site, share_struct* share);
	void DecompFUN_CORPSE(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share);
	void UpdateFUN_CORPSEBelowground(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share, double TMult, double WMult);
	void UpdateCORPSECompartment(const char* label, FUN_CORPSE_Compartment* comp, double soilT_K, double theta, double* NH4, double* NO3, double dt_years, FUN_CORPSE_Params* params, double microbialCN, double decompScalar);

	// FUN-CORPSE helper functions for microbial stoichiometry and plant N demand
	double ComputeDynamicMicrobialCN(share_struct* share);
	double ComputeAdditionalNacqCAllocation(veg_struct* veg, share_struct* share);
	double ComputePlantNDeficit(veg_struct* veg, share_struct* share);
	void EstimateFUNAllocationDemand(veg_struct* veg, clim_struct* clim, int rstep, share_struct* share, int CN_Mode);	// FUN-only pre-allocation demand estimator: cache today's potential growth C and N demand for later FUN uptake
	double ComputeFUNPassiveNUptake(site_struct* site, veg_struct* veg, share_struct* share, double NdemandTotal);	// FUN-only passive N uptake using the legacy PnET root-uptake formulation, capped at 80 percent of currently available mineral N
	void ApplyFUNCostBasedNUptake(site_struct* site, veg_struct* veg, share_struct* share, int doy, double Cavailable, double NdemandTotal);	// FUN-only N acquisition optimizer: solve cost-based external N uptake for the current day's growth C and N demand before final allocation
	double GetFUNRhizoFractionForDOY(share_struct* share, int doy);	// return the rhizosphere fraction forcing used for a given day of year
	void ApplyFUNRhizoSeasonality(share_struct* share, int doy);	// redistribute mineral FUN-CORPSE pools to the target rhizosphere fraction for the current day

	// FUN-CORPSE transfer, synchronization, and diagnostic helpers
	void TransferOrganicLayerToMineralSoil(share_struct* share);
	void SyncLegacySoilPoolsFromFUN(share_struct* share);
	void ResetFUN_CORPSEFluxDiagnostics(share_struct* share);

	// PnET inorganic N uptake retained as the active plant uptake pathway
	void PnETActiveRootUptake(site_struct* site, veg_struct* veg, share_struct* share, double TMult);

	//int netcdf_sfc_pres_temp_rd();
};
