#include "pnet_model.h"
#include <netcdf>
#include <iostream>
#include <mpi.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

void pnet_model::GetRegionPath(char* PathProj)
{
	char inputname[400];
	FILE* f_input_load;
	char  note[400];

	sprintf(inputname, "%s%s%s", exePath, PathInput, "input.txt");

	f_input_load = fopen(inputname, "r");

	if (f_input_load == NULL)
	{
		printf("Unable to open input file %s!\n", inputname);
		exit(1);
	}


	for (int i = 1; i <= 5; i++)fgets(note, 200, f_input_load);

	fscanf(f_input_load, "%s\n", inputname); fgets(note, 200, f_input_load);

	fclose(f_input_load);

	sprintf(PathProj, "%s%s%s%s", exePath, PathRegion, inputname, sep);


}

void pnet_model::ReadRegionFile()
{

	char inputname[400];
	FILE* f_input_load;
	char  note[400];

	sprintf(inputname, "%s%s%s", exePath, PathInput, "Region.rgn");

	f_input_load = fopen(inputname, "r");

	if (f_input_load == NULL)
	{
		printf("Unable to open input file %s!\n", inputname);
		exit(1);
	}


	fgets(note, 200, f_input_load);

	for (int i = 1; i <= 5; i++)fgets(note, 200, f_input_load);

	fscanf(f_input_load, "%s\n", inputname); fgets(note, 200, f_input_load);

	fclose(f_input_load);


}


void pnet_model::pnet_region(int argc, char** argv)
{
	/*****************************************************************************************************************
	NOTE:

	1, use keyword "DEFINE" to find variables that users have to set for different simulations
	2,


	***************************************************************************************************************/

	char inputname[450];
	char note[200];
	char climtmax[450], climtmin[450], climpar[450], climppt[450];  //input climate file location
	char pathCO2[450], pathO3[450], pathNH4[450], pathNO3[450];	 //path of air chemistries
	int CN_Mode = 1;											 //if running PnET-CN
	int ystep = 1;												 //unlike matlab code, ystep is initialized here instead of in initvars.c
	char ndep[200];

	FILE* ndephb;												 // hb n deposition
	int rstep;
	int ncycle, ryrs;											 //ncycle, number of climate record cycles for spin up; ryrs, number of record years


	int dayobs; // days of observations
	int ndays;
	int doyEnd;  // end of year for doy, leap year:366, otherwise 365
	int yr;
	int nv = 0;  // number of variables
	int rec, nsteps;




	//============DEFINE parameters===================
	int ctl_clm, ctl_lu, ctl_out;
	int ctl_co2;

	ctl_clm = 2;	// comtemporary-Merra
	ctl_lu = 0;		// Backyardm,bkam
	ctl_out = 1;
	ctl_co2 = ctl_clm;


	char lu_dir[400];
	sprintf(lu_dir, "%s%s%s%s", exePath, sep, "organized_data", sep); //directory of regional data

	char pnetresult[400];
	sprintf(pnetresult, "%s%s%s%s", exePath, sep, "result", sep);

	char fileout_yr[800];

	char climdayfile[450];

	sprintf(climdayfile, "%s%s%s", exePath, PathInput, "climateday_preindustrial.txt");




	int yrstart, yrend;		// climate starting year and ending year
	int yrendhis, yrspin;    // end of historical data

	if (ctl_clm == 1 || ctl_clm == 0 || ctl_clm == 3)  // GFDL -OR- Huber
	{
		yrstart = 1980;//1979;   		// DEFINE by users, starting year of climate
		yrendhis = 2015;// 2014
		yrend = 2015;//2099;		//2014	// DEFINE by users, ending year of climate
		ncycle = 20;//10;			// DEFINE, number of cycling for spin up simulation
	}

	if (ctl_clm == 2)  // MERRA
	{
		yrstart = 1000;//1979; 		// DEFINE by users, starting year of climate (simulation)
		yrendhis = 1900;// 2014     // End of spin up year
		yrend = 2050;//2014	        // DEFINE by users, ending year of climate
		ncycle = 56;//10;			// DEFINE, number of cycling for spin up simulation
	}




	//================set spatial data name and path==========================

	if (ctl_clm == 2)
	{
		ctl_co2 = 1;
		//	sprintf(climtmax, "%s%s%s%s", Path, "agg_macav2metdata_tasmax_", scenario, "_2006_2099_CONUS_daily.nc");
		sprintf(climtmax, "%s%s", lu_dir, "tmax_daily.nc");
		sprintf(climtmin, "%s%s", lu_dir, "tmin_daily.nc");
		sprintf(climpar, "%s%s", lu_dir, "par_daily.nc");
		sprintf(climppt, "%s%s", lu_dir, "ppt_daily.nc");
		sprintf(pathCO2, "%s%s", lu_dir, "CO2_daily.nc");
		sprintf(pathO3, "%s%s", lu_dir, "O3_daily.nc");
		sprintf(pathNO3, "%s%s", lu_dir, "NO3_daily.nc");
		sprintf(pathNH4, "%s%s", lu_dir, "NH4_daily.nc");
	}

	//================     read landcover data    ==========================
	//define the scope of netCDF data
	static const int NLAT = 122;
	static const int NLON = 127;
	static const int NREC = 54750;
	static const int NC_ERR = 2;

	//create array to hold landcover and temporate variable

	double nlcd[NLAT][NLON], lats[NLAT];
	sprintf(inputname, "%s%s%s", exePath, PathInput, "landcover_final.nc");
	NcVar var, latVar; //ncvariable to hold regional data
	NcFile datafile(inputname, NcFile::read);
	var = datafile.getVar("var");
	latVar = datafile.getVar("lat");
	latVar.getVar(lats);

	//create vector to hold the file
	vector<size_t> startp, countp;
	startp.push_back(0);
	startp.push_back(0);
	startp.push_back(0);
	countp.push_back(1);
	countp.push_back(NLAT);
	countp.push_back(NLON);

	var.getVar(startp, countp, nlcd);
	datafile.close();


	NcVar tmaxVar, tminVar, precVar, parVar, nh4Var, no3Var, co2Var, o3Var;
	//=====================read regional data===========================
	// read regional tmax data
	NcFile df1(climtmax, NcFile::read);
	tmaxVar = df1.getVar("tmax");


	// read regional tmin data
	NcFile df2(climtmin, NcFile::read);
	tminVar = df2.getVar("tmin");

	// read regional par data
	NcFile df3(climpar, NcFile::read);
	parVar = df3.getVar("par");

	// read regional ppt data
	NcFile df4(climppt, NcFile::read);
	precVar = df4.getVar("ppt");

	// read regional CO2 data
	NcFile df5(pathCO2, NcFile::read);
	co2Var = df5.getVar("CO2");

	// read regional O3 data
	NcFile df6(pathO3, NcFile::read);
	o3Var = df6.getVar("O3");

	// read regional NH4 data
	NcFile df7(pathNH4, NcFile::read);
	nh4Var = df7.getVar("NH4");

	// read regional NO3 data
	NcFile df8(pathNO3, NcFile::read);
	no3Var = df8.getVar("NO3");
	//	printf("read netcdf data to var variables \n");


		// ============initialize for regional running========================
		// Initialize MPI
	MPI_Init(&argc, &argv);

	int rank, size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Comm node_comm;

	// Determine the number of nodes and CPUs available
	int num_nodes = -1, num_cpus = -1;
	char* str_num_nodes = getenv("SLURM_NNODES");
	char* str_num_cpus = getenv("SLURM_CPUS_ON_NODE");
	if (str_num_nodes) {
		num_nodes = atoi(str_num_nodes);
	}
	if (str_num_cpus) {
		num_cpus = atoi(str_num_cpus);
	}

	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &node_comm);
	MPI_Comm_size(node_comm, &num_nodes);

	// Assign works to different CPUs on each node
	long unsigned int lat_start = rank * NLAT / size;
	long unsigned int lat_end = (rank + 1) * NLAT / size - 1;


	if (lat_end >= NLAT) {
		lat_end = NLAT - 1;
	}

	//==================================================================
	//=====================read regional data===========================
	int nlat_per_cpu = lat_end - lat_start + 1;


	clim_struct climhis, climday;	 //structure to hold the input climate data

	if (ctl_clm == 2)  // MERRA
	{
		yrstart = 1000;//1979; 		// DEFINE by users, starting year of climate (simulation)
		yrendhis = 1900;// 2014     // End of spin up year
		yrend = 2050;//2014	        // DEFINE by users, ending year of climate
		ncycle = 56;//10;			// DEFINE, number of cycling for spin up simulation
	}



	//calculate historical days
	yrspin = yrendhis - yrstart + 1;   // years of spinup
	climhis.timestep = 1;								// PnET historical time step
	climhis.length = yrspin * 365;							// for daily running,12784
	memset_climate(&climhis);							//allocate memory for the climate structure


	ryrs = yrend - yrendhis;			// years of records
	climday.timestep = 1;
	climday.length = ryrs * 365;
	memset_climate(&climday);


	//	printf("yrstart %d, yrendhis %d, yrend %d \n", yrstart, yrendhis, yrend);
	ReadClimDay(&climhis, climdayfile); //read historical data


	// allocate date to climday
	sprintf(inputname, "%s%s%s", exePath, PathInput, "date.txt");
	ndephb = fopen(inputname, "r");
	fgets(note, 200, ndephb);

	dayobs = climday.length;
	for (rec = 1; rec <= dayobs; rec++)
	{

		fscanf(ndephb, "%d %d \n", &climday.year[rec], &climday.doy[rec]);

		if (is_leapyear(climday.year[rec]))
			climday.ifleap[rec] = 1;
		else climday.ifleap[rec] = 0;

	}
	fclose(ndephb);


	out_struct out;
	memset_out(ryrs, &out);				//allocate memory for the output




	// ============initialize for structures========================
	site_struct site;
	veg_struct veg;
	share_struct share0;
	sprintf(inputname, "%s%s%s", exePath, PathInput, "input.txt");
	ReadInput(&site,&veg,&share0,inputname);		// read in input data from input files
	initvars(&site,&veg,&share0);


	// ============disturbance=============================================
	Disturbance(&site, &veg, &share0, 0);  // to allocate memory for disturbances
	Disturbance(&site, &veg, &share0, 2);  // default disturbance as HB

//	printf("set up model \n");

	for (int lat = lat_start; lat <= lat_end; lat++)
	{

		site.Lat = lats[lat];  // updating latitude

//		printf("executing lat is %d \n", lat);
		for (int lon = 0; lon < NLON; lon++) // HB
		{

			double nl = nlcd[lat][lon];
			//			site.ISA0 = isa0[lat][lon];
			//			site.ISA1 = (isa1[lat][lon] - isa0[lat][lon]) / 120;
			//			site.ISA2 = isa2[lat][lon] / 30;

			if (!isnan(nl))
			{
				share_struct share;
				share = share0;
				initvars(&site, &veg, &share);
				init_out(ryrs, &out);
				printf("executing lat is %d lon is %d \n", lat, lon);

				site.Lat = lats[lat];  // updating latitude
//				printf("executing lat is %d lon is %d \n", lat, lon);
//				printf("executing lat is %d \n", site.Lat);

				nsteps = climhis.length;
				//				printf("executing historical length is %d \n", nsteps);



				for (rstep = 1; rstep <= nsteps; rstep++)
				{
					AtmEnviron(&site, &climhis, rstep, &share);
					Phenology(&veg, &climhis, rstep, &share, 1);
					Photosyn(&site, &veg, &climhis, rstep, &share);
					Waterbal(&site, &veg, &climhis, rstep, &share);
					AllocateMo(&veg, &share, rstep, CN_Mode);
					Phenology(&veg, &climhis, rstep, &share, 2);
					CNTrans(&site, &veg, &climhis, rstep, &share);
					Decomp(&site, &veg, &climhis, rstep, &share);
					Leach(&share);

					//End-of-year activity
					if (climhis.doy[rstep] == 365)
					{
						AllocateYr(&site, &veg, &climhis, rstep + 1, &share, CN_Mode); // Note the rstep+1, not rstep
						YearInit(&share);
					}
				}// end of the time loop
//				write_clim(&climhis, nl);
				//create vector to hold the file
//				printf("finish historical run \n");

				//create vector to hold the file
				startp = { 0, lat, lon };
				countp = { 1, 1, 1 };

				// Read data from netCDF variables
				ndays = climday.length;
				for (rec = 0; rec < ndays; rec++)
				{
					// Read the data one record at a time.
					startp[0] = rec;
					tmaxVar.getVar(startp, countp, &climday.tmax[rec + 1]);
					tminVar.getVar(startp, countp, &climday.tmin[rec + 1]);
					parVar.getVar(startp, countp, &climday.par[rec + 1]);
					precVar.getVar(startp, countp, &climday.prec[rec + 1]);
					co2Var.getVar(startp, countp, &climday.CO2[rec + 1]);
					o3Var.getVar(startp, countp, &climday.O3[rec + 1]);
					nh4Var.getVar(startp, countp, &climday.NH4dep[rec + 1]);
					no3Var.getVar(startp, countp, &climday.NO3dep[rec + 1]);
				}


				//				write_clim(&climday, nl);
								// ===================       daily running      ==================
				//				printf("executing daily run length is %d \n", ndays);
				ystep = 1;
				for (rstep = 1; rstep <= ndays; rstep++)
				{
					AtmEnviron(&site, &climday, rstep, &share);
					Phenology(&veg, &climday, rstep, &share, 1);
					Photosyn(&site, &veg, &climday, rstep, &share);
					Waterbal(&site, &veg, &climday, rstep, &share);
					AllocateMo(&veg, &share, rstep, CN_Mode);
					Phenology(&veg, &climday, rstep, &share, 2);
					CNTrans(&site, &veg, &climday, rstep, &share);
					Decomp(&site, &veg, &climday, rstep, &share);
					Leach(&share);

					if (climday.doy[rstep] >= 365)
					{
						//						printf("the in share is %7.1f \n", share.WoodMass);
						AllocateYr(&site, &veg, &climday, rstep + 1, &share, CN_Mode); // Note the rstep+1, not rstep
						storeoutput(&veg, &share, &out, rstep + 1, &ystep, 1);//Linghui Meng 06302021
						YearInit(&share);
						//						printf("the deadwood in outfile is %7.1f \n", out.woodm[ystep]);
					}
				}

				//define the name of annual result
				sprintf(fileout_yr, "%s%1f%s", pnetresult, nl, "_Output_annual.csv");
				WriteoutYr_netcdf(&climday, &out, fileout_yr);
				//				printf("write out\n");

			}


		} // end of loop


	}


	//	printf("release memory\n");
	//=======================================================================//
	MPI_Barrier(MPI_COMM_WORLD);

	Disturbance(&site, &veg, &share0, 100);  // to free memory for disturbances
	memfree_out(&out);
	memfree_climate(&climday);



	// Terminate the sub-communicator
	MPI_Comm_free(&node_comm);
	// Finalize MPI
	MPI_Finalize();

	memfree_climate(&climhis);
	df1.close();
	df2.close();
	df3.close();
	df4.close();
	df5.close();
	df6.close();
	df7.close();
	df8.close();

	// Release my_var on rank 0
	if (rank == 0) {

		printf("   ===============================================\n");
		printf("             Model run ends\n");



	}


}

void pnet_model::StoreoutRgn(int grid, int nyrs, out_struct* out, out_struct* outrgn)
{

	//	int nyrs = clim->year[clim->length]-clim->year[1]+1;
	int avgyrs; //number of years to average the output
	avgyrs = 30;


	//	for (int i=1;i<nyrs+1;i++)
	for (int i = nyrs; i > nyrs - avgyrs; i--)
	{

		outrgn->folm[grid] += out->folm[i] / avgyrs;
		outrgn->deadwoodm[grid] += out->deadwoodm[i] / avgyrs;
		outrgn->woodm[grid] += out->woodm[i] / avgyrs;
		outrgn->rootm[grid] += out->rootm[i] / avgyrs;

		outrgn->hom[grid] += out->hom[i] / avgyrs;
		outrgn->hon[grid] += out->hon[i] / avgyrs;
		/**/
		outrgn->nppfol[grid] += out->nppfol[i] / avgyrs;
		outrgn->nppwood[grid] += out->nppwood[i] / avgyrs;
		outrgn->npproot[grid] += out->npproot[i] / avgyrs;
		outrgn->psn[grid] += out->psn[i] / avgyrs;		// total net psn
		outrgn->nep[grid] += out->nep[i] / avgyrs;
		outrgn->gpp[grid] += out->gpp[i] / avgyrs;


		outrgn->prec[grid] += out->prec[i] / avgyrs;
		outrgn->evap[grid] += out->evap[i] / avgyrs;
		outrgn->trans[grid] += out->trans[i] / avgyrs;
		outrgn->et[grid] += out->et[i] / avgyrs;
		outrgn->drain[grid] += out->drain[i] / avgyrs;
		outrgn->waterstress[grid] += out->waterstress[i] / avgyrs;

		outrgn->budc[grid] += out->budc[i] / avgyrs;
		outrgn->plantc[grid] += out->plantc[i] / avgyrs;
		outrgn->woodc[grid] += out->woodc[i] / avgyrs;
		outrgn->rootc[grid] += out->rootc[i] / avgyrs;


		outrgn->ndep[grid] += out->ndep[i] / avgyrs;
		outrgn->plantnYr[grid] += out->plantnYr[i] / avgyrs;
		outrgn->budn[grid] += out->budn[i] / avgyrs;
		outrgn->ndrain[grid] += out->ndrain[i] / avgyrs;
		outrgn->netnmin[grid] += out->netnmin[i] / avgyrs;
		outrgn->netnitrif[grid] += out->netnitrif[i] / avgyrs;

		outrgn->nratio[grid] += out->nratio[i] / avgyrs;
		outrgn->foln[grid] += out->foln[i] / avgyrs;
		outrgn->litm[grid] += out->litm[i] / avgyrs;
		outrgn->litn[grid] += out->litn[i] / avgyrs;

		outrgn->rmresp[grid] += out->rmresp[i] / avgyrs;
		outrgn->rgresp[grid] += out->rgresp[i] / avgyrs;
		outrgn->decresp[grid] += out->decresp[i] / avgyrs;

		outrgn->v1[grid] += out->v1[i] / avgyrs;
		outrgn->v2[grid] += out->v2[i] / avgyrs;
		/**/
	}

	init_out(nyrs, out);

}


void pnet_model::WriteoutRgn(int ngrid, out_struct* out, char* filename)
{

	FILE* fileoutY;

	fileoutY = fopen(filename, "w");
	if (fileoutY == NULL)
	{
		printf("Unable to open the Output_annual.txt !\n");
		exit(1);
	}


	fprintf(fileoutY, " grid, folm, deadwoodm, livewoodm, rootm, som, son,"
		"nppfol, nppwood, npproot, psn, nep, gpp,"
		"prec, evap, trans, et, drain, wtrstress,"
		"plantc, budc, woodc, rootc, ndep, plantnYr, budn, ndrain, netnmin, netnitrif, nratio, foln,"
		"litm, litn, rmresp, rgresp, decresp,"
		"temp, par"
		"\n");
	fprintf(fileoutY, " ,g/m2,g/m2,g/m2,g/m2,g/m2,gN/m2,"
		"g/m2,g/m2,g/m2,gC/m2,gC/m2,gC/m2,"
		"cm,cm,cm,cm,cm, ,"
		"gC/m2,gC/m2,gC/m2,gC/m2,gN/m2,gN/m2,gN/m2,gN/m2,gN/m2,gN/m2, , %%,"
		"gC/m2, gN/m2,gC/m2,gC/m2,gC/m2,"
		"oC, umol/m2/s  "
		"\n");




	for (int i = 1; i < ngrid + 1; i++)
	{

		fprintf(fileoutY, "%5d, ", i);  // sequence of grid	double*	folm;

		fprintf(fileoutY, "%8.2f, %8.2f, %8.2f", out->nppfol[i], out->nppwood[i], out->npproot[i]);// C ballance

		fprintf(fileoutY, "%8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f,", out->folm[i], out->deadwoodm[i],
			out->woodm[i], out->rootm[i], out->hom[i], out->hon[i]);// C ballance

		fprintf(fileoutY, "%8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f,", out->nppfol[i], out->nppwood[i],
			out->npproot[i], out->psn[i], out->nep[i], out->gpp[i]);// C ballance

		fprintf(fileoutY, "%8.2f, %8.2f, %8.2f, %8.2f, %8.2f, %8.2f,", out->prec[i], out->evap[i],
			out->trans[i], out->et[i], out->drain[i], out->waterstress[i]);// water cycle

		fprintf(fileoutY, "%8.2f, %8.2f, %8.2f, %8.2f, ", out->plantc[i], out->budc[i], out->woodc[i], out->rootc[i]);// carbon cycle

		fprintf(fileoutY, "%8.2f, %8.2f, %8.2f, %8.2f, %8.2f, ", out->ndep[i], out->plantnYr[i], out->budn[i], out->ndrain[i], out->netnmin[i]);// N cycle

		fprintf(fileoutY, "%7.2f, %7.2f, %7.2f, ", out->netnitrif[i], out->nratio[i], out->foln[i]);// N cycle

		fprintf(fileoutY, "%8.2f, %8.2f, %8.2f, %8.2f, %8.2f", out->litm[i], out->litn[i], out->rmresp[i], out->rgresp[i], out->decresp[i]);// N cycle

		fprintf(fileoutY, "\n");// N cycle


	}


	fclose(fileoutY);


}






double pnet_model::CO2_Daily(int year, int doy, int ctl)
{

	// year : xxxx
	// doy: 1-366
	// ctl: CO2 concentration scenarios


	double co2;
	int i,j;
	const static int years[]={1970,1980,1990,2000,2010,2020,2030,2040,2050,2060,2070,2080,2090,2100};
	const static double co2scn[2][14]={
			325,337,353,369,389,417,455,504,567,638,716,799,885,970,  //ISAMS_A1FI , ppm, scenario =0;
			325,337,353,369,388,412,437,463,488,509,525,537,545,549	//ISAMS_B1 , ppm, scenario =1;
										};
	j =0;
	for (i = 0; i<14; i++)
	{
		if (year<years[i]) break;

	}

	if (year == 2100) i =13 ;


	if (i==0)  // <1970
	{
		co2 = 282.23 + exp((year+doy/365.0)/51.35)*1.03*1.0e-15;
	}
	else
	{
		co2 = (co2scn[ctl][i]-co2scn[ctl][i-1])/(double(years[i]-years[i-1]))* (double(year)+double(doy)/365.-double(years[i-1]))+co2scn[ctl][i-1];
	}


	return co2;

}

