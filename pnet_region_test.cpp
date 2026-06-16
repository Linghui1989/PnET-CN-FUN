#include "pnet_model.h"
#include <iostream>
#include <string>
#include <netcdf>
#include <stdio.h>
#include <numeric>
#include <mpi.h>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;



void pnet_model::pnet_region_test(int argc, char** argv)
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

	float ISA0[NLAT][NLON], ISA1[NLAT][NLON], ISA2[NLAT][NLON];

	vector<size_t> start, count;
	start = { 0, 0 };
	count = { NLAT, NLON };

	sprintf(inputname, "%s%s%s", exePath, PathInput, "ISA_1900.nc");
	datafile.open(inputname, NcFile::read);
	var = datafile.getVar("ISA");
	var.getVar(start, count, ISA0);
	datafile.close();


	sprintf(inputname, "%s%s%s", exePath, PathInput, "ISA_2020.nc");
	datafile.open(inputname, NcFile::read);
	var = datafile.getVar("ISA");
	var.getVar(start, count, ISA1);
	datafile.close();

	sprintf(inputname, "%s%s%s", exePath, PathInput, "ISA_2050.nc");
	datafile.open(inputname, NcFile::read);
	var = datafile.getVar("ISA");
	var.getVar(start, count, ISA2);
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
	ReadInput(&site,&veg,&share0,inputname);	// read in input data from input files
	initvars(&site,&veg,&share0);
	


	// ============disturbance=============================================
	Disturbance(&site, &veg, &share0, 0);  // to allocate memory for disturbances
	Disturbance(&site, &veg, &share0, 2);  // default disturbance as HB

//	printf("set up model \n");
	
	for (int lat = lat_start; lat <= lat_end; lat++)
	{
		site.Lat = lats[lat];  // updating latitude

		for (int lon = 0; lon < NLON; lon++) // HB
		{

			float nl = nlcd[lat][lon];
			
			if (!isnan(nl))
			{
				nl = 1000 * lat + lon;
				share_struct share;
				share = share0;
				initvars(&site, &veg, &share);
				init_out(ryrs, &out);

				

				site.Lat = lats[lat];  // updating latitude
//				printf("executing lat is %d lon is %d \n", lat, lon);
//				printf("executing lat is %d \n", site.Lat);

				nsteps = climhis.length;

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

				share.ISA = ISA0[lat][lon];
				share.ISA1 = (ISA1[lat][lon] - ISA0[lat][lon]);
				share.ISA2 = ISA2[lat][lon];


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

				// ===================       daily running      ==================
				//printf("executing daily run length is %d \n", ndays);
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
						AllocateYr(&site, &veg, &climday, rstep + 1, &share, CN_Mode); // Note the rstep+1, not rstep
						storeoutput(&veg, &share, &out, rstep + 1, &ystep, 1);//Linghui Meng 06302021
						YearInit(&share);
					}
				}

				//define the name of annual result
				sprintf(fileout_yr, "%s%f%s", pnetresult, nl, "_Output_annual.csv");
				WriteoutYr_netcdf(&climday, &out, fileout_yr);
				//printf("write out\n");

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

