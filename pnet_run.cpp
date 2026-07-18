#include "pnet_model.h"

void pnet_model::pnet_run(int wins)
{
	int i;
	char inputname[400];
	char  note[200];
	int argc;
	char** argv;

	FILE * f_input_load;	

	PathDefine(wins);

	sprintf(inputname,"%s%s%s",exePath,PathInput,"input.txt");

	f_input_load = fopen(inputname, "r");	

	if(f_input_load==NULL)	
	{
		printf("Unable to open input file %s!\n",inputname);
		exit(1) ;
	}	

	for (i=1; i<=3; i++)fgets(note, 200, f_input_load);

	fscanf(f_input_load, "%d\n", &modelmode);fgets(note, 200, f_input_load);	
	fscanf(f_input_load, "%d\n", &modeltype);fgets(note, 200, f_input_load);

	fclose(f_input_load);

	if (modelmode==1)  pnet_site();

//	if (modelmode==2)  pnet_region();

//	if (modelmode==3)  pnetii_region();

//	if (modelmode==4)  pnetcn_region();

//	if (modelmode==5)	epscor_region_merrimack_future();

//	if (modelmode == 6)  pnet_region(argc, argv);

//	if (modelmode == 7) pnet_region_test(argc, argv);

}


void  pnet_model::pnet_site()
{
	
	if (0==modeltype)  pnet_day(); // PnET-Day
	
	if (1==modeltype)  pnet_ii();// PnET-II

	if (2==modeltype)	pnet_cn(); // PnET-CN

	if (3==modeltype)	pnet_daily(); // PnET-Daily

	if (4==modeltype)	pnet_fun_corpse(); // PnET-FUN-CORPSE

}

void pnet_model::PathDefine(int wins)
{

	if (wins==0)sprintf(sep,"%s","/");	 //for linux path
	else sprintf(sep,"%s","\\");	 // for windows path

	sprintf(PathInput,"%sInput%s",sep,sep);	
	sprintf(PathLib,"%sLibrary%s",sep,sep);	
	sprintf(PathInter,"%sInter%s",sep,sep);	
	sprintf(PathRegion,"%sRegion%s",sep,sep);	
	sprintf(PathOutSite,"%sResult%s",sep,sep,sep);	
	sprintf(PathOutRegion,"%sResult%s",sep,sep,sep);	

//	sprintf(RegionClim,"%sInput%sRegionClim%s",sep,sep,sep);


}




