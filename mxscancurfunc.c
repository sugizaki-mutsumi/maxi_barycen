#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "fitsio.h"
#include "atFunctions.h"

#include "maxi_time.h"
//#include "maxi_orb.h"

#include "maxi_att.h"


#include "maxi_coord.h"
#include "maxi_coord_opt.h"

#include "maxi_colea.h"
//#include "maxi_cea.h"
#include "maxi_caldb.h"
#include "maxi_caldb_opt.h"

#include "maxi_gti.h"

#include "maxi_issanc.h"
//#include "obsjud.h"
//#include "obsjudex2.h"
//#include "maxi_fitsutil_obscur.h"

#include "obscol.h"

#include "gsc_caldb.h"

#include "maxi_piband.h"
#include "maxi_scan.h"


//#include "pil.h"
#include "headas.h"


#define MAX_NROWS 10000000

/*
int mxscancurfunc(char *instr, char *outfilename, 
		  double ra, double dec, double mjd_start, double mjd_stop, double dt, double coltha_max, 
		  char *attlistfile, char *obsjud_dname, char *sarjfnamelist, char *pmfnamelist, 
		  int gtiskip, char **gtifilelist){
*/
/* NEW_CODE_2017 Pass new calibration file paramter to mxscancurfunc */
int mxscancurfunc(char *instr, char *outfilename, 
      double ra, double dec, double mjd_start, double mjd_stop, double dt, double coltha_max,
      char *attlistfile, char *obscol_dname, char *sarjlistfname,  char *pmlistfname, 
      int gtiskip, char **gtifilename, char *leapsecfile, char *ssdockinfo, char *teldefgsc, 
      char *teldefssc, char *coleagsc, char *coleassc){
/* END NEW_CODE_2017 */   

  //set_toolname("mxscancur");
  //set_toolversion("2.0");

  int    status; 
  char   toolname[30], toolver[20], creator[30];
  TScanFitsHdrKeys fscankeys;

  int    instr_id;
  int    num_camera, cameraid;
  double colphi_max, sin_colphi_max, sin_coltha_max; //coltha_max, 

  //char   instr[20];
  //char   attlistfile[FILENAME_MAX];
  //char   outfilename[FILENAME_MAX];
  //char   gtifilename[MAX_NUM_CAM][FILENAME_MAX];
  //char   sarjlistfname[FILENAME_MAX], pmlistfname[FILENAME_MAX], obsjud_dname[FILENAME_MAX];

  //double dt, ra, dec;
  //double mjd_start, mjd_stop;
  char   *leapfile;
  char   *mxdatadir, *mxsoftdir;

  fitsfile *fptr;
  AtQuat qpar; 
  AtRotMat attrm;
  AtVect vec_ecs;
  AtVect vec_det;
  AtVect vec_maxi;
  
  double mxtime; //, mjd; 
  double col_theta, col_phi; //, det_theta, det_phi
  double *ontime, *eatime; 
  double interval;
  double area;

  double tstart, tstop;
  char   extname[80];

  int num_activecam; //gtiskip, 

  /// GTI
  TGti *gti[MAX_NUM_CAM];
  unsigned char status_gti;

  /// Obsjud
  //OBSJUD_DATAFILES obsjudFiles;
  //OBSJUD_DATA      *obsjudData[MAX_NUM_CAM];
  //int igeom, use_obsjud;
  //unsigned char status_obs, status_obsjud_camid[MAX_NUM_CAM];
  
  /// Obscol
  PDLDAT*   obscol_pdlDat;
  int use_obscol;
  unsigned char status_obs, status_obscol_camid[MAX_NUM_CAM];

  TIssAncSarjFileList sarjlist;
  TIssAncPmFileList   pmlist;
  TSsAlpha ssalpha;
  TSsBeta  ssbeta;
  int status_alpha, status_beta;


  /// Attitude
  ATTset *attset;
  TELDEF    **pteldef; 
  TColEAMap **pcolea;


  // SSINFO
  char *ssdockinfo_filename;
  GSC_SSDOCKINFO *ssdockinfo_gsc;
  int ssdockfov;
  

  long *firstrow, nrows; 
  TScancur *scancurdata[MAX_NUM_CAM];


  char *username = NULL;
  int  hdunum, hdutype;

  //status = PILGetString("instr", instr);
  if (strncmp(instr, "GSC", 3)==0 || strncmp(instr, "gsc", 3)==0){
    instr_id = INSTR_ID_GSC;
    num_camera = NUM_GSC_COUNTER;
    colphi_max = GSC_COLPHI_LIMIT;
  } else if (strncmp(instr, "SSC", 3)==0 || strncmp(instr, "ssc", 3)==0){
    instr_id = INSTR_ID_SSC;
    num_camera = NUM_SSC_COUNTER;
    colphi_max = SSC_COLPHI_LIMIT;
  } else {
    fprintf(stderr, "Error: INSTRUME has to be GSC or SSC\n");
    return -1;
  }
  
  //status = PILGetFname("outfname", outfilename);
  //status = PILGetReal("ra",     &ra);
  //status = PILGetReal("dec",    &dec);
  //status = PILGetReal("tstart", &mjd_start);
  //status = PILGetReal("tstop",  &mjd_stop);
  //status = PILGetReal("dt",     &dt);
  ////status = PILGetReal("colphi_max",  &colphi_max);  //colphi_max = 42.0;
  //status = PILGetReal("coltha_max",  &coltha_max);  //coltha_max = 5.0;

  //status = PILGetFname("attlist", attlistfile);
  //// Obsjud
  //status = PILGetFname("obsjud_dname", obsjud_dname);
  //status = PILGetFname("sarjfnamelist", sarjlistfname);
  //status = PILGetFname("pmfnamelist", pmlistfname);

  /// skip out except for GTI
  //status = PILGetInt("gtiskip", &gtiskip);

  /*
  status = PILGetFname("gtifile0",  gtifilename[0]);
  status = PILGetFname("gtifile1",  gtifilename[1]);
  if (instr_id==INSTR_ID_GSC){
    status = PILGetFname("gtifile2",  gtifilename[2]);
    status = PILGetFname("gtifile3",  gtifilename[3]);
    status = PILGetFname("gtifile4",  gtifilename[4]);
    status = PILGetFname("gtifile5",  gtifilename[5]);
    status = PILGetFname("gtifile6",  gtifilename[6]);
    status = PILGetFname("gtifile7",  gtifilename[7]);
    status = PILGetFname("gtifile8",  gtifilename[8]);
    status = PILGetFname("gtifile9",  gtifilename[9]);
    status = PILGetFname("gtifile10", gtifilename[10]);
    status = PILGetFname("gtifile11", gtifilename[11]);
  } else {
    status = PILPutFname("gtifile2",  "NONE");
    status = PILPutFname("gtifile3",  "NONE");
    status = PILPutFname("gtifile4",  "NONE");
    status = PILPutFname("gtifile5",  "NONE");
    status = PILPutFname("gtifile6",  "NONE");
    status = PILPutFname("gtifile7",  "NONE");
    status = PILPutFname("gtifile8",  "NONE");
    status = PILPutFname("gtifile9",  "NONE");
    status = PILPutFname("gtifile10", "NONE");
    status = PILPutFname("gtifile11", "NONE");
  }
  */

  /// init TELDEF
  /* NEW_CODE_2017 Add support for teldef file list */
  pteldef = malloc(sizeof(**pteldef)*num_camera);
  if (instr_id == INSTR_ID_GSC){
    init_allgsc_teldef(pteldef,teldefgsc);
  } else {
    init_allssc_teldef(pteldef,teldefssc);
  }
  /* END NEW_CODE_2017 */
  
  /// init COLEA
  /* NEW_CODE_2017 Add support for coleagsc and coleassc file names */
  pcolea = malloc(sizeof(**pcolea)*num_camera);
  /* NEW_CODE_2017 Add support for coleagsc and coleassc file names */
  if (instr_id == INSTR_ID_GSC){
    init_allgsc_colea(pcolea,coleagsc);
  } else {
    init_allssc_colea(pcolea,coleassc);
  }
  /* END NEW_CODE_2017 */
  printf("gtifilename[0]=%s\n", gtifilename[0]);
  
  /// init GTI
  for(cameraid=0; cameraid<num_camera; cameraid++){
    if (strcmp(gtifilename[cameraid],"NONE")==0 || strcmp(gtifilename[cameraid],"none")==0
	|| strcmp(gtifilename[cameraid],"DEF")==0 || strcmp(gtifilename[cameraid],"def")==0
	){
      gti[cameraid]=NULL;
    } else {
      //gti[cameraid] = maxi_gti_init(gtifilename[cameraid]);
      gti[cameraid] = maxi_gti_readfitsfile(gtifilename[cameraid], "STDGTI");
    }
    if(gti[cameraid]==NULL) {
      fprintf(stderr, "Warning: GSC_%X GTI file is not set.  It is ignored\n", cameraid);
    }
  }
  
  /// init Time
  /* NEW_CODE_2017 Addeding refdata and filename support for leapfile */ 
  if (strcasecmp(leapsecfile,"CALDB") == 0) {
    leapfile = maxi_caldb_find_leapfile("AUTO");
  } else if (strcasecmp(leapsecfile,"REFDATA") == 0) {
    char * refdata = getenv("LHEA_DATA");
    leapfile = strcat(refdata,"/leapsec.fits");
  } else {
    leapfile = leapsecfile;
  }
  headas_printf("LEAPFILE = %s\n", leapfile);
  mission_time_init(leapfile);
  /* END NEW_CODE_2017 */

  headas_printf("Start Time (MJD) = %lf\n", mjd_start);
  headas_printf("Stop Time (MJD) = %lf\n",  mjd_stop);
  
  tstart = mjd2maxi(mjd_start);
  tstop  = mjd2maxi(mjd_stop);
  

  /// get MXDATA, MXSOFT
  mxdatadir = getenv("MXDATA");
  mxsoftdir = getenv("MXSOFT");

  /// init ATT
  if( strcasecmp(attlistfile,"DEF")==0 ) {
    if (mxdatadir == NULL) {
      headas_printf("Error: MXDATA environment variable is not set up.\n");
      return -1;
    }
    sprintf(attlistfile, "%s/trend/attlist.fits", mxdatadir);
  }

  attset = maxi_att_init(attlistfile);
  headas_printf("ATT:TSTART = %lf\n", attset->attfileList->tstart);
  headas_printf("ATT:TSTOP  = %lf\n", attset->attfileList->tstop);
  interval = tstop-tstart;


  /// init Obsjud Data
  /// init Obscol Data
  /* NEW_CODE_2017 Switch to using strcasecmp instead of multiple strcmp checks, 
                   Using REFDATA as default value for obscol_dname and searching 
                   the enviroment varible LHEA_DATA  */
  if (strcasecmp(obscol_dname, "NONE")==0){
    use_obscol = 0;
  } else {

    use_obscol = 1;
    if (strcasecmp(obscol_dname,"REFDATA") == 0){
      obscol_dname = getenv("LHEA_DATA");
    }   
   /* END NEW_CODE_2017 */
    printf("obscol_dname = %s\n", obscol_dname);
    obscol_pdlDat = initPdlDat(obscol_dname);
    printf("just called initPdlDat\n");

    //igeom = 0;

    /*
    if( strcmp(obsjud_dname,"DEF")==0 ) {
      if (mxsoftdir == NULL) {
	headas_printf("Error: MXSOFT environment variable is not set up.\n");
	return -1;
      }
      sprintf(obsjud_dname, "%s/refdata/maxi/obsjud", mxsoftdir);
    }

    
    status = obsjud_get_data_file_names(obsjud_dname, &obsjudFiles);
    printf("obsjud_get_data_file_names status = %d\n", status);
    if (status==0){
      return status;
    }


    for(cameraid=0; cameraid<num_camera; cameraid++){
      obsjudData[cameraid] = malloc(sizeof(*obsjudData[cameraid]));
      status = obsjud_initialize(&obsjudFiles, obsjudData[cameraid]);
      printf("obsjud_initialize camid=%d status = %d\n", cameraid, status);
      if (status==0){
	return status;
      }
    }
    */



    /// init IssAncSarj
    if( strcmp(sarjlistfname,"DEF")==0 ) {
      if (mxdatadir == NULL) {
	headas_printf("Error: MXDATA environment variable is not set up.\n");
	return -1;
      }
      sprintf(sarjlistfname, "%s/trend/isalist.fits", mxdatadir);
    }
    status = maxi_issanc_sarjlist_init(&sarjlist, sarjlistfname);
    printf("sarjlist.tstart = %lf\n", sarjlist.tstart);
    printf("sarjlist.tstop  = %lf\n", sarjlist.tstop);
    printf("sarjlist.nfile  = %d\n",  sarjlist.nfile);
    printf("sarjlist.current_idx  = %d\n",  sarjlist.current_idx);
    
    /// init IssAncPm
    if( strcmp(pmlistfname,"DEF")==0 ) {
      if (mxdatadir == NULL) {
	headas_printf("Error: MXDATA environment variable is not set up.\n");
	return -1;
      }
      sprintf(pmlistfname, "%s/trend/isplist.fits", mxdatadir);
    }
    status = maxi_issanc_pmlist_init(&pmlist, pmlistfname);
    printf("pmlist.tstart = %lf\n", pmlist.tstart);
    printf("pmlist.tstop  = %lf\n", pmlist.tstop);
    printf("pmlist.nfile  = %d\n",  pmlist.nfile);
    printf("pmlist.current_idx  = %d\n",  pmlist.current_idx);
  }


  /* NEW_CODE_2017 Adding file name option for space shuttle docking info file */
  //// SSDOCK info 
  if (strcasecmp(ssdockinfo,"CALDB") == 0) {
    ssdockinfo_filename = maxi_caldb_find("GSC", "SSDOCKINFO", "CALDB");
    ssdockinfo_gsc = init_gsc_ssdockinfo(ssdockinfo_filename);
    headas_printf("SSDOCKINFO filename=%s\n",  ssdockinfo_filename);
    free(ssdockinfo_filename);
  } else {
    ssdockinfo_filename = ssdockinfo;
    ssdockinfo_gsc = init_gsc_ssdockinfo(ssdockinfo_filename);
    headas_printf("SSDOCKINFO filename=%s\n",  ssdockinfo_filename);
  }
  /* END NEW_CODE_2017 */
  
  //printf("### status=%d\n", status);
  headas_printf("StartTime = %lf\n", tstart);
  headas_printf("EndTime   = %lf\n", tstop);
  headas_printf("Delta_T   = %lf\n", dt);
  headas_printf("Interval  = %lf [sec] = %.3lf [hr]\n", interval, interval/3600.0);
  headas_printf("Source (ra,dec) = (%lf,%lf)\n", ra, dec);
  headas_printf("\n");


  ////
  ////  initialize Header Key parameters
  ////
  get_toolname(toolname);
  get_toolversion(toolver);
  sprintf(creator, "%s_%s", toolname, toolver);
  /* NEW_CODE_2017 Make sure user keyword is set to task name */
  username = "mxscancurfunc";
  /* END NEW_CODE_2017 */


  /// init cameraid param
  ontime = malloc(sizeof(*ontime)*num_camera);
  eatime = malloc(sizeof(*eatime)*num_camera);
  firstrow = malloc(sizeof(*firstrow)*num_camera);
  
  ///// initialize data
  for(cameraid=0; cameraid<num_camera; cameraid++) {
    ontime[cameraid]=0.0;
    eatime[cameraid]=0.0;
    firstrow[cameraid]=1;
    scancurdata[cameraid] = maxi_scancur_init(MAX_NROWS);
  }


  /// create output file
  headas_printf("output filename = %s\n", outfilename);  
  headas_clobberfile(outfilename);
  status=0;
  if ( fits_create_file(&fptr, outfilename, &status) ) 
    return status;
  
  for(cameraid=0; cameraid<num_camera; cameraid++) {
    sprintf(extname, "SCANCUR_CAMID%X", cameraid);
    maxi_scancurhdu_inittab(fptr, extname);
  }



  status_gti = 0;
  status_obs = 0;

  // vec_ecs
  atPolDegToVect(1.0, ra, dec, vec_ecs);
  sin_coltha_max = sin(coltha_max*DEG2RAD);
  sin_colphi_max = sin(colphi_max*DEG2RAD);

  printf("sin_coltha_max=%lf, sin_colphi_max=%lf\n", sin_coltha_max, sin_colphi_max);

  mxtime=tstart; 
  maxi_att_qpar(attset, mxtime, qpar);
  atQuatToRM(qpar, attrm);


  for(mxtime=tstart; mxtime<tstop; mxtime+=dt) {

    /// check gtiskip or not
    if(gtiskip==1){
      num_activecam=0;
      for(cameraid=0; cameraid<num_camera; cameraid++){
	if(maxi_gti_getstat(gti[cameraid], mxtime)==1){
	  num_activecam++;
	}
      }
      if(num_activecam==0){
	continue;
      }
    }


    /// get ATT RM
    //maxi_att_qpar(attset, tstart, qpar);
    maxi_att_qpar(attset, mxtime, qpar);
    atQuatToRM(qpar, attrm);
    //atRotVect(attrm, vec_ecs, vec_maxi);
    maxi_vececs2vecmaxi(attrm, vec_ecs, vec_maxi);

    /// get Paddle Angle
    if (use_obscol==1){
      status_alpha = maxi_issanc_sarjlist_get_ssalpha(&sarjlist, mxtime, ssalpha);
      status_beta  = maxi_issanc_pmlist_get_ssbeta(&pmlist, mxtime, ssbeta);    
    }

    for(cameraid=0; cameraid<num_camera; cameraid++){

      /// GTI
      if(gti[cameraid]!=NULL){
	status_gti = (unsigned char)maxi_gti_getstat(gti[cameraid], mxtime);
      } else {
	status_gti = -1;
      }
      if(gtiskip==1 && status_gti!=1){
	continue;
      }

      /// col angle
      maxi_vecmaxi2vecdet(pteldef[cameraid], vec_maxi, vec_det);
      if( vec_det[2]<0.0 ) {
	continue;
      }
      if( !( fabs(vec_det[1])<sin_coltha_max && fabs(vec_det[0])<sin_colphi_max ) ) {
	continue;
      }
      maxi_vecdet2colsp(vec_det, &col_theta, &col_phi);

      //printf("col_theta=%lf, col_phi=%lf\n", col_theta, col_phi);  
      if( fabs(col_theta)>coltha_max || fabs(col_phi)>colphi_max ) {
        continue;
      }

      //// shuttle dock 
      ssdockfov = gsc_caldb_ssdock_checkfov(cameraid, mxtime, col_theta, col_phi, ssdockinfo_gsc);

      //// col area
      area = get_colea(pcolea[cameraid], col_theta, col_phi);

      ///// obsjud //////
      status_obs = 0;
      if (use_obscol==1){
	if ( status_alpha==1 && status_beta==1 ){
	  //status_obsjud = (unsigned char)obsjud_radec(&obsjudData, qpar, pteldef[cameraid], igeom, ra, dec, ssalpha, ssbeta);

	  status_obscol_camid[cameraid] = (unsigned char)abcolAll(obscol_pdlDat, cameraid, col_phi, ssalpha, ssbeta);
	  } else {
	    status_obscol_camid[cameraid] = 20;
	}
	status_obs = status_obscol_camid[cameraid];
      } 


      /// store data
      nrows = scancurdata[cameraid]->nrows;

      scancurdata[cameraid]->time[nrows] = mxtime;
      scancurdata[cameraid]->coltha[nrows] = col_theta;
      scancurdata[cameraid]->colphi[nrows] = col_phi;
      scancurdata[cameraid]->colea[nrows]  = area;

      scancurdata[cameraid]->gtiflag[nrows] = status_gti;
      scancurdata[cameraid]->obsflag[nrows] = status_obs + ssdockfov*100;

      scancurdata[cameraid]->nrows+=1;


      /// ontime and eatime
      if(status_gti==1 ){ //&& status_obsjud==0) {
	scancurdata[cameraid]->ontime+=dt;
	scancurdata[cameraid]->eatime+=(area*dt);
      }


      /// flush scancur to file
      if(scancurdata[cameraid]->nrows>=MAX_NROWS){
	hdunum = cameraid+2;
	fits_movabs_hdu(fptr, hdunum, &hdutype, &status);
	maxi_scancurhdu_append(fptr, firstrow[cameraid], scancurdata[cameraid]);

	ontime[cameraid] += scancurdata[cameraid]->ontime;
	eatime[cameraid] += scancurdata[cameraid]->eatime;
	firstrow[cameraid] += scancurdata[cameraid]->nrows;
	maxi_scancur_clear(scancurdata[cameraid]);
      }


    }
  }
 
  /// flush scancur to file
  for(cameraid=0; cameraid<num_camera; cameraid++){
    hdunum = cameraid+2;
    fits_movabs_hdu(fptr, hdunum, &hdutype, &status);
    maxi_scancurhdu_append(fptr, firstrow[cameraid], scancurdata[cameraid]);
    ontime[cameraid] += scancurdata[cameraid]->ontime;
    eatime[cameraid] += scancurdata[cameraid]->eatime;
    maxi_scancur_clear(scancurdata[cameraid]);
  }

  maxi_att_close(attset);


  //maxi_orb_close();
  if(instr_id == INSTR_ID_GSC){
    free_allgsc_teldef(pteldef);
    free_allgsc_colea(pcolea);
  } else {
    free_allssc_teldef(pteldef);
    free_allssc_colea(pcolea);
  }
  free(pteldef);
  free(pcolea);


  // default key parameter
  //maxi_obscur_fitsHdrKeys_init(&fobscurkeys, instr, creator, username, tstart, tstop, dt, 0.0, ra, dec);
  maxi_scan_fitsHdrKeys_init(&fscankeys, instr, creator, username, tstart, tstop, dt, ra, dec);

  for(cameraid=0; cameraid<num_camera; cameraid++){
    // update instr keyword
    if(instr_id==INSTR_ID_GSC) {
      sprintf(fscankeys.instrume, "GSC_%X", cameraid);
    } else {
      if(cameraid==0) strcpy(fscankeys.instrume, "SSC_H");
      else strcpy(fscankeys.instrume, "SSC_Z");
    }
    //fobscurkeys.exposure = exposure[cameraid];

    hdunum = cameraid+2;
    fits_movabs_hdu(fptr, hdunum, &hdutype, &status);

    maxi_scan_fitsHdrKeys_update(fptr, &fscankeys);
    HDpar_stamp(fptr, 0, &status);
    fits_write_chksum(fptr, &status);
    fits_write_date(fptr, &status);

    
  }

  fits_close_file(fptr, &status);


  /// free file pinter /// this should not be done
  //if(leapfile) free(leapfile);
  //if(mxdatadir) free(mxdatadir);
  //if(mxsoftdir) free(mxsoftdir);

  /// free memory allocation 
  if(ontime) free(ontime);
  if(eatime) free(eatime);
  if(firstrow) free(firstrow);


  for(cameraid=0; cameraid<num_camera; cameraid++){
   if(gti[cameraid]) maxi_gti_free(gti[cameraid]);
    //if(obscolData[cameraid]) free(obscolData[cameraid]);
    if(scancurdata[cameraid]) maxi_scancur_free(scancurdata[cameraid]);
  }

  // if(obscol_pdlDat) freePdlDat(obscol_pdlDat);

  //return status;

  return 0;
}
  
