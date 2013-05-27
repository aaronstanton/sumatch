/* Copyright (c) University of Alberta, 2013.*/
/* All rights reserved.                       */
/* sumatchamp  :  $Date: Feb     2013- Last version Feb    2013  */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <time.h>

#ifndef MARK
#define MARK fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__);fflush(stderr);
#endif

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   						       ",
  " SMATCHAMP  Match trace-by-trace RMS amplitude of two input ",
  "            datasets                                        ",
  " 	   						       ",
  " Example of use:                                            ",
  " 	   						       ",
  " suplane taper=0 > tmp.su                                   ",
  " suplane taper=1 > tmp_ref.su                               ",
  " sumatchamp in=tmp.su ref=tmp_ref.su out=tmp_matched.su     ",
  "                                                            ",
NULL};
/* Credits:
 * Aaron Stanton.
 * Trace header fields accessed: ns, dt, ntr,
 * Last changes: Feb : 2013 
 */
/**************** end self doc ***********************************/

segy tr;
FILE* fpin; 
FILE* fpref;
FILE* fpout;

int main(int argc, char **argv)
{
  int verbose;
  time_t start,finish;
  double elapsed_time;
  int it,ih;
  float **datain  = 0;
  float **dataref = 0;
  float **dataout = 0;
  int nt, ntr, nh; 
  cwp_String in, out, ref;
  float trace_rms_in = 0;
  float trace_rms_ref = 0;
  float trace_rms_ratio = 0;
  float dt;

  if (verbose) fprintf(stderr,"*******SUMATCHAMP*********\n");
  /* Initialize */ 
  initargs(argc, argv);
  requestdoc(1);

  start=time(0);    
  if (!getparstring("in",&in)) err("in required."); 
  if (!getparstring("out",&out)) err("out required."); 
  if (!getparstring("ref",&ref)) err("ref required."); 
  if (!getparint("verbose", &verbose))  verbose =0;
  if (!getparint("ntr",&ntr)){ 
    ntr = 1000000;
    if (verbose) fprintf(stderr,"warning: ntr paramater not set; using ntr=%d\n",ntr);
  }
  /***********************************************************************
   begin reading the input data
  ***********************************************************************/  
  fpin = efopen(in, "r");
  ih = 0;	
  if(!fgettra(fpin,&tr,0)) err("can't read first trace of first input file");

  if (!tr.dt) err("dt header field must be set");
  if (!tr.ns) err("ns header field must be set");
  dt = ((float) tr.dt)/1000000.0;
  nt = (int) tr.ns;
  if (verbose) fprintf(stderr,"nt=%d\n",nt);

  fpref = efopen(ref, "r");
  ih = 0;	
  if(!fgettra(fpref,&tr,0)) err("can't read first trace of second input file");

  /* Allocate memory for data */
  nh = ntr;
  datain   = ealloc2float(nt,nh);

  do {
    memcpy((void *) datain[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
  } while (fgettr(fpin,&tr));
  nh=ih;
  if (verbose) fprintf(stderr,"done reading first input file: ih=%d\n",ih);

  fpref = efopen(ref, "r");
  ih = 0;	
  if(!fgettra(fpref,&tr,0)) err("can't read first trace of second input file");
	

  /* Allocate memory for data */
  dataref   = ealloc2float(nt,nh);

  do {
    memcpy((void *) dataref[ih],(const void *) tr.data,nt*sizeof(float));
    ih++;
    if (ih > nh) err("Number of traces > %d\n",nh); 
  } while (fgettr(fpref,&tr));
  if (verbose) fprintf(stderr,"done reading second input file: ih=%d\n",ih);
  if(ih != nh) err("input files must have same number of traces");
  fclose(fpref);
  /***********************************************************************
   end reading the input data
  ***********************************************************************/

  if (verbose) fprintf(stderr,"processing %d traces \n", nh);

  /* Allocate memory for data */
  dataout   = ealloc2float(nt,nh);

  for (ih=0;ih<nh;ih++){
    trace_rms_in = 0;
    trace_rms_ref = 0;
    for (it=0;it<nt;it++){
      trace_rms_in  = trace_rms_in  + datain[ih][it]*datain[ih][it];
      trace_rms_ref = trace_rms_ref + dataref[ih][it]*dataref[ih][it];
    }
    trace_rms_in  = sqrt(trace_rms_in);
    trace_rms_ref = sqrt(trace_rms_ref);
    if (abs(trace_rms_in) > 0.0001){
      trace_rms_ratio = trace_rms_ref/trace_rms_in;
    }
    else{
      trace_rms_ratio = trace_rms_ref/(trace_rms_in+0.0001);  
    }
    for (it=0;it<nt;it++){
      dataout[ih][it] = datain[ih][it]*trace_rms_ratio;
    }
  }
  /***********************************************************************
   outputting the data
  ***********************************************************************/
  rewind(fpin);
  fpout = efopen(out,"w");
  for (ih=0;ih<nh;ih++){ 
    fgettr(fpin,&tr);
    memcpy((void *) tr.data,(const void *) dataout[ih],nt*sizeof(float));
    fputtr(fpout,&tr);
  }
  fclose(fpout);

 /******** End of output **********/
 
 finish=time(0);
 elapsed_time=difftime(finish,start);
 if (verbose) fprintf(stderr,"Total time required: %f \n", elapsed_time);
 
 return EXIT_SUCCESS;
}





