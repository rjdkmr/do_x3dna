/*
 * This file is part of do_x3dna
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2014-2018  Rajendra Kumar
 *
 * do_x3dna uses 3DNA package (http://x3dna.org).
 * Please cite the original publication of the 3DNA package:
 * Xiang-Jun Lu & Wilma K. Olson (2003)
 * 3DNA: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures
 * Nucleic Acids Res. 31(17), 5108-21.
 *
 * do_x3dna is a free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * do_x3dna is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with do_x3dna.  If not, see <http://www.gnu.org/licenses/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "gromacs/commandline/viewit.h"
#include "gromacs/utility/futil.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/cmdlineinit.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/index.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"


#include "../ExtractData.h"
#include "do_x3dna.h"

void CopyRightMsg() {

    const char *copyright[] = {
            "                                                                        ",
            "                  :-)  do_x3dna (-:                                     ",
            "                                                                        ",
            "               Author: Rajendra Kumar                                   ",
            "                                                                        ",
            "         Copyright (C) 2014-2017  Rajendra Kumar                        ",
            "                                                                        ",
            "do_x3dna uses 3DNA package (http://x3dna.org).                          ",
            "Please cite the original publication of the 3DNA package:               ",
            "Xiang-Jun Lu & Wilma K. Olson (2003)                                    ",
            "3DNA: a software package for the analysis, rebuilding and visualization ",
            "of three-dimensional nucleic acid structures.                           ",
            "Nucleic Acids Res. 31(17), 5108-21.                                     ",
            "                                                                        ",
            "do_x3dna is a free software: you can redistribute it and/or modify      ",
            "it under the terms of the GNU General Public License as published by    ",
            "the Free Software Foundation, either version 3 of the License, or       ",
            "(at your option) any later version.                                     ",
            "                                                                        ",
            "do_x3dna is distributed in the hope that it will be useful,             ",
            "but WITHOUT ANY WARRANTY; without even the implied warranty of          ",
            "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           ",
            "GNU General Public License for more details.                            ",
            "                                                                        ",
            "You should have received a copy of the GNU General Public License       ",
            "along with do_x3dna.  If not, see <http://www.gnu.org/licenses/>.       ",
            "                                                                        ",
            "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS     ",
            "\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT     ",
            "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR   ",
            "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT    ",
            "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,   ",
            "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED",
            "TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR  ",
            "PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF  ",
            "LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING    ",
            "NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      ",
            "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            ",
            "                                                                        ",
            "                           :-)  do_x3dna (-:                            ",
            "                                                                        ",
            "                                                                        "
    };
    int i = 0;
    char *str;
    for(i=0; i<42; i++) {
        str = strdup(copyright[i]);
        fprintf(stderr,"%s\n", str);
    }
}



int calculate_avg_std(double *data, int n, double filter, double *avg, double *std)	{
	int i=0, count=1;
	double sum=0, sum_sq=0;
	for(i=0;i<n;i++)	{
		if(data[i]==filter)
			continue;
		sum += data[i];
		sum_sq += (data[i]*data[i]);

		count++;
	}

	*avg = sum/count;
	*std = sqrt(  (sum_sq - (sum*sum)/count)/count  );

	return 0;
}

int local_base_pair_step_out(gmx_bool bAvg, const char *fn_avg_out, const char *property[], const char *ComName[], int nframe, int **max_bps, int max_num_bps, char *fn_inp_base_pair,char *fn_inp_lbps, const gmx_output_env_t *oenv)	{
	int i,j,m, n;
	double dum = 999;
	FILE *f_inp_base_pair, *f_inp_lbps;
	real *time;
	int timer1=0, timer2=0;
	char *buff1, *buff2;
	char **lines1, **lines2, **time_temp;
	int num1=0, num2=0, *tmp_bp1, *tmp_bp2;
	double *lbps;

	t_local_bp_step_param *bp_param;

	snew(time,nframe);
	snew(bp_param,max_num_bps);
	for(n=0;n<(max_num_bps);n++)
	{
		bp_param[n].bp11 = max_bps[n][0];
		bp_param[n].bp12 = max_bps[n][1];
		bp_param[n].bp21 = max_bps[n][2];
		bp_param[n].bp22 = max_bps[n][3];

		snew(bp_param[n].property,6);
		for(i=0;i<6;i++)	{
			snew(bp_param[n].property[i],nframe);
			for(j=0;j<nframe;j++)
				bp_param[n].property[i][j] = dum;
		}
	}


	//Reading base-pairs and parameters from X3DNA output
	fprintf(stderr,"\nReading parameters from %s file.....", fn_inp_lbps);
	f_inp_base_pair = fopen(fn_inp_base_pair, "r");
	f_inp_lbps = fopen(fn_inp_lbps, "r");



	//Initialization First Frame
	while(1)	{
		buff1 = get_line(f_inp_base_pair);
		if(strstr(buff1,"# Time")!=NULL)
			break;
	}

	while(1)	{
		buff2 = get_line(f_inp_lbps);

		if((strstr(buff2,"# Time")!=NULL) && (buff2 != NULL))	{
			break;
		}
	}


	time_temp = split_by_char(buff2, "=", NULL);
	time[0] = strtof(time_temp[1],NULL);

	//Extracting data from all frames
	while(1)	{

		lines1 = get_block_lines(f_inp_base_pair,"# Time",&num1);
		lines2 = get_block_lines(f_inp_lbps,"# Time",&num2);

		for(i=0;i<num1-2;i++)	{

			  if(!is_first_numeric(lines1[i]))
				  continue;

			  if(strstr(lines2[i],"----")!=NULL)
				  continue;

			  tmp_bp1 = extract_column_integer(lines1[i],1,2);
			  tmp_bp2 = extract_column_integer(lines1[i+1],1,2);
			  lbps = extract_column_double(lines2[i],1,6);

			  //printf("%d % d%d %d %f %f %f %f %f %f\n",tmp_bp1[0],tmp_bp1[1], tmp_bp2[0],tmp_bp2[1], lbps[0], lbps[1], lbps[2], lbps[3], lbps[4], lbps[5]);

			  for(j=0;j<max_num_bps;j++)	{

				  if	(
							(
								((tmp_bp1[0] == max_bps[j][0] && tmp_bp1[1] == max_bps[j][1]) || (tmp_bp1[1] == max_bps[j][0] && tmp_bp1[0] == max_bps[j][1]))
							&&
								((tmp_bp2[0] == max_bps[j][2] && tmp_bp2[1] == max_bps[j][3]) || (tmp_bp2[1] == max_bps[j][2] && tmp_bp2[0] == max_bps[j][3]))
							)
						||
							(
								((tmp_bp2[0] == max_bps[j][0] && tmp_bp2[1] == max_bps[j][1]) || (tmp_bp2[1] == max_bps[j][0] && tmp_bp2[0] == max_bps[j][1]))
							&&
								((tmp_bp1[0] == max_bps[j][2] && tmp_bp1[1] == max_bps[j][3]) || (tmp_bp1[1] == max_bps[j][2] && tmp_bp1[0] == max_bps[j][3]))
							 )
					  	 )
				  {
					  for(m=0;m<6;m++)	{
						  bp_param[j].property[m][timer1] = lbps[m];

						  //printf("%f", bp_param[j].property[m][timer1] );
					  }
					  //printf("\n");
				  }
			  }

		}

		timer1++;
		timer2++;


		time_temp = split_by_char(lines1[num1], "=", NULL);
		time[timer1] = strtof(time_temp[1],NULL);


		//Processing for last frame
		if(timer1==nframe-1)	{
			lines1 = get_all_lines(f_inp_base_pair,&num1);
			lines2 = get_all_lines(f_inp_lbps,&num2);

			for(i=0;i<num1-1;i++)	{
				  if(!is_first_numeric(lines1[i]))
					  continue;

				  if(strstr(lines2[i],"----")!=NULL)
					  continue;

				  tmp_bp1 = extract_column_integer(lines1[i],1,2);
				  tmp_bp2 = extract_column_integer(lines1[i+1],1,2);
				  lbps = extract_column_double(lines2[i],1,6);

				  for(j=0;j<max_num_bps;j++)	{
					  if	(
								(
									((tmp_bp1[0] == max_bps[j][0] && tmp_bp1[1] == max_bps[j][1]) || (tmp_bp1[1] == max_bps[j][0] && tmp_bp1[0] == max_bps[j][1]))
								&&
									((tmp_bp2[0] == max_bps[j][2] && tmp_bp2[1] == max_bps[j][3]) || (tmp_bp2[1] == max_bps[j][2] && tmp_bp2[0] == max_bps[j][3]))
								)
							||
								(
									((tmp_bp2[0] == max_bps[j][0] && tmp_bp2[1] == max_bps[j][1]) || (tmp_bp2[1] == max_bps[j][0] && tmp_bp2[0] == max_bps[j][1]))
								&&
									((tmp_bp1[0] == max_bps[j][2] && tmp_bp1[1] == max_bps[j][3]) || (tmp_bp1[1] == max_bps[j][2] && tmp_bp1[0] == max_bps[j][3]))
								 )
						  	 )
					  {
						  //printf("%d\t",timer1);
						  for(m=0;m<6;m++)	{
							  bp_param[j].property[m][timer1] = lbps[m];
							 //printf("%f", bp_param[j].property[m][timer1] );
						  }
						  //printf("\n");
					  }
				  }
			}

			break;
		}

	}


	fclose(f_inp_base_pair);
	fclose(f_inp_lbps);
	// Output Files Generation
	FILE *fOut[6], *fLBPM=NULL;
	char **leg,**setname;
	char *fnOut[6], *fnLBPM;
	//char *property[6]=	{ "Shear", "Stretch", "Stagger", "Buckle", "Propeller", "Opening" };
	char title[256];

	snew(leg,max_num_bps);

	for(i=0;i<max_num_bps;i++)	{
		snew(leg[i],256);
		sprintf(leg[i], "%d=%d/%d=%d",max_bps[i][0],max_bps[i][1],max_bps[i][2],max_bps[i][3]);
	}

	for(m=0;m<6;m++)	{
		snew(fnOut[m],256);
		sprintf(fnOut[m],"%s_%s.xvg",property[m],ComName[0]);
		sprintf(title,"Local base pair parameter: %s",property[m]);
		fOut[m] = xvgropen(fnOut[m],title,output_env_get_time_label(oenv),property[m], oenv);
		xvgr_legend(fOut[m],max_num_bps,(const char **) leg,oenv);
	}

	for(n=0;n<nframe;n++)	{

		for(m=0;m<6;m++)
			fprintf(fOut[m],"\n%12.7f",time[n]);

		for(i=0;i<max_num_bps;i++)
			for(m=0;m<6;m++)
				fprintf(fOut[m], "%10.3f",bp_param[i].property[m][n]);
	}

	fprintf(stderr,"\rFinished reading parameters from %s file\n", fn_inp_lbps);
	fprintf(stderr,"Following output files generated:\n");
	for(m=0;m<6;m++)
		fprintf(stderr,"\t%s\n",fnOut[m]);

	for(m=0;m<6;m++)
		xvgrclose(fOut[m]);

	if(bAvg)	{
		fprintf(stderr,"\rCalculating Average and Standard Deviation of local base pair parameters...");

		snew(fnLBPM,256);
		sprintf(fnLBPM,"%s_%s.xvg",fn_avg_out, ComName[0]);
		fLBPM = xvgropen_type(fnLBPM,"Local base pair parameter","Base Pair","Parameters Value",exvggtXYDY, oenv);
		for(j=0;j<max_num_bps;j++)
			fprintf(fLBPM,"# %d\t%d=%d/%d=%d\n",j+1,max_bps[j][0],max_bps[j][1],max_bps[j][2],max_bps[j][3]);
		for(i=0;i<max_num_bps;i++)	{
			snew(bp_param[i].avg,6);
			snew(bp_param[i].std,6);
			for(j=0;j<6;j++)	{
				calculate_avg_std(bp_param[i].property[j], nframe , dum, &bp_param[i].avg[j], &bp_param[i].std[j]);
			}
		}

		snew(setname,6);
		for(i=0;i<6;i++)
			snew(setname[i],256);
		snew(setname[1],256);
		sprintf(setname[0], "%s", property[0]);
		sprintf(setname[1], "%s", property[1]);
		sprintf(setname[2], "%s", property[2]);
		sprintf(setname[3], "%s", property[3]);
		sprintf(setname[4], "%s", property[4]);
		sprintf(setname[5], "%s", property[5]);
		xvgr_new_dataset(fLBPM,0,6,(const char **) setname,oenv);

		for(i=0;i<6;i++)	{
			for(j=0;j<max_num_bps;j++)
				fprintf(fLBPM,"%d %10.3f %10.3f\n",(j+1),bp_param[j].avg[i],bp_param[j].std[i]);
			fprintf(fLBPM,"&\n\n");
		}

		fprintf(stderr,"\rCalculated Average and Standard Deviation of local base pair parameters. Output file: %s\n",fnLBPM);
	}

	xvgrclose(fLBPM);
	sfree(bp_param);
	free(buff1);
	free(buff2);
	free(lines1);
	free(lines2);
	free(time_temp);
	free(time);
	return 0;
}

int local_base_pair_out(gmx_bool bAvg, const char *ComName[], int nframe, int **max_bp, int max_num_bp, char *fn_inp_base_pair,char *fn_inp_lbp, const gmx_output_env_t *oenv)	{
	int i,j,m, n;
	double dum = 999;
	FILE *f_inp_base_pair, *f_inp_lbp;
	real *time;
	int timer1=0, timer2=0;
	char *buff1, *buff2;
	char **lines1, **lines2, **time_temp;
	int num1=0, num2=0, *tmp_bp;
	double *lbp;

	t_local_bp_param *bp_param;

	snew(time,nframe);
	snew(bp_param,max_num_bp);
	for(n=0;n<max_num_bp;n++)
	{
		bp_param[n].bp1 = max_bp[n][0];
		bp_param[n].bp2 = max_bp[n][1];
		snew(bp_param[n].property,6);
		for(i=0;i<6;i++)	{
			snew(bp_param[n].property[i],nframe);
			for(j=0;j<nframe;j++)
				bp_param[n].property[i][j] = dum;
		}
	}


	//Reading base-pairs and parameters from X3DNA output
	fprintf(stderr,"\nReading parameters from %s file.....", fn_inp_lbp);
	f_inp_base_pair = fopen(fn_inp_base_pair, "r");
	f_inp_lbp = fopen(fn_inp_lbp, "r");



	//Initialization First Frame
	while(1)	{
		buff1 = get_line(f_inp_base_pair);
		if(strstr(buff1,"# Time")!=NULL)
			break;
	}

	while(1)	{
		buff2 = get_line(f_inp_lbp);

		if((strstr(buff2,"# Time")!=NULL) && (buff2 != NULL))	{
			break;
		}
	}


	time_temp = split_by_char(buff2, "=", NULL);
	time[0] = strtof(time_temp[1],NULL);


	//Extracting data from all frames
	while(1)	{

		lines1 = get_block_lines(f_inp_base_pair,"# Time",&num1);
		lines2 = get_block_lines(f_inp_lbp,"# Time",&num2);

		for(i=0;i<num1;i++)	{

			  if(!is_first_numeric(lines1[i]))
				  continue;

			  tmp_bp = extract_column_integer(lines1[i],1,2);
			  lbp = extract_column_double(lines2[i],1,6);

			  //printf("%d %d %f %f %f %f %f %f\n",tmp_bp[0],tmp_bp[1], lbp[0], lbp[1], lbp[2], lbp[3], lbp[4], lbp[5]);

			  for(j=0;j<max_num_bp;j++)	{
				  if( ((max_bp[j][0]==tmp_bp[0]) && (max_bp[j][1]==tmp_bp[1]) ) || ( (max_bp[j][0]==tmp_bp[1]) && (max_bp[j][1]==tmp_bp[0]) ))	{
					  //printf("%d\t",timer1);
					  for(m=0;m<6;m++)	{
						  bp_param[j].property[m][timer1] = lbp[m];

						  //printf("%f", bp_param[j].property[m][timer1] );
					  }
					  //printf("\n");
				  }
			  }
		}


		timer1++;
		timer2++;


		time_temp = split_by_char(lines1[num1], "=", NULL);
		time[timer1] = strtof(time_temp[1],NULL);


		//Processing for last frame
		if(timer1==nframe-1)	{
			lines1 = get_all_lines(f_inp_base_pair,&num1);
			lines2 = get_all_lines(f_inp_lbp,&num2);

			for(i=0;i<num1;i++)	{
				  if(!is_first_numeric(lines1[i]))
					  continue;

				  tmp_bp = extract_column_integer(lines1[i],1,2);
				  lbp = extract_column_double(lines2[i],1,6);

				  for(j=0;j<max_num_bp;j++)	{
					  if( ((max_bp[j][0]==tmp_bp[0]) && (max_bp[j][1]==tmp_bp[1]) ) || ( (max_bp[j][0]==tmp_bp[1]) && (max_bp[j][1]==tmp_bp[0]) ))	{
						  //printf("%d\t",timer1);
						  for(m=0;m<6;m++)	{
							  bp_param[j].property[m][timer1] = lbp[m];
							 //printf("%f", bp_param[j].property[m][timer1] );
						  }
						  //printf("\n");
					  }
				  }
			}

			break;
		}

	}

	fclose(f_inp_base_pair);
	fclose(f_inp_lbp);
	// Output Files Generation
	FILE *fOut[6], *fLBPM=NULL;
	char **leg,**setname;
	char *fnOut[6], *fnLBPM;
	const char *property[6]=	{ "Shear", "Stretch", "Stagger", "Buckle", "Propeller", "Opening" };
	char title[256];

	snew(leg,max_num_bp);

	for(i=0;i<max_num_bp;i++)	{
		snew(leg[i],256);
		sprintf(leg[i], "%d=%d",max_bp[i][0],max_bp[i][1]);
	}

	for(m=0;m<6;m++)	{
		snew(fnOut[m],256);
		sprintf(fnOut[m],"%s_%s.xvg",property[m],ComName[0]);
		sprintf(title,"Local base pair parameter: %s",property[m]);
		fOut[m] = xvgropen(fnOut[m],title,output_env_get_time_label(oenv),property[m], oenv);
		xvgr_legend(fOut[m],max_num_bp,(const char **) leg,oenv);
	}

	for(n=0;n<nframe;n++)	{

		for(m=0;m<6;m++)
			fprintf(fOut[m],"\n%12.7f",time[n]);

		for(i=0;i<max_num_bp;i++)
			for(m=0;m<6;m++)
				fprintf(fOut[m], "%10.3f",bp_param[i].property[m][n]);
	}

	fprintf(stderr,"\rFinished reading parameters from %s file\n", fn_inp_lbp);
	fprintf(stderr,"Following output files generated:\n");
	for(m=0;m<6;m++)
		fprintf(stderr,"\t%s\n",fnOut[m]);

	for(m=0;m<6;m++)
		xvgrclose(fOut[m]);

	if(bAvg)	{
		fprintf(stderr,"\rCalculating Average and Standard Deviation of local base pair parameters...");

		snew(fnLBPM,256);
		sprintf(fnLBPM,"Avg_Local_BP_param_%s.xvg",ComName[0]);
		fLBPM = xvgropen_type(fnLBPM,"Local base pair parameter","Base Pair","Parameters Value",exvggtXYDY, oenv);
		for(j=0;j<max_num_bp;j++)
			fprintf(fLBPM,"# %d\t%d=%d\n",j+1,max_bp[j][0],max_bp[j][1]);
		for(i=0;i<max_num_bp;i++)	{
			snew(bp_param[i].avg,6);
			snew(bp_param[i].std,6);
			for(j=0;j<6;j++)	{
				calculate_avg_std(bp_param[i].property[j], nframe , dum, &bp_param[i].avg[j], &bp_param[i].std[j]);
			}
		}

		snew(setname,6);
		for(i=0;i<6;i++)
			snew(setname[i],256);
		snew(setname[1],256);
		sprintf(setname[0],"Shear");
		sprintf(setname[1],"Stretch");
		sprintf(setname[2],"Stagger");
		sprintf(setname[3],"Buckle");
		sprintf(setname[4],"Propeller");
		sprintf(setname[5],"Opening");
		xvgr_new_dataset(fLBPM,0,6,(const char **) setname,oenv);

		for(i=0;i<6;i++)	{
			for(j=0;j<max_num_bp;j++)
				fprintf(fLBPM,"%d %10.3f %10.3f\n",(j+1),bp_param[j].avg[i],bp_param[j].std[i]);
			fprintf(fLBPM,"&\n\n");
		}

		fprintf(stderr,"\rCalculated Average and Standard Deviation of local base pair parameters. Output file: %s\n",fnLBPM);
	}

	xvgrclose(fLBPM);
	sfree(bp_param);
	free(buff1);
	free(buff2);
	free(lines1);
	free(lines2);
	free(time_temp);
	free(time);
	return 0;
}


void hbond_process(int nframe, int **max_bp, int max_num_bp, char *fn_inp_base_pair,char *fn_inp_hbond,int NFILE, const char *ComName[],t_filenm fnm[], const gmx_output_env_t *oenv)	{
	FILE *f_inp_base_pair, *f_inp_hbond, *fout_log, *fout_map;
	const char *fn_out_map, *fn_out_log;
	int i,j,n;
	real *time;
	int timer1=0, timer2=0;
	char *buff1, *buff2, **Bp_Map;
	char **lines1, **lines2, **time_temp;
	int num1=0, num2=0, *tmp_bp, *hbond;


	fn_out_map = opt2fn("-map",NFILE,fnm);
	fn_out_log = opt2fn("-g",NFILE,fnm);

	snew(Bp_Map,max_num_bp);
	for(i=0;i<max_num_bp;i++)	{
		snew(Bp_Map[i],nframe+1);
		for(n=0;n<nframe;n++)
			Bp_Map[i][n] = '0';
	}
	snew(time,nframe);

	//Writing base pair log file
	fout_log = gmx_ffopen(fn_out_log,"w");
	fprintf(fout_log,"Index\t\tBase Pair\n");
	for(i=0;i<max_num_bp;i++)	{
		fprintf(fout_log,"%d\t\t%d\t%d\n",i,max_bp[i][0],max_bp[i][1]);
	}
	gmx_ffclose(fout_log);

	//Reading base-pairs and hydrogen bond number from X3DNA output
	fprintf(stderr,"Reading file.....");
	f_inp_base_pair = fopen(fn_inp_base_pair, "r");
	f_inp_hbond = fopen(fn_inp_hbond, "r");

	//Initialization First Frame
	while(1)	{
		buff1 = get_line(f_inp_base_pair);
		remove_leading_white_space(buff1);
		if(strstr(buff1,"# Time")!=NULL)
			break;
	}

	while(1)	{
		buff2 = get_line(f_inp_hbond);
		remove_leading_white_space(buff1);
		if(strstr(buff2,"# Time")!=NULL)
			break;
	}


	time_temp = split_by_char(buff1, "=", NULL);
	time[0] = strtof(time_temp[1],NULL);


	//Extracting data from all frames
	while(1)	{

		lines1 = get_block_lines(f_inp_base_pair,"# Time",&num1);
		lines2 = get_block_lines(f_inp_hbond,"# Time",&num2);


		for(i=0;i<num1;i++)	{
			  if(!is_first_numeric(lines1[i]))
				  continue;

			  tmp_bp = extract_column_integer(lines1[i],1,2);
			  hbond = extract_column_integer(lines2[i],1,1);

			  for(j=0;j<max_num_bp;j++)	{
				  if( ((max_bp[j][0]==tmp_bp[0]) && (max_bp[j][1]==tmp_bp[1]) ) || ( (max_bp[j][0]==tmp_bp[1]) && (max_bp[j][1]==tmp_bp[0]) ))	{
					  Bp_Map[j][timer1] = (char)(((int)'0')+hbond[0]);
					  //printf("%d   %d   %d\n", j, (timer1), hbond[0]);
				  }
			  }
		}

		timer1++;
		timer2++;


		time_temp = split_by_char(lines1[num1], "=", NULL);
		time[timer1] = strtof(time_temp[1],NULL);


		//Processing for last frame
		if(timer1==nframe-1)	{
			lines1 = get_all_lines(f_inp_base_pair,&num1);
			lines2 = get_all_lines(f_inp_hbond,&num2);

			for(i=0;i<num1;i++)	{
				  if(!is_first_numeric(lines1[i]))
					  continue;

				  tmp_bp = extract_column_integer(lines1[i],1,2);
				  hbond = extract_column_integer(lines2[i],1,1);

				  for(j=0;j<max_num_bp;j++)	{
					  if( ((max_bp[j][0]==tmp_bp[0]) && (max_bp[j][1]==tmp_bp[1]) ) || ( (max_bp[j][0]==tmp_bp[1]) && (max_bp[j][1]==tmp_bp[0]) ))	{
						  Bp_Map[j][timer1] = (char)(((int)'0')+hbond[0]);
						  //printf("%d   %d   %d\n", j, (timer1), hbond[0]);
					  }
				  }
			}

			break;
		}

	}

	fclose(f_inp_base_pair);
	fclose(f_inp_hbond);

	//Map matrix generation
	real *bp_index;
	t_matrix matrix;
	t_mapping *map;
	real *r, *g, *b;
	t_xpmelmt c;
	char temp[STRLEN];


	snew(bp_index,max_num_bp);
	for(i=0;i<max_num_bp;i++)
		bp_index[i]=(i+1);
	matrix.axis_x = time;
	matrix.axis_y = bp_index;
	matrix.nx=nframe;
	matrix.ny=max_num_bp;
	sprintf(matrix.title,"Base Pairs with respect to Time");
	matrix.legend[0] =0;
	sprintf(matrix.label_x, "%s",output_env_get_time_label(oenv));
	sprintf(matrix.label_y,"Base Pair Index");
	matrix.bDiscrete=TRUE;

	//Creating color map for XPM file
	snew(map,4);
	snew(r,4);
	snew(g,4);
	snew(b,4);
	r[0] = 1;	g[0] = 1;	b[0] = 1;
	r[1] = 1;	g[1] = 0;	b[1] = 0;
	r[2] = 0;	g[2] = 0;	b[2] = 1;
	r[3] = 0;	g[3] = 1;	b[3] = 0;
	for(i=0;i<4;i++)	{
		map[i].code.c1 = (char)(((int)'0')+i);
		map[i].code.c2=0;
		temp[0] = (char)(((int)'0')+i);		temp[1] = '\0';
		map[i].desc = strdup(strcat(temp," Hydrogen Bond"));
		map[i].rgb.r = r[i];
		map[i].rgb.g = g[i];
		map[i].rgb.b = b[i];
	}
	matrix.map = map;
	matrix.nmap = 4;

	//Creating main matrix for XPM
	snew(matrix.matrix,nframe);
	for(i=0;i<nframe;i++)
		snew(matrix.matrix[i],max_num_bp);
	c.c2=0;
	for(i=0;i<max_num_bp;i++){
		for(j=0;j<nframe;j++){
			c.c1=Bp_Map[i][j];
			matrix.matrix[j][i] = std::max( (short)0, searchcmap(matrix.nmap,matrix.map,c));
			//printf("%d",matrix.matrix[j][i]);
		}
	}

	fn_out_map = opt2fn("-map",NFILE,fnm);
	matrix.flags=0;
	fout_map =gmx_ffopen(fn_out_map ,"w");
	write_xpm_m(fout_map ,matrix);
	gmx_ffclose(fout_map );

	fprintf(stderr,"\rFinished reading file and Hydrogen bond map generated.\n");

	sfree(matrix.matrix);
	free(buff1);
	free(buff2);
	free(Bp_Map);
	free(lines1);
	free(lines2);
	free(time_temp);
	free(time);
}

int** get_max_base_pairs(char *filename, int *num_bp)	{
	int **bp=NULL, *tmp_bp=NULL, max_num_bp=0;
	char *buffer=NULL;
	gmx_bool bFirst=TRUE, bToAdd=FALSE;
	int i=0, timer=0;
	FILE *fInput;

	fInput = fopen(filename,"r");

	snew(bp,1);
	fprintf(stderr,"No.\tResidue Numbers\n");

	while(!feof(fInput))	{

		  buffer = get_line(fInput);
		  remove_leading_white_space(buffer);

		  if(buffer==NULL)
			  continue;

		  if(strstr(buffer,"# Time")!=NULL)	{
			  timer++;
			  continue;
		  }

		  if(timer==2)
			  bFirst=FALSE;

		  if((bFirst) && (timer==1))	{

			  if(!is_first_numeric(buffer))
				  continue;

			  max_num_bp += 1;
			  srenew(bp, max_num_bp);
			  bp[max_num_bp-1] = extract_column_integer(buffer,1,2);
			  fprintf(stderr,"%d:\t%d = %d\n",max_num_bp, bp[max_num_bp-1][0],bp[max_num_bp-1][1]);
		  }

		  if(timer > 1)		{

			  if(!is_first_numeric(buffer))
				  continue;
			  bToAdd = TRUE;
			  tmp_bp = extract_column_integer(buffer,1,2);
			  for(i=0;i<max_num_bp;i++)	{
				  if(     ((tmp_bp[0]== bp[i][0]) && (tmp_bp[1]==bp[i][1])) \
					  ||  ((tmp_bp[0]== bp[i][1]) && (tmp_bp[1]==bp[i][0])) )		{
					  bToAdd = FALSE;
					  break;
				  }

			  }
			  if(bToAdd)	{
				  max_num_bp += 1;
				  srenew(bp,max_num_bp);
				  bp[max_num_bp-1] = extract_column_integer(buffer,1,2);
				  fprintf(stderr,"%d:\t%d = %d\n",max_num_bp, bp[max_num_bp-1][0],bp[max_num_bp-1][1]);
			  }

		  }
	}

	fclose(fInput);
	free(buffer);
	//free(tmp_bp);

	*num_bp = max_num_bp;
	return bp;
}

int** get_max_base_pairs_step(char *filename, int *num_bp)	{
	int **bps=NULL, *tmp_bps_prev=NULL, *tmp_bps_curr=NULL, max_num_bps=0;
	char *buffer=NULL;
	gmx_bool bFirst=TRUE, bToAdd=FALSE, bFirstBp=TRUE;
	int i=0, timer=0;
	FILE *fInput;

	fInput = fopen(filename,"r");

	snew(bps,1);
	fprintf(stderr,"No.\tResidue Numbers\n");

	while(!feof(fInput))	{

		  buffer = get_line(fInput);
		  remove_leading_white_space(buffer);

		  if(buffer==NULL)
			  continue;

		  if(strstr(buffer,"# Time")!=NULL)	{
			  timer++;
			  bFirstBp = TRUE;
			  continue;
		  }

		  if(timer==2)
			  bFirst=FALSE;

		  if((bFirst) && (timer==1))	{

			  if(!is_first_numeric(buffer))
				  continue;

			  if(bFirstBp)	{
				  tmp_bps_prev = extract_column_integer(buffer,1,2);
				  bFirstBp = FALSE;
				  continue;
			  }

			  max_num_bps += 1;
			  srenew(bps, max_num_bps);
			  snew(bps[max_num_bps-1],4);

			  tmp_bps_curr = extract_column_integer(buffer,1,2);
			  bps[max_num_bps-1][0] = tmp_bps_prev[0];
			  bps[max_num_bps-1][1] = tmp_bps_prev[1];
			  bps[max_num_bps-1][2] = tmp_bps_curr[0];
			  bps[max_num_bps-1][3] = tmp_bps_curr[1];
			  printf("%d:\t%d = %d || %d = %d\n",max_num_bps,bps[max_num_bps-1][0],bps[max_num_bps-1][1],bps[max_num_bps-1][2],bps[max_num_bps-1][3]);
		  }

		  if(timer > 1)		{

			  if(!is_first_numeric(buffer))
				  continue;
			  bToAdd = TRUE;

			  if(bFirstBp)	{
				  tmp_bps_prev = extract_column_integer(buffer,1,2);
				  bFirstBp = FALSE;
				  continue;
			  }

			  tmp_bps_curr = extract_column_integer(buffer,1,2);

			  for(i=0;i<max_num_bps;i++)	{
				  if	(
							(
								((tmp_bps_prev[0] == bps[i][0] && tmp_bps_prev[1] == bps[i][1]) || (tmp_bps_prev[1] == bps[i][0] && tmp_bps_prev[0] == bps[i][1]))
							&&
								((tmp_bps_curr[0] == bps[i][2] && tmp_bps_curr[1] == bps[i][3]) || (tmp_bps_curr[1] == bps[i][2] && tmp_bps_curr[0] == bps[i][3]))
							)
						||
							(
								((tmp_bps_curr[0] == bps[i][0] && tmp_bps_curr[1] == bps[i][1]) || (tmp_bps_curr[1] == bps[i][0] && tmp_bps_curr[0] == bps[i][1]))
							&&
								((tmp_bps_prev[0] == bps[i][2] && tmp_bps_prev[1] == bps[i][3]) || (tmp_bps_prev[1] == bps[i][2] && tmp_bps_prev[0] == bps[i][3]))
							 )
					  	 )
				  {
					  bToAdd = FALSE;
					  break;
				  }

			  }
			  if(bToAdd)	{
				  max_num_bps += 1;
				  srenew(bps,max_num_bps);
				  snew(bps[max_num_bps-1],4);
				  bps[max_num_bps-1][0] = tmp_bps_prev[0];
				  bps[max_num_bps-1][1] = tmp_bps_prev[1];
				  bps[max_num_bps-1][2] = tmp_bps_curr[0];
				  bps[max_num_bps-1][3] = tmp_bps_curr[1];
				  printf("%d:\t%d = %d || %d = %d\n",max_num_bps,bps[max_num_bps-1][0],bps[max_num_bps-1][1],bps[max_num_bps-1][2],bps[max_num_bps-1][3]);

			  }

		  }

		  if(!bFirstBp)
			  tmp_bps_prev=tmp_bps_curr;

	}

	fclose(fInput);
	free(buffer);

	*num_bp = max_num_bps;
	return bps;
}

void write_time(real t, FILE *f_cum_data[])	{
	int i=0;
	for(i=0;i<10;i++)	{
		fprintf(f_cum_data[i],"\n# Time = %15.5f\n",t);
	}
}

int add_data_to_files(char *fn_input, FILE *f_cum_data[])		{
	int i = 0, j= 0, bp = 0;
	gmx_bool bProp[10], bS1=FALSE, bS2=FALSE;
	FILE *f_input;
	int number=0, elem=0;
	char **data=NULL, **two_strand=NULL, *tmp_str, *bp1, *bp2;
	char **SplitData = NULL;
	int num_bp=0;

	f_input = fopen(fn_input, "r");
	data = get_all_lines(f_input, &number);
	fclose(f_input);

	for(i=0;i<10;i++)
		bProp[i]=FALSE;

	for(i=0;i<number;i++)	{
		if (data[i]==NULL)
			continue;

		if(strstr(data[i],"Number of base-pairs")!=NULL)	{
			SplitData = split_by_char(data[i], ":", NULL);
			num_bp = atoi(SplitData[1]);
			free(SplitData);
		}

		if(strstr(data[i],"****************")!=NULL)	{
			bProp[eBasePairs] = FALSE;
			bProp[eHbond] = FALSE;
			bProp[eMgroove] = FALSE;
			bProp[eHelixRad] = FALSE;

			if (bProp[eBBnDihedral])	{
				for(bp=0;bp<num_bp;bp++)
					fprintf(f_cum_data[eBBnDihedral],"%s\n", two_strand[bp]);
				bS2 = FALSE;
				bProp[eBBnDihedral] = FALSE;
				bp = 0;
				for(j=0; j<num_bp; j++)
						free(two_strand[j]);
				free(two_strand);
			}

			if (bProp[eSugarConf])	{
				for(bp=0;bp<num_bp;bp++)
					fprintf(f_cum_data[eSugarConf],"%s\n", two_strand[bp]);
				bS2 = FALSE;
				bProp[eSugarConf] = FALSE;
				bp = 0;
				for(j=0; j<num_bp; j++)
						free(two_strand[j]);
				free(two_strand);
			}


			continue;
		}

		if(strstr(data[i],"~~~~~~~~~~~~~~~~")!=NULL)	{
			bProp[eLBP] = FALSE;
			bProp[eLBPS] = FALSE;
			bProp[eLBPH] = FALSE;
			continue;
		}

		if(strstr(data[i],"RMSD of the bases")!=NULL)	{
			bProp[eBasePairs] = TRUE;
			continue;
		}
		if(strstr(data[i],"Detailed H-bond information")!=NULL)		{
			bProp[eHbond] = TRUE;
			continue;
		}
		if(strstr(data[i],"Local base-pair parameters")!=NULL)		{
			bProp[eLBP] = TRUE;
			continue;
		}
		if(strstr(data[i],"Local base-pair step parameters")!=NULL)		{
			bProp[eLBPS] = TRUE;
			continue;
		}
		if(strstr(data[i],"Local base-pair helical parameters")!=NULL)		{
			bProp[eLBPH] = TRUE;
			continue;
		}
		if(strstr(data[i],"Position (Px, Py, Pz) and local helical axis")!=NULL)		{
			bProp[eHelAxis] = TRUE;
			continue;
		}
		if(strstr(data[i],"Minor and major groove widths")!=NULL)		{
			bProp[eMgroove] = TRUE;
			continue;
		}

		if(strstr(data[i],"Helix radius")!=NULL)		{
			bProp[eHelixRad] = TRUE;
			continue;
		}


		if(strstr(data[i],"Main chain and chi torsion angles")!=NULL)		{
			bProp[eBBnDihedral] = TRUE;
			two_strand = (char **) malloc (num_bp * sizeof(char*));
			for(j=0; j<num_bp; j++)
				two_strand[j] = (char *) malloc (160 * sizeof(char));
			continue;
		}

		if(strstr(data[i],"Sugar conformational parameters")!=NULL)		{
			bProp[eSugarConf] = TRUE;
			two_strand = (char **) malloc (num_bp * sizeof(char*));
			for(j=0; j<num_bp; j++)
				two_strand[j] = (char *) malloc (200 * sizeof(char));
			continue;

		}

		if (bProp[eBBnDihedral])	{
			if ((strstr(data[i], "Strand I")!=NULL) && (bS1==FALSE))	{
				bS1 = TRUE;
				bS2 = FALSE;
				bp = 0;
				continue;
			}
			if((strstr(data[i], "Strand II")!=NULL) && (bS2==FALSE))	{
				bS1 = FALSE;
				bS2 = TRUE;
				bp = 0;
				continue;
			}
		}

		if (bProp[eSugarConf])	{
			if ((strstr(data[i], "Strand I")!=NULL) && (bS1==FALSE))	{
				bS1 = TRUE;
				bS2 = FALSE;
				bp = 0;
				continue;
			}
			if((strstr(data[i], "Strand II")!=NULL) && (bS2==FALSE))	{
				bS1 = FALSE;
				bS2 = TRUE;
				bp = 0;
				continue;
			}
		}


		if(!is_first_numeric(data[i]))
				continue;

		if(bProp[eBasePairs])	{
			SplitData = split_by_char(data[i], ":", &elem);
			bp1 = extract_digits(SplitData[1]);
			bp2 = extract_digits(SplitData[3]);
			fprintf(f_cum_data[eBasePairs],"%s   %s\n",bp1,bp2);
			//for (j=0; j<elem; j++)
				//free(SplitData[j]);
			free(SplitData);
			free(bp1);
			free(bp2);
		}

		if(bProp[eHbond])	{
			SplitData = split_by_space(data[i], NULL);
			bp1 = extract_digits(SplitData[2]);
			fprintf(f_cum_data[eHbond],"%s\n", bp1);
			free(SplitData);
			free(bp1);
		}

		if(bProp[eLBP])	{
			SplitData = split_by_space(data[i], NULL);
			fprintf(f_cum_data[eLBP],"%s   %s   %s   %s   %s   %s\n",SplitData[2],SplitData[3],SplitData[4],SplitData[5],SplitData[6],SplitData[7]);
			free(SplitData);
		}

		if(bProp[eLBPS])	{
			SplitData = split_by_space(data[i], NULL);
			fprintf(f_cum_data[eLBPS],"%s   %s   %s   %s   %s   %s\n",SplitData[2],SplitData[3],SplitData[4],SplitData[5],SplitData[6],SplitData[7]);
			free(SplitData);
		}

		if(bProp[eLBPH])	{
			SplitData = split_by_space(data[i], NULL);
			fprintf(f_cum_data[eLBPH],"%s   %s   %s   %s   %s   %s\n",SplitData[2],SplitData[3],SplitData[4],SplitData[5],SplitData[6],SplitData[7]);
			free(SplitData);
		}

		if(bProp[eMgroove])	{
			SplitData = split_by_space(data[i], NULL);
			fprintf(f_cum_data[eMgroove],"%s   %s   %s   %s\n",SplitData[2],SplitData[3],SplitData[4],SplitData[5]);
			free(SplitData);
		}

		if(bProp[eHelAxis])	{
			SplitData = split_by_space(data[i], NULL);
			fprintf(f_cum_data[eHelAxis],"%s   %s   %s   %s   %s   %s\n",SplitData[2],SplitData[3],SplitData[4],SplitData[5],SplitData[6],SplitData[7]);
			free(SplitData);
		}

		if(bProp[eHelixRad])	{
			SplitData = split_by_space(data[i], NULL);
			fprintf(f_cum_data[eHelixRad],"%s   %s   %s   %s   %s   %s\n",SplitData[2],SplitData[3],SplitData[4],SplitData[5],SplitData[6],SplitData[7]);
			free(SplitData);
		}

		if(bProp[eBBnDihedral])	{
			SplitData = split_by_space(data[i], NULL);
			tmp_str = (char *) malloc (80 * sizeof(char));
			if (bS1)	{
				sprintf(tmp_str, "%s   %s   %s   %s   %s   %s  %s  ",SplitData[2],SplitData[3],SplitData[4],SplitData[5],SplitData[6],SplitData[7],SplitData[8]);
				strcpy(two_strand[bp], tmp_str);
				bp += 1;
			}

			if (bS2)	{
				sprintf(tmp_str, "%s   %s   %s   %s   %s   %s  %s",SplitData[2],SplitData[3],SplitData[4],SplitData[5],SplitData[6],SplitData[7],SplitData[8]);
				strcat(two_strand[bp], tmp_str);
				bp += 1;
			}
			free(SplitData);
			free(tmp_str);
		}

		if(bProp[eSugarConf])	{
			SplitData = split_by_space(data[i], NULL);
			tmp_str = (char *) malloc (100 * sizeof(char));
			if (bS1)	{
				sprintf(tmp_str, "%s   %s   %s   %s   %s   %s  %s  %s  ",SplitData[2],SplitData[3],SplitData[4],SplitData[5],SplitData[6],SplitData[7],SplitData[8],SplitData[9]);
				strcpy(two_strand[bp], tmp_str);
				bp += 1;
			}

			if (bS2)	{
				sprintf(tmp_str, "%s   %s   %s   %s   %s   %s  %s  %s",SplitData[2],SplitData[3],SplitData[4],SplitData[5],SplitData[6],SplitData[7],SplitData[8],SplitData[9]);
				strcat(two_strand[bp], tmp_str);
				bp += 1;
			}
			free(SplitData);
			free(tmp_str);
		}

	}

	for(i=0;i<number+1;i++)
		free(data[i]);
	free(data);
	return num_bp;
}

int gmx_3dna(int argc,char *argv[])
{

  const char *desc[] = {
		  "do_x3dna uses 3DNA package to calculate several structural descriptors or parameters ",
		  "of DNA/RNA using the GROMACS MD trajectory. It extracts output of the 3DNA package,"
		  "and saves these parameters to external output files with function of time.\n\n",

		  "To execute do_x3dna, 3DNA package should be installed and $X3DNA environment ",
		  "variable (Detail is given in 3DNA manual) should be defined.\n\n",

		  "NOTE 1: Before running do_x3dna, make sure 3DNA package is correctly working for ",
		  "the input DNA/RNA structure. To check it, run do_x3dna with a input gro/pdb file containing",
		  "only one frame instead of xtc/trr trajectory file as follows : \n",
		  "do_x3dna -s topol.tpr -f check.gro -n index.ndx -noisy\n",
		  "\"-noisy\" option switch on the output message from the 3DNA package on the display, "
		  "that would be necessary to troubleshoot problems related to the 3DNA package. ",
		  "Most common problem is the residue names mismatch in input DNA/RNA structure "
		  " and the 3DNA package dictionary.\n\n"

		  "NOTE 2: Only PBC corrected trajectory and tpr files should be used as inputs. ",
		  "PBC corrected PDB/GRO file can be used in place of tpr file. Index file SHOULD BE",
		  "used to select ONLY DNA/RNA in the input trajectory file.\n\n",

		  "NOTE 3: If \"-ref\" option is used, all base-pairs/steps parameters will be calculated on",
		  "the basis of the structure of the DNA/RNA present in the input tpr/pdb file with \"-s\" option."
		  " \"-ref\" option SHOULD BE USED if output files are further used as input files in Python APIs or",
		  "Scripts that are provided with do_x3dna package. To analyze the formation or breaking of ",
		  "base-pairs during MD simulations, either \"-noref\" could be used or do not include \"-ref\" ",
		  "option because this option is switched off by default .\n\n"

		  "NOTE 4: If \"-fit\" is enabled, during fitting procedure the DNA/RNA is translated to origin ",
		  "such that its center of mass is located at the origin. Most of the parameters are unaffected ",
		  "by this fitting, however coordinates of the local helical axis could mismatch with the input ",
		  "coordinates with the DNA/RNA.\n\n"

		  "\"-hbond\" option extracts hydrogen bonds for each base pair. ",
		  "A map.log (\"-g\") file is generated containing the base pair information as per index ",
		  "of the the hydrogen bond map (\"-map\").\n\n",

		  "\"-lbpm\" option calculates Local Base Pair Parameters (Shear, Stretch, Stagger, Buckle, ",
		  "Propeller and Opening) with function of time, and average (with \"-avg\") of these parameters ",
		  "with function of the base-pairs.\n\n"

		  "\"-lbpsm\" option calculates Local Base Pair-Step Parameters (Shift, Slide, Rise, Tilt, Roll ",
		  "and Twist) with function of time, and average (with \"-avg\") of these parameters with function",
		  "of the base-steps.\n\n",

		  "\"-lbphm\" option calculates Local Base Pair-Helical Parameters (X-displacement, Y-displacement,",
		  "H-rise, Inclination, Tip and H-twist) with function of time, and average (with \"-avg\") of these",
		  "parameters with function of the base-steps.\n\n",

		  "Also, all the above parameters including local helical axis, major and minor grooves, local helical",
		  "radius, backbone dihedral angles (alpha, beta, gamma, delta, epsilon, zeta and chi) of both strands",
		  "are calculated using 3DNA package for each frame and written in separate files with function of time.\n"
		  "OUTPUT FILE LIST:\n",
		  "    base_pairs_g.dat             => Base-pairs\n",
		  "    h-bond_g.dat                 => Hydrogen bonds between base-pairs\n",
		  "    L-BP_g.dat                   => Base-pairs parameters\n",
		  "    L-BPS_g.dat                  => Base-steps parameters\n",
		  "    L-BPH_g.dat                  => Helical Base-steps parameters\n",
		  "    HelAxis_g.dat                => Local helical axis coordinates\n",
		  "    MGroove_g.dat                => Major and Minor grooves\n",
		  "    HelixRad_g.dat               => Local helical radius\n",
		  "    BackBoneCHiDihedrals_g.dat   => Backbone dihederal angles including Chi-dihedral\n",
		  "    SugarDihedrals_g.dat         => Sugar dihederal angles including puckring type\n"
		  "Name of these files could be change by setting different suffix instead of \"g\" using \"-name\" option.",
		  "These files could be used with the Python APIs or scripts for further analysis.\n"
  };


  gmx_bool bLBPM=FALSE, bRef=FALSE, bVerbose=FALSE, bHbond=FALSE, bLBPSM=FALSE, bLBPHM=FALSE, bAvg=TRUE;
  gmx_bool bFit=TRUE, bM=TRUE, bAnalyzeC_Option=FALSE;
  static const char *ComName[]  = { "g" };
  gmx_output_env_t *oenv;

  t_pargs pa[] = {
		  { "-noisy", FALSE, etBOOL, {&bVerbose}, "Generate information from the X3DNA package" },
		  { "-hbond", FALSE, etBOOL, {&bHbond}, "Hydrogen bond map for base pairs"},
		  { "-ref", FALSE, etBOOL, {&bRef}, "Base pair parameters will be calculated from base pair of the reference frame"},
		  { "-name", FALSE, etSTR, {&ComName}, "Output file names will be suffixed by this word after \"_\"" },
		  { "-fit", TRUE, etBOOL, {&bFit}, "Fitting frames on reference structure" },
		  { "-mwa", TRUE, etBOOL, {&bM},  "Mass weighted fitting" },
		  { "-lbpm", FALSE, etBOOL, {&bLBPM}, "To calculate local base pair parameters" },
		  { "-lbpsm", FALSE, etBOOL, {&bLBPSM}, "To calculate local base-pair step parameters" },
		  { "-lbphm", FALSE, etBOOL, {&bLBPHM}, "To calculate local base-pair helical parameters" },
		  { "-avg", TRUE, etBOOL, {&bAvg}, "Average and Standard Deviation over all the frames" },
      { "-c", TRUE, etBOOL, {&bAnalyzeC_Option}, "Output structural parameters between helical regions (\"----\" by default). It will invoke \"-c\" option with 3DNA analyze command." }
    };

  t_filenm   fnm[] = {
     { efTRX, "-f",   NULL,      ffREAD },
     { efTPS, NULL,   NULL,      ffREAD },
     { efNDX, NULL,   NULL,      ffOPTRD },
     { efXVG, NULL,  "BP_count", ffWRITE },
     { efXPM, "-map", "BP_map",   ffOPTWR},
     { efLOG, "-g",	  "map",     ffOPTWR},
   };

  #define NFILE asize(fnm)
  int npargs;
  CopyRightMsg();
  npargs = asize(pa);

  if ( ! parse_common_args(&argc,argv, PCA_CAN_TIME | PCA_TIME_UNIT | PCA_CAN_VIEW , NFILE,fnm,npargs,pa, asize(desc),desc,0,NULL,&oenv) )
  {
   	return 0;
  }

  //GROMACS stuffs
  t_trxstatus *status;
  t_topology top;
  int        ePBC;
  t_atoms    *atoms;
  int    *index, *ifit;
  matrix     box;
  gmx_rmpbc_t  gpbc=NULL;
  int        nres,nr0;
  real       t, *w_rls=NULL;
  int        gnx, nfit;
  int 		*nres_index;
  char       *grpnm, *fitname;
  rvec       *xref,*x, x_shift;
  int        i,natoms, nframe=0, num_bp=0;


  FILE       *tapein;
  FILE       *fnum_bp,*tmpf;
  const char *fnBP_Count;
  char       suffix_name[32], pdbfile[32],inpfile[32],x3dna_out_file[32], title[256];
  char       find_pair_cmd[256], analyze_cmd[256];
  const char *dptr, *leg[1] = {"No. of BP"};
  char		  fn_cum_data[10][32];
  FILE		  *f_cum_data[10];

  fnBP_Count= opt2fn("-o",NFILE,fnm);

  read_tps_conf(ftp2fn(efTPS,NFILE,fnm), &top, &ePBC, &xref, NULL, box, FALSE);
  atoms=&(top.atoms);

  if (bFit)		{
	  if (fn2bTPX(ftp2fn(efTPS,NFILE,fnm)) == FALSE)	{
		  printf("\n\nWARNING: Not a TPR file.... Switching off mass weighted fitting...\n\n");
		  bM = FALSE;
	  }

	  printf("\nChoose a group for the least squares fit\n");
	  get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&nfit,&ifit,&fitname);
	  if (nfit < 3)
		  gmx_fatal(FARGS,"Need >= 3 points to fit!\n");

  }		else
	  nfit=0;

  printf("\nChoose a group for the X3DNA analysis\n");
  get_index(atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&gnx,&index,&grpnm);

   //Initialization of fitting stuff
  if(bFit)		{
	  snew(w_rls,atoms->nr);
    for(i=0; (i<nfit); i++) {
      if(bM)  {
        w_rls[ifit[i]]=atoms->atom[ifit[i]].m;
      }
	    else {
        w_rls[ifit[i]]=1.0;
      }
    }
  }

  //To get number and array of residues in selected group......
  nres=0;
  nr0=-1;
  snew(nres_index,1);
  for(i=0; (i<gnx); i++) 	{
	  if (atoms->atom[index[i]].resind != nr0) 		{
		  nr0=atoms->atom[index[i]].resind;
		  srenew(nres_index,nres+1);
		  nres_index[nres]=atoms->resinfo[nr0].nr;
		  //printf("%d %d\n",nres,nres_index[nres]);
		  nres++;
	  }
  }
  fprintf(stderr,"There are %d residues in your selected group\n",nres);


  //Creating name for temporary PDB file
  strcpy(suffix_name,"ddXXXXXX");
  sprintf(suffix_name,"ddXXXXXX");
  gmx_tmpnam(suffix_name);
  remove(suffix_name);
  sprintf(pdbfile,"%s.pdb",suffix_name);

  if ((tmpf = fopen(pdbfile,"w")) == NULL)
	  gmx_fatal(FARGS,"Can not open pdb file %s",pdbfile);
  fclose(tmpf);
  remove(pdbfile);

  //Creating name for temporary input file for find_pair program
  sprintf(inpfile,"%s.inp",suffix_name);
  if ((tmpf = fopen(inpfile,"w")) == NULL)
	  gmx_fatal(FARGS,"Can not open inp file %s",inpfile);
  fclose(tmpf);
  remove(inpfile);

  //Creating name for output file of analyze program
  sprintf(x3dna_out_file,"%s.out",suffix_name);

  //Creating variable for executing command of find_pair from X3DNA
  dptr=getenv("X3DNA");
  if (!gmx_fexist(dptr))
	  gmx_fatal(FARGS,"$X3DNA environment not found", dptr);
  sprintf(find_pair_cmd,"$X3DNA/bin/find_pair %s %s %s",pdbfile, inpfile, bVerbose?"":"2> /dev/null");
  if (bVerbose)
	  fprintf(stderr,"find_pair command='%s'\n",find_pair_cmd);

  //Creating variable for executing command of analyze from X3DNA
  if (bAnalyzeC_Option)
    sprintf(analyze_cmd,"$X3DNA/bin/analyze -c %s %s", inpfile,bVerbose?"":"2> /dev/null");
  else
    sprintf(analyze_cmd,"$X3DNA/bin/analyze %s %s", inpfile,bVerbose?"":"2> /dev/null");
  if (bVerbose)
	  fprintf(stderr,"find_pair command='%s'\n",analyze_cmd);

  //Reading first frame
  natoms=read_first_x(oenv,&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);

  // Fitting first frame to initial structure
  if (bFit)	{
    // Translate the reference structure to origin and store the original position
    copy_rvec(xref[index[0]], x_shift);
    reset_x(nfit,ifit,top.atoms.nr,NULL,xref,w_rls);
    rvec_dec(x_shift, xref[index[0]]);

    // Translate the frame to origin, fit (by rotating) it to reference structure
    // and translate again the frame to original position of reference structure
    reset_x(nfit,ifit,top.atoms.nr,NULL,x,w_rls);
    do_fit(natoms,w_rls,xref,x);
    for (i = 0; (i < natoms); i++)
         rvec_inc(x[i], x_shift);
  }

  if (natoms > atoms->nr)
		  gmx_fatal(FARGS,"\nTrajectory does not match topology!");
  if (gnx > natoms)
	  gmx_fatal(FARGS,"\nTrajectory does not match selected group!");


  //OUTPUT FILES HANDLING
  fnum_bp = xvgropen(fnBP_Count,"Number of base Pairs",output_env_get_time_label(oenv),"No. of BP",oenv);
  xvgr_legend(fnum_bp,1,leg,oenv);

  //Cumulative Data Files Opening
  	sprintf(fn_cum_data[eBasePairs],"base_pairs_%s.dat",ComName[0]);
  	sprintf(fn_cum_data[eHbond],"h-bond_%s.dat",ComName[0]);
  	sprintf(fn_cum_data[eLBP],"L-BP_%s.dat",ComName[0]);
  	sprintf(fn_cum_data[eLBPS],"L-BPS_%s.dat",ComName[0]);
  	sprintf(fn_cum_data[eLBPH],"L-BPH_%s.dat",ComName[0]);
  	sprintf(fn_cum_data[eHelAxis],"HelAxis_%s.dat",ComName[0]);
  	sprintf(fn_cum_data[eMgroove],"MGroove_%s.dat",ComName[0]);
  	sprintf(fn_cum_data[eHelixRad],"HelixRad_%s.dat",ComName[0]);
  	sprintf(fn_cum_data[eBBnDihedral],"BackBoneCHiDihedrals_%s.dat",ComName[0]);
  	sprintf(fn_cum_data[eSugarConf],"SugarDihedrals_%s.dat",ComName[0]);

  	f_cum_data[eBasePairs] = gmx_ffopen(fn_cum_data[eBasePairs],"w");

    f_cum_data[eHbond] = gmx_ffopen(fn_cum_data[eHbond],"w");

    f_cum_data[eLBP] = gmx_ffopen(fn_cum_data[eLBP],"w");
    fprintf(f_cum_data[eLBP],"#Shear    Stretch   Stagger    Buckle  Propeller  Opening\n");

    f_cum_data[eLBPS] = gmx_ffopen(fn_cum_data[eLBPS],"w");
    fprintf(f_cum_data[eLBPS],"#Shift     Slide      Rise      Tilt      Roll     Twist\n");

    f_cum_data[eLBPH] = gmx_ffopen(fn_cum_data[eLBPH],"w");
    fprintf(f_cum_data[eLBPH],"#X-disp    Y-disp   h-Rise     Incl.       Tip   h-Twist\n");

    f_cum_data[eHelAxis] = gmx_ffopen(fn_cum_data[eHelAxis],"w");
    fprintf(f_cum_data[eHelAxis],"#Position (Px, Py, Pz) and local helical axis vector (Hx, Hy, Hz)for each dinucleotide step\n");

    f_cum_data[eMgroove] = gmx_ffopen(fn_cum_data[eMgroove],"w");
    fprintf(f_cum_data[eMgroove],"#Minor Groove        Major Groove\n");
    fprintf(f_cum_data[eMgroove],"#P-P     Refined     P-P     Refined\n");

    f_cum_data[eHelixRad] = gmx_ffopen(fn_cum_data[eHelixRad],"w");
    fprintf(f_cum_data[eHelixRad],"#   Strand I Atoms                Strand II Atoms \n");
    fprintf(f_cum_data[eHelixRad],"#P        O4'       C1'        P        O4'        C1'\n");

    f_cum_data[eBBnDihedral] = gmx_ffopen(fn_cum_data[eBBnDihedral],"w");
    fprintf(f_cum_data[eBBnDihedral],"#Strand I                                                    Strand II \n");
    fprintf(f_cum_data[eBBnDihedral],"#alpha    beta   gamma   delta  epsilon   zeta    chi   |||  alpha    beta   gamma   delta  epsilon   zeta    chi\n");

    f_cum_data[eSugarConf] = gmx_ffopen(fn_cum_data[eSugarConf],"w");
    fprintf(f_cum_data[eSugarConf],"#Strand I                                                              Strand II \n");
    fprintf(f_cum_data[eSugarConf],"#v0      v1      v2      v3      v4      tm       P    Puckering  |||  v0      v1      v2      v3      v4      tm       P    Puckering\n");


    for (i=0;i<10;i++)
    	setbuf(f_cum_data[i],NULL);

  //=======================================================================================
  gpbc = gmx_rmpbc_init(&top.idef,ePBC,natoms);
  if(bRef)		{
	  gmx_rmpbc(gpbc,natoms,box,x);
	  tapein=gmx_ffopen(pdbfile,"w");
	  write_pdbfile_indexed(tapein,NULL,&top.atoms,xref,-1,box,' ',-1,gnx,index,NULL,TRUE);
	  gmx_ffclose(tapein);

	  if(0 != system(find_pair_cmd))
		  gmx_fatal(FARGS,"Failed to execute command: %s",find_pair_cmd);
	  remove(pdbfile);
   }

  //***************************************************************************************
  //START of Trajectory loop
  //***************************************************************************************
  do 	{
	  gmx_rmpbc(gpbc,natoms,box,x);

	  //Fitting the frame to reference structure
	  if (bFit)	{
      // Translate the frame to origin, fit (by rotating) it to reference structure
      // and translate again the frame to original position of reference structure
      reset_x(nfit,ifit,top.atoms.nr,NULL,x,w_rls);
      do_fit(natoms,w_rls,xref,x);
      for (i = 0; (i < natoms); i++)
           rvec_inc(x[i], x_shift);
	  }

	  tapein=gmx_ffopen(pdbfile,"w");
	  write_pdbfile_indexed(tapein,NULL,atoms,x,ePBC,box,' ',-1,gnx,index,NULL,TRUE);
	  gmx_ffclose(tapein);

	  //Executing program $X3DNA/bin/find_pair
	  if(!bRef)
		  if(0 != system(find_pair_cmd))
			  gmx_fatal(FARGS,"Failed to execute command: %s",find_pair_cmd);

	  //Executing program $X3DNA/bin/analyze
	  if(0 != system(analyze_cmd))
		  gmx_fatal(FARGS,"Failed to execute command: %s",analyze_cmd);

	  write_time(t, f_cum_data);
	  num_bp = add_data_to_files(x3dna_out_file, f_cum_data);
	  fprintf(fnum_bp,"%12.7f   %d \n",t, num_bp);

	  remove(pdbfile);
	  nframe++;
  }		while(read_next_x(oenv,status,&t,x,box));

  remove(inpfile);
  remove("bestpairs.pdb"); remove("hel_regions.pdb"); remove("col_helices.scr"); remove("col_chains.scr"); remove("bp_step.par");
  remove("bp_helical.par"); remove("stacking.pdb"); remove("hstacking.pdb"); remove("cf_7methods.par"); remove("auxiliary.par");
  remove("bp_order.dat"); remove("ref_frames.dat");

  for (i=0;i<7;i++);
  	  fclose(f_cum_data[i]);

  //============================
  //Post Processing if enabled
  //============================

  fprintf(stderr,"\nFinished Trajectory Reading and X3DNA Executions....\n");


  if(bHbond || bLBPM || bLBPSM || bLBPHM)	{

	  fprintf(stderr,"========================================================\n");
	  fprintf(stderr,"                        Post Processing                 \n");
	  fprintf(stderr,"========================================================\n");

	  fprintf(stderr,"\rGetting all base pairs....\n");
	  int **bp=NULL, max_num_bp=0;
	  bp = get_max_base_pairs(fn_cum_data[eBasePairs],&max_num_bp);
	  fprintf(stderr,"\nMaximum number of base pairs: %d\n\n", max_num_bp);

	  if(bHbond)	{
		  fprintf(stderr,"=====================Hydrogen Bonds=====================\n");
		  hbond_process(nframe,bp,max_num_bp,fn_cum_data[eBasePairs],fn_cum_data[eHbond],NFILE, ComName, fnm, oenv);
		  fprintf(stderr,"==========================END===========================\n");
	  }

	  if(bLBPM)	{
		  fprintf(stderr,"\n\n========================================================\n");
		  fprintf(stderr,"             Local Base-Pair Parameters                 \n");
		  fprintf(stderr,"========================================================\n");
		  local_base_pair_out(bAvg, ComName,nframe, bp, max_num_bp, fn_cum_data[eBasePairs], fn_cum_data[eLBP],oenv);
		  fprintf(stderr,"=======================END==============================\n");
	  }

	  if(bLBPSM || bLBPHM)	{
		  int **bps=NULL, max_num_bps=0;

		  fprintf(stderr,"\n\nGetting all base pair steps....\n");
		  bps = get_max_base_pairs_step(fn_cum_data[eBasePairs], &max_num_bps);
		  fprintf(stderr,"\nMaximum number of base pair steps: %d\n\n", max_num_bps);

		  if(bLBPSM)	{
			  const char *property[6]=	{ "Shift", "Slide", "Rise", "Tilt", "Roll", "Twist" };
			  fprintf(stderr,"\n\n========================================================\n");
			  fprintf(stderr,"             Local Base-Pair Step Parameters                 \n");
			  fprintf(stderr,"========================================================\n");

			  local_base_pair_step_out(	bAvg, "Avg_bp_step_param", property, ComName, nframe, bps, max_num_bps,fn_cum_data[eBasePairs],fn_cum_data[eLBPS], oenv);

			  fprintf(stderr,"=======================END==============================\n");
		  }
		  if(bLBPHM)	{
			  const char *property[6]=	{ "X-displacement", "Y-displacement", "H-Rise", "Inclination", "Tip", "H-twist" };
			  fprintf(stderr,"\n\n========================================================\n");
			  fprintf(stderr,"             Local Base-Steps Helical Parameters                 \n");
			  fprintf(stderr,"========================================================\n");

			  local_base_pair_step_out(	bAvg, "Avg_bp_helical_param", property, ComName, nframe, bps, max_num_bps,fn_cum_data[eBasePairs],fn_cum_data[eLBPH], oenv);

			  fprintf(stderr,"=======================END==============================\n");
		  }


	  }

  }

  fprintf(stdout, "Thanks for using do_x3dna!!!\n");
  fprintf(stdout, "\n++++ PLEASE READ AND CITE THE FOLLOWING REFERENCE ++++\n");
  fprintf(stderr, "-------- -------- ------------------- -------- ----------\n");
  fprintf(stderr, "Xiang-Jun Lu & Wilma K. Olson (2003)\n");
  fprintf(stderr, "3DNA: a software package for the analysis, rebuilding and visualization\n");
  fprintf(stderr, "of three-dimensional nucleic acid structures.\n");
  fprintf(stderr, "Nucleic Acids Res. 31(17), 5108-21.\n");
  fprintf(stderr, "-------- -------- ------------------- -------- ----------\n");

  return 0;
}


int main(int argc, char *argv[])
{
  gmx_run_cmain(argc, argv, &gmx_3dna);
  return 0;
}
