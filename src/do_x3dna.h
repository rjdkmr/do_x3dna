/*
 * This file is part of do_x3dna
 *
 * Author: Rajendra Kumar
 * Copyright (C) 2014, 2015  Rajendra Kumar
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

#ifndef DO_X3DNA_H_
#define DO_X3DNA_H_

enum { eBasePairs, eHbond, eLBP, eLBPS, eLBPH, eHelAxis, eMgroove, eHelixRad, eBBnDihedral, eSugarConf };

typedef struct {
	int bp1, bp2;
	double **property;
	double *avg, *std;
} t_local_bp_param;

typedef struct {
	int bp11,bp12, bp21, bp22;
	double **property;
	double *avg, *std;
} t_local_bp_step_param;


int calculate_avg_std(double *data, int n, double filter, double *avg, double *std);
int local_base_pair_step_out(gmx_bool bAvg, char *fn_avg_out, char *property[], char *ComName[], int nframe, int **max_bps, int max_num_bps, char *fn_inp_base_pair,char *fn_inp_lbps, const output_env_t oenv);
int local_base_pair_out(gmx_bool bAvg, char *ComName[], int nframe, int **max_bp, int max_num_bp, char *fn_inp_base_pair,char *fn_inp_lbp, const output_env_t oenv);
void hbond_process(int nframe, int **max_bp, int max_num_bp, char *fn_inp_base_pair,char *fn_inp_hbond,int NFILE, char *ComName[],t_filenm fnm[], const output_env_t oenv);
int** get_max_base_pairs(char *filename, int *num_bp);
int** get_max_base_pairs_step(char *filename, int *num_bp);
void write_time(real t, FILE *f_cum_data[]);
int gmx_3dna(int argc,char *argv[]);


#endif /* DO_X3DNA_H_ */
