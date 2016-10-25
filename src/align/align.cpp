#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "levmar.h"
#pragma comment(lib, "../src/align/levmar.lib")

#include "align.h"
#include "candide3_model.h"


void normalize_point2s( 
	float point2s[],	// (j,i) ==> (x,y) = (0,0)~(1,1);
	int n,				// number of points = sizeof(points)/2
	int w,				// width  = ncols
	int h,				// height = nrows
	bool b_rowcol)
{
	float x,y;
	for (int k=0; k<n; k++)
	{
		if (b_rowcol)
		{
			x = point2s[k*2+1];
			y = point2s[k*2];
		}else
		{
			x = point2s[k*2];
			y = point2s[k*2+1];
		}
		x = x/(w-1);
		y = (h-1-y)/(h-1);
		assert(0<=x&&x<=1 && 0<=y&&y<=1);
		{
			point2s[k*2] = x;
			point2s[k*2+1] = y;
		}
	}
}
void denormalize_point2s( 
	float point2s[],	// (j,i) <== (x,y) = (0,0)~(1,1);
	int n,				// number of points = sizeof(points)/2
	int w,				// width  = ncols
	int h,				// height = nrows
	bool b_rowcol)
{
	float x,y;
	for (int k=0; k<n; k++)
	{
		{
			x = point2s[k*2];
			y = point2s[k*2+1];
		}
		x = x*(w-1);
		y = (1-y)*(h-1);
		assert(0<=x&&x<=w-1 && 0<=y&&y<=h-1);
		if (b_rowcol)
		{
			point2s[k*2+1] = x;
			point2s[k*2]   = y;
		}else
		{
			point2s[k*2]   = x;
			point2s[k*2+1] = y;
		}
	}
}
void flip_normalized_point2s(
	float point2s[],	// (x,y)
	int n				// number of points = sizeof(point2s)/2
	)
{
	for (int k=0; k<n; k++)
	{
		//point2s[k*2];
		point2s[k*2+1] = 1-point2s[k*2+1];
	}
}


////////////////////////////////////////
// 1-based corr_table
const int corr_table[][3] = { // i_asm, i_candide, weight,
	//outline
	1,		63,		1,//left
	7,		44,		1,//bottom
	13,		30,		1,//right
	//mouth
	60,		90,		1,//left
	63,		8,		1,//top
	66,		89,		1,//right
	75,		42,		1,//bottom
	//left-eye
	31,		57,		1,//right
	33,		53,		1,//top
	35,		54,		1,//left
	37,		58,		1,//bottom
	//right-eye
	41,		24,		1,//left
	43,		20,		1,//top
	45,		21,		1,//right
	47,		25,		1,//bottom
};
const int num_corr = sizeof(corr_table)/sizeof(int[3]);

const char* lm_stop_reason[]={"0",
	"1 - stopped by small ||J^T e||_inf",
	"2 - stopped by small ||dp||_2",
	"3 - stopped by max # iterations",
	"4 - singular matrix. Restart from current p with increased mu",
	"5 - no further error reduction is possible. Restart with increased mu",
	"6 - stopped by small ||e||_2",
	"7 - stopped by invalid (i.e. NaN or Inf) func values, a user error",
};

void transform(double *p, double &x, double &y)
{
	double tx = p[0];
	double ty = p[1];
	double sx = p[2];
	double sy = p[3];
	double r = p[4];
	double xx =   (sx)*x + (-r*sy)*y + tx;
	double yy = (r*sx)*x +    (sy)*y + ty;
	x = xx;
	y = yy;
	//printf("by (%g %g %g %g %g)\n transfer to(%g %g)\n",
	//	 tx, ty, sx, sy, r, xx,yy);
}

float* asm_points=0;

void obj_func(double *p, double *f, int m, int n, void *data)
{
	register int i;

	//printf("obj_func, m=%i, n=%i\n", m,n);
	for(i=0; i<n; i++)
	{
		int i_asm = corr_table[i][0]-1;
		double asm_x = asm_points[i_asm*2];
		double asm_y = asm_points[i_asm*2+1];
		//printf("asm[%i] =(%g %g)\n",i_asm+1,asm_x,asm_y);

		int i_cand = corr_table[i][1]-1;
		double cand_x = candide3_v[i_cand*3];
		double cand_y = candide3_v[i_cand*3+1];
		//printf("cand[%i] =(%g %g)\n",i_cand+1,cand_x,cand_y);
		
		transform(p, cand_x, cand_y);
		//printf("by (%g %g %g %g %g)\n transfer to cand[%i] =(%g %g)\n",
		//	p[0], p[1], p[2], p[3], p[4], i_cand+1,cand_x,cand_y);
		
		double dx = asm_x-cand_x;
		double dy = asm_y-cand_y;
		f[i]= sqrt(dx*dx + dy*dy) *corr_table[i][2];
		//printf("f[%i] = %g\n", i,f[i]);
	}
}

bool align_candide3_to_normalized_stasm(
	float normalized_landmarks[],	// (x,y) = (0,0)~(1,1)	//
	float candide3_tex_coords[],	// (u,v) = (0,0)~(1,1)	//texture 2d coordinates
	float candide3_scales[3],								// (sx,sy,sz)			//scale factors
	float candide3_center[2]								// (cx,cy)
	)
{
	asm_points = normalized_landmarks;
	
	const int max_iter_num = 1000;
	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	opts[0]=LM_INIT_MU;
	opts[1]=1E-15; opts[2]=1E-15; opts[3]=1E-20;
	opts[4]= LM_DIFF_DELTA;		// forward differencing 
	//opts[4]=-LM_DIFF_DELTA;	// central differencing; more accurate but more expensive
	
	const int m= 5;
	const int n= num_corr;
	double p[m]; // parameters
	double lb[m], ub[m]; // bounds
	double f[n]; // evaluations
	double dr = 3.1416/8;
	p[0]= 0.5;	lb[0]=0;	ub[0]=1;	//tx
	p[1]= 0.5;	lb[1]=0;	ub[1]=1;	//ty
	p[2]= 1;	lb[2]=0.1;	ub[2]=10;	//sx
	p[3]= 1;	lb[3]=0.1;	ub[3]=10;	//sy
	p[4]= 0;	lb[4]=-dr;	ub[4]=dr;	//r
	for(int i=0; i<n; i++)
		f[i]=0.0;
	
	assert(asm_points);
	int niter = dlevmar_bc_dif(obj_func, p, f, m, n,
		lb, ub, 0,
		max_iter_num, opts, info, 0, 0, 0);  // dif Jacobian	
	
	printf("\n  %g in %i iter, stop reason: %s\n", info[5], max_iter_num, lm_stop_reason[(int)info[6]]);
	printf("Solution: tx ty sx sy r \n  ");
	for(int i=0; i<m; ++i)
		printf("%.7g ", p[i]);
	printf("\n");
	printf("||e0||, ||e||, ||J^T e||, ||dp||, mu/max[J^T J]_ii, #iter, stop, #f, #J, #lin :\n  ");
	for(int i=0; i<LM_INFO_SZ; ++i)
		printf("%g ", info[i]);
	printf("\n\n");


	//////////////////////////////
	for (int i=0; i<ALIGN_VERTEXS_NUM; i++)
	{
		double cand_x = candide3_v[i*3];
		double cand_y = candide3_v[i*3+1];

		transform(p, cand_x, cand_y);

		candide3_tex_coords[i*2]   = cand_x; 
		candide3_tex_coords[i*2+1] = cand_y;
	}

	candide3_scales[0] = p[2];
	candide3_scales[1] = p[3];
	candide3_scales[2] = 1;

	candide3_center[0] = p[0];
	candide3_center[1] = p[1];

	return (niter>0);
}
