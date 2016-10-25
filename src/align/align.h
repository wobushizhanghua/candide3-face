#ifndef __ALIGN_H__  
#define __ALIGN_H__

void normalize_point2s( 
	float point2s[],	// (j,i) ==> (x,y) = (0,0)~(1,1);
	int n,				// number of points = sizeof(point2s)/2
	int w,				// width  = ncols
	int h,				// height = nrows
	bool b_rowcol=false);
void denormalize_point2s( 
	float point2s[],	// (j,i) <== (x,y) = (0,0)~(1,1);
	int n,				// number of points = sizeof(point2s)/2
	int w,				// width  = ncols
	int h,				// height = nrows
	bool b_rowcol=false);
void flip_normalized_point2s(
	float point2s[],	// (x,y)
	int n				// number of points = sizeof(point2s)/2
	);

const int ALIGN_LANDMARKS_NUM = 77; //=stasm_NLANDMARKS
const int ALIGN_VERTEXS_NUM = 113;

bool align_candide3_to_normalized_stasm(
	float normalized_landmarks[ALIGN_LANDMARKS_NUM],	// (x,y) = (0,0)~(1,1)	//
	float candide3_tex_coords[ALIGN_VERTEXS_NUM],		// (u,v) = (0,0)~(1,1)	//texture 2d coordinates
	float candide3_scales[3],							// (sx,sy,sz)			//scale factors
	float candide3_center[2]							// (cx,cy)
	);
#endif //__ALIGN_H__