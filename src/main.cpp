#include <iostream>
#include <direct.h>
#include "CDApp.h"
#include <stdio.h>
#include <stdlib.h>
#include "opencv/highgui.h"
#include "opencv2/opencv.hpp"

#include "stasm/stasm_lib.h"

#include "align/align.h"


using namespace cv;
using namespace std;

//函数功能：根据vec中存储的点的坐标拟合曲线；
//vec为为存储点坐标的容器，index为存储拟合好曲线的数组，len既为拟合曲线的次数，也为数组index的长度；
void fittingCurve(vector<CvPoint2D64f> vec,double *index,int len)
{
	double *px=new double [vec.size()];
	double *py=new double [len*vec.size()];
	int i=0;
	for(vector<CvPoint2D64f>::iterator itr=vec.begin();itr!=vec.end();++itr)
	{
		px[i]=(*itr).x;
		int j=0;
		while (j<len)
		{
			py[len*i+j]=pow((*itr).y,double(j));
			j++;
		}
		i++;
	}

	CvMat xMat=cvMat(vec.size(),1,CV_64FC1,px);
	CvMat yMat=cvMat(vec.size(),len,CV_64FC1,py);
	CvMat *yTransposedMat=cvCreateMat(yMat.cols,yMat.rows,CV_64FC1);
	cvTranspose(&yMat,yTransposedMat);//求yMat的转置

	double *a=new double [len*len];
	for(int i=0;i<len*len;++i)
	{
		a[i]=0;
	}
	CvMat invMat1=cvMat(len,len,CV_64FC1,a);
	cvGEMM(yTransposedMat,&yMat,1,NULL,0,&invMat1,0);//yMat的转置与yMat矩阵相乘
	cvInvert(&invMat1,&invMat1,0);//求invMat的逆矩阵

	double *b=new double [len];
	for(int i=0;i<len;++i)
	{
		b[i]=0;
	}
	CvMat invMat2=cvMat(len,1,CV_64FC1,b);
	cvGEMM(yTransposedMat,&xMat,1,NULL,0,&invMat2,0);//求yMat的转置矩阵与xMat矩阵相乘


	cvGEMM(yTransposedMat,&xMat,1,NULL,0,&invMat2,0);//求yTransposedMat矩阵与xMat 矩阵的乘积
	CvMat indexMat=cvMat(len,1,CV_64FC1,index);
	cvGEMM(&invMat1,&invMat2,1,NULL,0,&indexMat,0);


	cvReleaseMat(&yTransposedMat);
	delete [] a;
	delete [] b;
	delete [] px;
	delete [] py;
}

int mai1n()
{
	vector<CvPoint2D64f> vec;
	int x[]={406,402,401,400,398,397,397,396,396,395,394,393,388,391,391,391};
	int y[]={402,428,454,480,506,532,558,584,610,636,662,688,1078,1104,1130,1156};
	for (int i=0;i<16;i++)
	{
		vec.push_back(cvPoint2D64f(x[i],y[i]));
	}
	double index[3];
	fittingCurve(vec,index,3);
	for (int i=0;i<3;++i)
	{
		printf("%f\n", index[i]);
	}
	cv::waitKey();
	return 0;
}

//三次贝塞尔曲线  
float bezier3funcX(float uu,CvPoint *controlP){  
	float part0 = controlP[0].x * uu * uu * uu;  
	float part1 = 3 * controlP[1].x * uu * uu * (1 - uu);  
	float part2 = 3 * controlP[2].x * uu * (1 - uu) * (1 - uu);  
	float part3 = controlP[3].x * (1 - uu) * (1 - uu) * (1 - uu);     
	return part0 + part1 + part2 + part3;   
}      
float bezier3funcY(float uu,CvPoint *controlP){  
	float part0 = controlP[0].y * uu * uu * uu;  
	float part1 = 3 * controlP[1].y * uu * uu * (1 - uu);  
	float part2 = 3 * controlP[2].y * uu * (1 - uu) * (1 - uu);  
	float part3 = controlP[3].y * (1 - uu) * (1 - uu) * (1 - uu);     
	return part0 + part1 + part2 + part3;   
} 

void createCurve(CvPoint *originPoint,int originCount,vector<CvPoint> &curvePoint){  
	//控制点收缩系数 ，经调试0.6较好，CvPoint是opencv的，可自行定义结构体(x,y)  
	float scale = 0.6;  
	vector<CvPoint>  midpoints(originCount);

	//生成中点       
	for(int i = 0 ;i < originCount ; i++){      
		int nexti = (i + 1) % originCount;  
		midpoints[i].x = (originPoint[i].x + originPoint[nexti].x)/2.0;  
		midpoints[i].y = (originPoint[i].y + originPoint[nexti].y)/2.0;  
	}      

	//平移中点  
	vector<CvPoint> extrapoints(originCount * 2);   
	for(int i = 0 ;i < originCount ; i++){  
		int nexti = (i + 1) % originCount;  
		int backi = (i + originCount - 1) % originCount;  
		CvPoint midinmid;  
		midinmid.x = (midpoints[i].x + midpoints[backi].x)/2.0;  
		midinmid.y = (midpoints[i].y + midpoints[backi].y)/2.0;  
		int offsetx = originPoint[i].x - midinmid.x;  
		int offsety = originPoint[i].y - midinmid.y;  
		int extraindex = 2 * i;  
		extrapoints[extraindex].x = midpoints[backi].x + offsetx;  
		extrapoints[extraindex].y = midpoints[backi].y + offsety;  
		//朝 originPoint[i]方向收缩   
		int addx = (extrapoints[extraindex].x - originPoint[i].x) * scale;  
		int addy = (extrapoints[extraindex].y - originPoint[i].y) * scale;  
		extrapoints[extraindex].x = originPoint[i].x + addx;  
		extrapoints[extraindex].y = originPoint[i].y + addy;  

		int extranexti = (extraindex + 1)%(2 * originCount);  
		extrapoints[extranexti].x = midpoints[i].x + offsetx;  
		extrapoints[extranexti].y = midpoints[i].y + offsety;  
		//朝 originPoint[i]方向收缩   
		addx = (extrapoints[extranexti].x - originPoint[i].x) * scale;  
		addy = (extrapoints[extranexti].y - originPoint[i].y) * scale;  
		extrapoints[extranexti].x = originPoint[i].x + addx;  
		extrapoints[extranexti].y = originPoint[i].y + addy;  

	}      

	CvPoint controlPoint[4];  
	//生成4控制点，产生贝塞尔曲线  
	for(int i = 0 ;i < originCount ; i++){  
		controlPoint[0] = originPoint[i];  
		int extraindex = 2 * i;  
		controlPoint[1] = extrapoints[extraindex + 1];  
		int extranexti = (extraindex + 2) % (2 * originCount);  
		controlPoint[2] = extrapoints[extranexti];  
		int nexti = (i + 1) % originCount;  
		controlPoint[3] = originPoint[nexti];      
		float u = 1;  
		while(u >= 0){  
			int px = bezier3funcX(u,controlPoint);  
			int py = bezier3funcY(u,controlPoint);  
			//u的步长决定曲线的疏密  
			u -= 0.005;  
			CvPoint tempP = cvPoint(px,py);  
			//存入曲线点   
			curvePoint.push_back(tempP);  
		}      
	}  
}  

float polynomial(float x, Mat &factor) {
	float r = 0.0;
	for (int i = 0 ; i < factor.rows; i++) {
		r += factor.at<float>(i, 0) * pow(x, i);
	}
	return r;
}

void polyfit(const Mat& src_x, const Mat& src_y, Mat& dst, int order)
{
	printf("order is %d", order);
	CV_Assert((src_x.rows>0)&&(src_y.rows>0)&&(src_x.cols==1)&&(src_y.cols==1)
		&&(dst.cols==1)&&(dst.rows==(order+1))&&(order>=1));
	Mat X;
	X = Mat::zeros(src_x.rows, order+1,CV_32FC1);
	Mat copy;
	for(int i = 0; i <=order;i++)
	{
		copy = src_x.clone();
		pow(copy,i,copy);
		Mat M1 = X.col(i);
		copy.col(0).copyTo(M1);
	}
	Mat X_t, X_inv;
	transpose(X,X_t);
	Mat temp = X_t*X;
	Mat temp2;
	invert (temp,temp2);
	Mat temp3 = temp2*X_t;
	Mat W = temp3*src_y;
	W.copyTo(dst);
}

#define __global_parameters__
///////////////////////////////////////////////////
float normalized_landmarks[2*stasm_NLANDMARKS];
float tex_coords[2*ALIGN_VERTEXS_NUM];
float model_scales[3];
float model_center[3];
float image_w, image_h;
cv::Mat img_tex;
CDFaceWindow* face_window;

int maptex_callback()
{
	//printf("\n maptex_callback\n");
	static bool first = true; 
	if (first)
	{
		printf("\n maptex_callback first\n");

		face_window->make_current();
		glEnable(GL_TEXTURE_2D);
		//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, 0x2003);	//GL_REPLACE=0x2003
		GLuint idtex;
		glGenTextures(1, &idtex);
		glBindTexture( GL_TEXTURE_2D, idtex ); 
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		//glTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, 512, 512, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, (img_tex.data));
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 512, 512, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, (img_tex.data));

		first = !first;
	}

	float r_hw = image_h/image_w;
	float r_yx = model_scales[1]/model_scales[0];
	//printf("(h,w)= (%g,%g)  r_yx= %g\n", image_h, image_w, r_yx);

	glScalef(1, r_hw, 1);
 	
	//background
	//glPushMatrix();
	glBegin(GL_QUADS);
 	{
 		glTexCoord2f(0,1);		glVertex3f(-1,-1, -10);
 		glTexCoord2f(1,1);		glVertex3f( 1,-1, -10);
 		glTexCoord2f(1,0);		glVertex3f( 1, 1, -10);
 		glTexCoord2f(0,0);		glVertex3f(-1, 1, -10);
 	}
 	glEnd();
	//glPopMatrix();
	
	glScalef(1, r_yx, 1);
	glTranslatef((model_center[0]-0.5)*2, (model_center[1]-0.5)*2, 0);

	glEnableClientState( GL_TEXTURE_COORD_ARRAY );
	glTexCoordPointer(2, GL_FLOAT, 0, tex_coords );


	//glDisableClientState( GL_TEXTURE_COORD_ARRAY );

	return 0;
}

#define __align_callback_begin__
int align_callback(int argc, const char** argv, CDFaceWindow* fw)
{
	face_window = fw;
	
	char path0[100] = "../data/201.jpg";
	
	char* path;
	if (argc>1) path = (char*)argv[1]; else path = path0;

	cv::Mat img0(cv::imread(path, CV_LOAD_IMAGE_COLOR));
	if (!img0.data)
	{
		printf("Cannot load %s\n", path);
		exit(1);
	}
	cv::Mat_<unsigned char> img1;
	cv::cvtColor(img0, img1, CV_RGB2GRAY);
	cv::Mat_<unsigned char> img;
	
	
	img1.copyTo(img);
	int foundface;
	float landmarks[2 * stasm_NLANDMARKS]; // x,y coords (note the 2)

	if (!stasm_search_single(&foundface, landmarks,
		(const char*)img.data, img.cols, img.rows, path, "../data"))
	{
		printf("Error in stasm_search_single: %s\n", stasm_lasterr());
		exit(1);
	}

	img1.copyTo(img);
	if (!foundface)
		printf("No face found in %s\n", path);
	else
	{
		// draw the landmarks on the image as white dots (image is monochrome)
		stasm_force_points_into_image(landmarks, img.cols, img.rows);
		for (int i = 0; i < stasm_NLANDMARKS; i++)
			img(cvRound(landmarks[i*2+1]), cvRound(landmarks[i*2])) = 255;
	}

	cv::imwrite("minimal.bmp", img);
	cv::namedWindow("stasm minimal", cv::WINDOW_NORMAL);
	cv::imshow("stasm minimal", img);


/*	vector<CvPoint> result_p;
	CvPoint input_p[16];
	for (int i = 0; i < 16; i++) {
		input_p[i].x = landmarks[i * 2 + 1];
		input_p[i].y = landmarks[i * 2];
	}

	createCurve(input_p, 16, result_p);

	for (int i = 0; i < result_p.size(); i++) {
		//	img(result_p[i].x, result_p[i].y) = 255;
	}
*/


#define __do_alignment__
	//对准开始//////////////////////////////////////

	for (int i=0; i<2*stasm_NLANDMARKS; i++)
	{
		normalized_landmarks[i] = landmarks[i];
	}

	normalize_point2s(normalized_landmarks, stasm_NLANDMARKS, img.cols, img.rows);
	if (!align_candide3_to_normalized_stasm(normalized_landmarks, tex_coords, model_scales, model_center))
	{
		printf("Error in align_candide3_to_normalized_stasm()\n");
		exit(1);
	}

	flip_normalized_point2s(tex_coords, ALIGN_VERTEXS_NUM);
	image_w = img.cols;
	image_h = img.rows;

	//对准结束///////////////////////////////////////
	

	img1.copyTo(img);
	float w = img.cols;
	float h = img.rows;
	for (int i = 0; i < ALIGN_VERTEXS_NUM; i++)
	{
		float x = tex_coords[i*2];
		float y = tex_coords[i*2+1];
		x = (x)*(w-1);
		y = (y)*(h-1); 
		img(cvRound(y), cvRound(x)) = 255;
	}
	cv::imwrite("candide_dot.bmp", img);
	cv::namedWindow("candide_dot", cv::WINDOW_NORMAL);
	cv::imshow("candide_dot", img);


#define __texure_mapping_begin__
	//映射纹理开始//////////////////////////////////////////
	
	//CDFaceWindow:Fl_Gl_Window 已构造完毕, 但是 context==0!!!
	{
		cv::resize(img0, img_tex, cv::Size(512,512), 0, 0, CV_INTER_LANCZOS4);
//		cv::namedWindow("imag_tex", cv::WINDOW_AUTOSIZE);
//		cv::imshow("imag_tex", img_tex);
		cv::namedWindow("original", cv::WINDOW_AUTOSIZE);
		cv::imshow("original", img0);
	}

	//加载纹理对象
	//映射纹理坐标
	face_window->faceData.set_maptex_callback(maptex_callback);

	//映射纹理结束//////////////////////////////////////////


	return 0;
}

int main(int argc, const char * argv[])
{
	//candide 模型读取、渲染部分

	CDApp candideApp(argc, argv);

	return candideApp.run(align_callback);

	
	cv::waitKey();
	return 0;
}
