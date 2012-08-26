/*
 * =====================================================================================
 *
 *       Filename:  GAUSSian.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  30/07/12 16:13:34
 *       Revision:  none
 *       Compiler:  gcc
 *          Email:  huaijin511@gmail.com
 *         Author:  Huaijin(Chao Chen) (GUN)
 *   Organization:  Wu Han University (in China)
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <cv.h>
#include <highgui.h>
#include <cxcore.h>
#include <stdio.h>
#include <string>
static IplImage* convert_to_gray32( IplImage* img );

int main(int arg, char** argv)
{
	IplImage* src = cvLoadImage(argv[1]);
	IplImage* src_gray32, *result; 
	//String src_color = (string)src->colorModel;
	//printf("-------%s-------\n", src->colorModel );
	//std::cout<<  "----"<<src_color<<std::endl;
	src_gray32 = convert_to_gray32( src );
	double sigma[4] ={1.6,2.0,2.6,3.2};
	//scanf("Enter sigma = %lf", &sigma )
	int i;
	result = cvCreateImage( cvGetSize(src), IPL_DEPTH_8U, 1 );
	for (i = 0;i < 4; i++ )
	{
		cvSmooth( src_gray32, src_gray32, CV_GAUSSIAN, 0, 0, sigma[i], sigma[i] );
		cvConvertScale(src_gray32, result, 255.0, 0 );
		char fname[15];
		sprintf( fname, "%lf_%s", sigma[i], argv[1] );
		cvSaveImage( fname, result );
		//cvReleaseImage( &result );
	}
	/*
	cvNamedWindow("scaleImage", 1);
	cvShowImage( "scaleImage", src_gray32 );
	if (cvWaitKey(0))
	{
		cvReleaseImage( &src );
		cvReleaseImage( &src_gray32 );
		
	}*/
		cvReleaseImage( &src );
		cvReleaseImage( &src_gray32 );
		cvReleaseImage( &result );

	return 0;
}

static IplImage* convert_to_gray32( IplImage* img )
{
  IplImage* gray8, * gray32;
  
  gray32 = cvCreateImage( cvGetSize(img), IPL_DEPTH_32F, 1 );
  if( img->nChannels == 1 )
	 gray8 = cvCloneImage( img );
  else
    {
      gray8 = cvCreateImage( cvGetSize(img), IPL_DEPTH_8U, 1 );
      cvCvtColor( img, gray8, CV_BGR2GRAY );
    }
  cvConvertScale( gray8, gray32, 1.0 / 255.0, 0 );

  cvReleaseImage( &gray8 );
  return gray32;
}
