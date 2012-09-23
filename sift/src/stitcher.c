/*

	Based Rob Hess'SIFT to detect features & and get the homography matrix

	Auther :Huaijin Chen 
    Email: huaijin511@gmail.com


*/

#include "sift.h"
#include "imgfeatures.h"
#include "kdtree.h"
#include "utils.h"
#include "xform.h"
#include "ransac.h"
//#include "stitcher.h"

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>

#include <stdio.h>

/* the maximum number of keypoint NN candidates to check during BBF search */
#define KDTREE_BBF_MAX_NN_CHKS 200

/* threshold on squared ratio of distances between NN and 2nd NN */
#define NN_SQ_DIST_RATIO_THR 0.49

#define MAX_MATCH_POINT 2000
/* functins Yuan Xing */
double max( double n1, double n2 );

double min( double n1,double n2 );

CvSize get_Stitched_Size(CvSize im1_size, CvSize im2_size,  CvMat* Homo_Mat, double XDATA[], double YDATA[] );

int main( int argc, char** argv )
{

	IplImage* img1, * img2, * stacked;
	struct feature* feat1, * feat2, * feat;
	struct feature** nbrs;
	struct kd_node* kd_root;
	CvPoint pt1, pt2;
	double d0, d1;
	int n1, n2, k, i, m = 0;
	int match_cnt = 0;

	//CvMat *im1_mask = cvCreateMat( img1->height, img1->width, CV_64FC1 );
	//CvMat *im2_mask = cvCreateMat( img2->height, img2->width, CV_64FC1 );
  
	//cvSet( im1_mask, cvScalar( 1, 0, 0, 0 ), NULL );
	//cvSet( im2_mask, cvScalar( 1, 0, 0, 0 ), NULL );

	if( argc != 3 )
		fatal_error( "usage: %s <img1> <img2>", argv[0] );
  
	img1 = cvLoadImage( argv[1], 1 );
	if( ! img1 )
		fatal_error( "unable to load image from %s", argv[1] );
	img2 = cvLoadImage( argv[2], 1 );
	if( ! img2 )
		fatal_error( "unable to load image from %s", argv[2] );
	stacked = stack_imgs( img1, img2 );

	fprintf( stderr, "Finding features in %s...\n", argv[1] );
	n1 = sift_features( img1, &feat1 );
	fprintf( stderr, "Finding features in %s...\n", argv[2] );
	n2 = sift_features( img2, &feat2 );
	kd_root = kdtree_build( feat2, n2 );
	for( i = 0; i < n1; i++ )
    {
		feat = feat1 + i;
		k = kdtree_bbf_knn( kd_root, feat, 2, &nbrs, KDTREE_BBF_MAX_NN_CHKS );
		if( k == 2 )
		{
			d0 = descr_dist_sq( feat, nbrs[0] );
			d1 = descr_dist_sq( feat, nbrs[1] );
			if( d0 < d1 * NN_SQ_DIST_RATIO_THR  )
			{
				if( m >= 2000 ) break;
				pt1 = cvPoint( cvRound( feat->x ), cvRound( feat->y ) );
				pt2 = cvPoint( cvRound( nbrs[0]->x ), cvRound( nbrs[0]->y ) );
				pt2.y += img1->height;
				cvLine( stacked, pt1, pt2, CV_RGB(255,0,255), 1, 8, 0 );
				m++;
				feat1[i].fwd_match = nbrs[0];
			 }
		}
		free( nbrs );
    }

	fprintf( stderr, "Found %d total matches\n", m );
	display_big_img( stacked, "Matches" );
	cvWaitKey( 0 );
	/*********************************************************************************************************/
     
	CvMat* H;
    IplImage* xformed;
    H = ransac_xform( feat1, n1, FEATURE_FWD_MATCH, lsq_homog, 4, 0.01, homog_xfer_err, 3.0, NULL, NULL );
	/* output H */
	printf("xform Homography Matrix is :\n"); 
	printMat( H );

	
	/* get the size of new image  */
	double XDATA[2];
	double YDATA[2];
	CvSize new_size;
	CvSize im1_size = cvGetSize( img1 );
	CvSize im2_size = cvGetSize( img2 );
		
	new_size = get_Stitched_Size( im1_size, im2_size, H, XDATA, YDATA );


	/*declare the mask*/
	CvMat *im1_mask = cvCreateMat( new_size.height, new_size.width, CV_64FC1 );
	CvMat *im2_mask = cvCreateMat( new_size.height, new_size.width, CV_64FC1 );

	CvMat *im1_tempMask = cvCreateMat( im1_size.height, im1_size.width, CV_64FC1 );
	CvMat *im2_tempMask = cvCreateMat( im2_size.height, im2_size.width, CV_64FC1 );
	cvSet( im1_tempMask, cvScalar(1,0,0,0), NULL );
	cvSet( im2_tempMask, cvScalar(1,0,0,0), NULL );

	/*  get translation Matrix for Aligning */
	CvMat *T = cvCreateMat( 3, 3, CV_64FC1 );
	double tx = 0;
	double ty = 0;
	if( XDATA[0] < 0 )
	{
		tx =  -XDATA[0] ;
		printf("tx = %f\n", tx );
	}
	if( YDATA[0] < 0 )
	{
		ty = -YDATA[0];
		printf("ty = %f\n", ty );
	}
	cvmSet( T, 0, 0, 1 );
	cvmSet( T, 0, 2, tx );
	cvmSet( T, 1, 1, 1 );
	cvmSet( T, 1, 2, ty );
	cvmSet( T, 2, 2, 1 );
	printf("T Matrix:\n");
	printMat( T );

	/* Transform and Align image2 */
	cvGEMM( T, H, 1, NULL, 0, H, 0 );
	printMat( H );
	xformed = cvCreateImage( new_size, IPL_DEPTH_8U, 3 );
	cvWarpPerspective( img1, xformed, H, CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS, cvScalarAll( 0 ) );
	cvNamedWindow( "Xformed1", 1 );
	cvShowImage( "Xformed1", xformed );
	cvWaitKey( 0 );
	cvDestroyWindow("Xformed1");
	//cvSaveImage("im2.png", xformed);
	cvWarpPerspective( im1_tempMask, im1_mask, H, CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS, cvScalarAll( 0 ) );
	//cvSaveImage("im1_mask.png", im1_mask);
	cvNamedWindow( "im1_mask", 1 );
	cvShowImage( "im1_mask", im1_mask );
	cvWaitKey( 0 );
	cvDestroyWindow("im1_mask");

	/* Align image1 to bound */
	cvWarpPerspective( im2_tempMask, im2_mask, T, CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS, cvScalarAll( 0 ) );
	//cvSaveImage( "im2_mask.png", im2_mask );
	cvNamedWindow( "im12_mask", 1 );
	cvShowImage( "im2_mask", im2_mask );
	cvWaitKey( 0 );
	cvDestroyWindow("im2_mask");
	cvSetImageROI( xformed, cvRect( tx, ty, img2->width, img2->height ) );
	cvCopy( img2, xformed, NULL );
	
	IplImage* huaijin = cvCreateImage( new_size, IPL_DEPTH_8U, 3 );
	cvWarpPerspective( img2, huaijin, T, CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS, cvScalarAll( 0 ) );
	cvNamedWindow( "im12_mask_i", 1 );
	cvShowImage( "im2_mask_i", huaijin);
	cvWaitKey( 0 );
	cvDestroyWindow("im2_mask_i");
	cvResetImageROI( xformed );
	//cvSaveImage( "re.png", xformed );


	/* composite */


	cvNamedWindow( "Xformed", 1 );
	cvShowImage( "Xformed", xformed );
	cvWaitKey( 0 );
	cvDestroyWindow("Xformed");
	/*  */


	cvReleaseImage( &xformed );
	cvReleaseMat( &H );
	cvReleaseMat( &im1_tempMask );
	cvReleaseMat( &im2_tempMask );
	cvReleaseMat( &T );
	cvReleaseMat ( &im1_mask );
	cvReleaseMat ( &im2_mask );
	cvReleaseImage( &stacked );
	cvReleaseImage( &img1 );
	cvReleaseImage( &img2 );
	kdtree_release( kd_root );
	free( feat1 );
	free( feat2 );

		
/****************************************************
 * using RANSAC algorithm get the Homegraphy Matrix
 * get T image1-->image2
 * 
 * n_pts		the number of  points for estimating parameters  
 *
 *
 *
 * create unique indics of matchs : get_randi( int j )
 *
 *
 * ******slove*********** 
 * 1.create the coff matrix
 * 2.Ah=0 -->decomposit A = UDV^T using gsl 
 * 3 
 *      [ v19 v29 v39 ... v99 ]
 * h = -------------------------
 *               v99
 *
 * ***************************************************/
	/*
	CvMat *H1 = cvCreateMat( 3, 3, CV_64FC1 );
	CvMat *inliers_mask = cvCreateMat( m, 1, CV_64FC1 );
	RANSAC_Homography( m, pts1, pts2, H1, inliers_mask );
	printf("my code  Homography Matrix is :\n"); 
	for ( i = 0; i < H->rows; i++ ){
		for( k = 0; k < H->cols; k++ ){
			printf("%f	",cvmGet( H1, i, k ));
		}
		printf("\n");
	}


	cvReleaseMat( &H1 );
	cvReleaseMat( &inliers_mask );*/
/***********************************************
 * composit image1 & image2 
 * 1) transform image2 to image
 *
 * 2)***stitched image bounds****
 * W = max( [size(im1,2) size(im1,2)-XDATA(1) size(im2,2) size(im2,2)+XDATA(1)] );
 * H = max( [size(im1,1) size(im1,1)-YDATA(1) size(im2,1) size(im2,1)+YDATA(1)] );
 *
 * 3)*** Align image1 to bound ***
 *
 * 4)*** Align image2 to bound ***
 * 
 * 5)*** Check  size of bounds *** 
 * 
 * 6)*** combine both images ***
 *		im1_mask
 *		im2_mask
 *		im1_part_mask
 *		im2_part_mask
 *		com_part_mask
 *		stitched_image
 * 
 * 7)****copy im2 transformed to ROI of stitching plan ( just a idear )
 *
 ************************************************/

  void compositImages( IplImage *im1, IplImage *im2, CvMat *H )
  {
	/*  1.create a plan																			*/
	/*  2.transform im2 & im2_mask																*/
	/*	3.emstime translation of im2 --T  cp im1 -->plan with cvRect( tx, ty)					*/
	/*	4.cp tranformed im2 to plan with im2_mask where im2_mask has the same size with plan	*/

  
  }

  return 0;
}
/* 
 * ***  FUNCTION  **************************************************************
 *         Name:  get_Stitched_Size
 *  Description:  get size of new imags 
 *		  input:  im1_size,im2_size, the sizes of images
 *				  Homo_Mat ,homography matrix
 *
 *		output:  the size of new image
 * *****************************************************************************
 */

CvSize get_Stitched_Size(CvSize im1_size, CvSize im2_size,  CvMat* Homo_Mat, double XDATA[], double YDATA[]  )/*{{{*/
{
	int Width  = 0;
	int Height = 0;
	double p1[3] = { 0, 0, 1 };
	double p2[3] = { im2_size.width-1 , 0, 1 };
	double p3[3] = { 0, im2_size.height-1, 1 };
	double p4[3] = { im2_size.width-1, im2_size.height-1, 1 };

	
	CvMat pp1 = cvMat( 3, 1,  CV_64FC1, p1 );
	CvMat pp2 = cvMat( 3, 1,  CV_64FC1, p2 );
	CvMat pp3 = cvMat( 3, 1,  CV_64FC1, p3 );
	CvMat pp4 = cvMat( 3, 1,  CV_64FC1, p4 );
	
	cvGEMM( Homo_Mat, &pp1, 1, NULL, 0, &pp1, 0 );
	cvGEMM( Homo_Mat, &pp2, 1, NULL, 0, &pp2, 0 );
	cvGEMM( Homo_Mat, &pp3, 1, NULL, 0, &pp3, 0 );
	cvGEMM( Homo_Mat, &pp4, 1, NULL, 0, &pp4, 0 );

	/* Normalization pp --->(x y 1) */
	
	
	double L1 = cvmGet( &pp1, 2, 0 );  
	double L2 = cvmGet( &pp2, 2, 0 );  
	double L3 = cvmGet( &pp3, 2, 0 );  
	double L4 = cvmGet( &pp4, 2, 0 );  

	XDATA[0] = min( min( min( cvmGet(&pp1,0,0)/L1, cvmGet(&pp2,0,0)/L2 ), cvmGet(&pp3,0,0)/L3 ),cvmGet(&pp4,0,0 )/L4 ); 
	YDATA[0] = min( min( min( cvmGet(&pp1,1,0)/L1, cvmGet(&pp2,1,0)/L2 ), cvmGet(&pp3,1,0)/L3 ),cvmGet(&pp4,1,0 )/L4 ); 

	XDATA[1] = max( max( max( cvmGet(&pp1,0,0)/L1, cvmGet(&pp2,0,0)/L2 ), cvmGet(&pp3,0,0)/L3 ),cvmGet(&pp4,0,0 )/L4 ); 
	YDATA[1] = max( max( max( cvmGet(&pp1,1,0)/L1, cvmGet(&pp2,1,0)/L2 ), cvmGet(&pp3,1,0)/L3 ),cvmGet(&pp4,1,0 )/L4 ); 
	
	Width  = max( max( max( im1_size.width,  XDATA[1] ), im1_size.width-XDATA[0]   ), XDATA[1]-XDATA[0] ) ;
	Height = max( max( max( im1_size.height, YDATA[1] ), im1_size.height-YDATA[0] ), YDATA[1]-YDATA[0] ) ;

	Width  = cvRound( Width );
	printf("New image's Width is %d\n", Width );
	Height = cvRound( Height );
	printf("New image's Height is %d\n", Height );

	return cvSize( Width, Height );

}/*}}}*/

/*{{{*/ /* utils functions */
double min( double n1,double n2 )
{
	return n1 <= n2 ? n1 : n2;
}
double max( double n1, double n2 )
{ 

	return n1 >= n2 ? n1 : n2;
}

/*}}}*/
