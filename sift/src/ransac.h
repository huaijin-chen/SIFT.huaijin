/*
 * =====================================================================================
 *
 *       Filename:  ransac.h
 *
 *    Description:  estimate 2D homography matrix by RANSAC 
 *
 *        Version:  1.0
 *        Created:  23/08/12 10:32:44
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Huaijin(Chao Chen) (GUN), huaijin511@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */

#ifndef _RANSAC_
#define _RANSAC_ 

#define FALSE 0
#define TRUE 1
#include "cv.h"
#include "cxcore.h"
#include <stdlib.h>
#include <stdio.h>
/*******************************************************************
 *	RANSAC algorithm
 *	intput:
 *			num   : the number of points pairs
 *			m1,m2 : pts pairs
 *	
 *	output:
 *			H     : best homography Matrix
 *			inlier_mask (indicate inlier pts pairs in (m1, m2) as 1; outlier: 0)
 *
 *
 ******************************************************************/
void RANSAC_Homography ( int mun, CvPoint2D64f *m1, CvPoint2D64f *m2, CvMat *H, CvMat *inliers_mask);

void compute_Homograpy_matrix (int n, CvPoint2D64f *pts1, CvPoint2D64f *pts2, CvMat *H );
int compute_Number_Inliers( int num, CvPoint2D64f *pts1,  CvPoint2D64f *pts2, CvMat *H, CvMat *inlier_mask, double *dist_std );
void Normalization(int num, CvPoint2D64f *p, CvMat *T);
int isColinear( int num, CvPoint2D64f *p );
#endif

