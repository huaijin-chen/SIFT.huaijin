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
#include "corex.h"

#ifndef _RANSAC_
#define _RANSAC_ 

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
void RANSAC_Homography ( int mun, CvPoint2D64f *m1, CvPoint2D *m2, CvMat *H, CvMat *inliers_mask);

void 
#endif

