/*
 * =====================================================================================
 *
 *       Filename:  ransac.c
 *
 *    Description: Using RANSAC algorithm to estimate homography matrix H. 
 *
 *        Version:  1.0
 *        Created:  23/08/12 
 *       Revision:  none
 *       Compiler:  g++
 *          Email:  huaijin511@gmail.com
 *         Author:  Huaijin(Chao Chen) (GUN)
 *   Organization:  Wu Han University (in China)
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include "ransac.h"

#define T_DIST 30

void RANSAC_Homography ( int mun, CvPoint2D64f *m1, CvPoint2D64f *m2, CvMat *H, CvMat *inliers_mask)/*{{{*/
{
	int i, j;
	int N = 1000;
	int s = 4;
	int sample_cnt = 0;
	int numinlier, MAX_num;
	double e, p = 0.99;
	double curr_dist_std, dist_std;
	bool iscolinear;

	CvPoint2D64f *curr_m1 = new CvPoint2D64f[s]; 
	CvPoint2D64f *curr_m2 = new CvPoint2D64f[s]; 

	int *curr_idx = new int[s];

	CvMat *curr_inlier_mask = cvCreateMat(num,1,CV_64FC1); 	//k-by-1 
 	CvMat *curr_H = cvCreateMat(3,3,CV_64FC1);  			//3-by-3
 	CvMat *T1 = cvCreateMat(3,3,CV_64FC1);  				//3-by-3
 	CvMat *T2 = cvCreateMat(3,3,CV_64FC1); 					//3-by-3
 	CvMat *invT2 = cvCreateMat(3,3,CV_64FC1); 				//3-by-3
	CvMat *tmp_pt = cvCreateMat(3,1,CV_64FC1); 				//3-by-1

	
	/* RANSAC algorithm (reject outliers and obtain the best H)  */
	//srand(134);
	MAX_num = -1;
	/* create s random & unique indecs  */
	while( N < sample_cnt )
	{
		iscolinear = true;
		while( iscolinear == true )
		{
			iscolinear = false;
			
			for ( i = 0; i < s; i++ ) 
			{
				curr_ind[i] = rand() % num;
				for ( j = 0; j < i; j++ )
				{
					if ( curr_idx[i] == curr_idx[j] ) 
					{
						iscolinear = true;
						break;
					}/* end for*/

				}/* end for */

				if(iscolinear == true) break;
				curr_m1[i].x = m1[curr_idx[i]].x; 
				curr_m1[i].y = m1[curr_idx[i]].y; 
				curr_m2[i].x = m2[curr_idx[i]].x; 
				curr_m2[i].y = m2[curr_idx[i]].y; 
			}
			/* check whether this points colinear  */
			if( iscolinear == false )
				iscolinear = isColinear( s, curr_m1 );
			
		}/* end while */
		/* noramlizte DLT */
		Noramlization( s, curr_m1, T1 );/* curr_m1 <- T1 * curr_m1 */
		Noramlization( s, curr_m2, T2 );/* curr_m1 <- T2 * curr_m2 */
		computer_homograpy_matrix( s, curr_m1, curr_m2, curr_H );
		cvInvert(T2, invT2);
		cvMatMul(invT2, curr_H, curr_H); // curr_H <- invT2 * curr_H
		cvMatMul(curr_H, T1, curr_H);
		curr_H <- curr_H * T1
		/* Calculate the distance for each putative correspondence
			and compute the number of inliers */
		numinlier = compute_Number_Inliers(num,m1,m2,curr_H,curr_inlier_mask,&curr_dist_std);
		// Update a better H
		if(numinlier > MAX_num || (numinlier == MAX_num && curr_dist_std <dist_std)){
			MAX_num = numinlier;
			cvCopy(curr_H, H);
			cvCopy(curr_inlier_mask, inlier_mask);
			dist_std = curr_dist_std;
		}/* end if */
		/* update number N by Algorithm 4.5 */
		e = 1 - (double)numinlier / (double)num;
		N = (int)(log(1-p)/log(1-pow(1-e,s)));
		sample_cnt++;



	}/* end while */

	delete curr_m1, curr_m2, curr_idx;
	cvReleaseMat(&curr_H);
	cvReleaseMat(&T1);
	cvReleaseMat(&T2);
	cvReleaseMat(&invT2);
	cvReleaseMat(&tmp_pt);
	cvReleaseMat(&curr_inlier_mask);

}/*end function  */
/*}}}*/ 

/****************************************************************
* 
* Compute the homography matrix H 
* solve the optimization problem min ||Ah||=0 s.t. ||h||=1 
* where A is 2n*9, h is 9*1 
* 
*	input :	n (number of pts pairs) 
*			p1, p2 (coresponded pts pairs x and x') 
*
*	output: 3*3 matrix H
* 
* ****************************************************************/
void compute_Homograpy_matrix (int n, CvPoint2D64f *pts1, CvPoint2D64f *pts2, cvMat *H )/*{{{*/
{
	int i;
	CvMat *A = cvCreateMat(2*n, 9, CV_64FC1); 
	CvMat *U = cvCreateMat(2*n, 2*n, CV_64FC1); 
	CvMat *D = cvCreateMat(2*n, 9, CV_64FC1);
	CvMat *V = cvCreateMat(9, 9, CV_64FC1); 

    cvZero(A); 

	/*create coff matrix A*/
	for(i=0; i<n; i++)
	{ 
		// 2*i row 
		cvmSet(A,2*i,3,-p1[i].x); 
		cvmSet(A,2*i,4,-p1[i].y); 
		cvmSet(A,2*i,5,-1); 
		cvmSet(A,2*i,6,p2[i].y*p1[i].x); 
		cvmSet(A,2*i,7,p2[i].y*p1[i].y); 
		cvmSet(A,2*i,8,p2[i].y); 
        // 2*i+1 row 
		cvmSet(A,2*i+1,0,p1[i].x); 
		cvmSet(A,2*i+1,1,p1[i].y); 
		cvmSet(A,2*i+1,2,1); 
		cvmSet(A,2*i+1,6,-p2[i].x*p1[i].x); 
		cvmSet(A,2*i+1,7,-p2[i].x*p1[i].y); 
		cvmSet(A,2*i+1,8,-p2[i].x); 
	}//endfor
	
	//SVD
	cvSVD( A, D, U, V, CV_SVD_U_T|CV_SVD_V_T );
	
	//take the last column of V^T, last row of V
	for(i=0; i<9; i++) 
	  cvmSet(H, i/3, i%3, cvmGet(V, 8, i));
		
	cvReleaseMat(&A); 
	cvReleaseMat(&U); 
	cvReleaseMat(&D); 
	cvReleaseMat(&V);
}
	
/*end function computer_homograpy_matrix */	
/*}}}*/

	
/***********************************************************
 *	comput the numbers of inlier by computing distance under
 *	particular H
 *	distance = d( Hx, x' ) + d( invHx', x )
 *
 *	intput : 
 *
 *		num :number point pairs 
 *		p1,p2 : conresponed pts1 & pts2 x--->x'
 *		H : the homography Matrix
 *
 *	output :
 *		inlier_mask
 *		dist_std : std of the distance among all the inliers
 *	retrun :
 *		number of inliers 
 *
 * ***************************************************/
int compute_Number_Inliers( int num, CvPoint2D64f *pts1,  CvPoint2D64f *pts2, cvMat *inlier_mask, double dist_std ) /*{{{*/
{
	int i, num_inlier;
	double curr_dist, sum_dist, mean_dist;
	CvPoint2D64f  tmp_pt; 
	CvMat *dist = cvCreateMat( num, 1, CV_64FC1 ); 
	CvMat *x  = cvCreateMat( 3, 1, CV_64FC1 ); 
	CvMat *xp = cvCreateMat( 3, 1, CV_64FC1 ); 
	CvMat *pt = cvCreateMat( 3, 1, CV_64FC1 ); 
	CvMat *invH = cvCreateMat( 3, 3, CV_64FC1); 

	cvInvert( H, invH );

	sum_dist = 0;
	num_inlier = 0;
	cvZore( inlier_mask );

	/* check points */

	for ( i = 0; i < num; i++ )
	{
		/* initial  point1 */
		cvmSet( x, 0, 0, pts1[i].x );
		cvmSet( x, 0, 1, pts1[i].y );
		cvmSet( x, 0, 2, 1 );

		/* initial x' */
		cvmSet( xp, 0, 0, pts2[i].x );
		cvmSet( xp, 0, 1, pts2[i].y );
		cvmSet( xp, 0, 2, 1 );

		/* compute distanc  d( Hx, x' ) */

		cvMatMul( H, x, pt );
		tmp_pt.x = (int)(cvmGet(pt,0,0)/cvmGet(pt,2,0)); 
  		tmp_pt.y = (int)(cvmGet(pt,1,0)/cvmGet(pt,2,0)); 
  		curr_dist = pow(tmp_pt.x-pts2[i].x, 2.0) + pow(tmp_pt.y-pts2[i].y, 2.0); 

		/* d(x, invH)*/
		cvMatMul( invH, xp, pt );
		tmp_pt.x = (int)(cvmGet(pt,0,0)/cvmGet(pt,2,0)); 
  		tmp_pt.y = (int)(cvmGet(pt,1,0)/cvmGet(pt,2,0)); 
  		curr_dist = pow(tmp_pt.x-pts1[i].x, 2.0) + pow(tmp_pt.y-pts1[i].y, 2.0); 
		

		if(curr_dist < T_DIST)
		{ 
			num_inlier++; 
			cvmSet(inlier_mask,i,0,1); 
			cvmSet(dist,i,0,curr_dist); 
			sum_dist += curr_dist; 

		}/*endif*/ 

	}/* endfor */

	/* Compute the standard deviation of the distance */
	mean_dist = sum_dist/(double)num_inlier;
	*dist_std = 0;
	for (i = 0; i < num; i++ )
	{
		if(cvmGet(inlier_mask,i,0) == 1)
			*dist_std += pow(cvmGet(dist,i,0)-mean_dist,2.0); 
	}

	*dist_std /= (double) (num_inlier -1); 

	cvReleaseMat(&dist); 
	cvReleaseMat(&x); 
	cvReleaseMat(&xp); 
	cvReleaseMat(&pt); 
	cvReleaseMat(&invH); 

	return num_inlier; 
	
}
	/*}}}*/





/**************************************************************************
 *
 * finding the normalization  x = T*x'
 * where T=
 *			| s 0 tx |
 *			| 0 s ty |
 *			| 0 0 1  |
 *	
 *  
 *
 *
 *************************************************************************/
void Normalization(int num, CvPoint2D64f *p, CvMat *T){ /*{{{*/ 
	double scale;
	double meanx = 0, meany = 0;
	double tx, ty;
	double value;
	int i;
	CvMat *x = cvCreateMat(3,1,CV_64FC1);
	CvMat *xp = cvCreateMat(3,1,CV_64FC1);

	for ( i = 0; i < num; i++ ){
		meanx += p[i].x;
		meany += p[i].y;
	}
	meanx /= (double)num;
	meany /= (double)num;

	value = 0;
	
	for ( i = 0; i < num; i++ ) {
		value = sqrt( pow( p[i].x - meanx, 2.0 ) + pow( p[i].y -meany, 2.0 ) );
	}
	value /= (double)num;
	scale = sqrt( 2 ) / value;
	tx = -scale * meanx;
	ty = -scale * meany;

    /* create T from tx ty scale */
	cvZero( T );
	cvmSet( 0, 0,scale );
	cvmSet( 0, 2, tx );
	cvmSet( 1, 1, scale );
	cvmSet( 1, 2, ty );
	cvmSet( 2, 2, 1.0 );
	
	/* transform x' = T*x */

	for ( i = 0; i <num; i++ ) {
		cvmSet( x, 0, 0, p[i].x );
		cvmSet( x, 1, 0, p[i].y );
		cvmSet( x, 2, 0, 1.0 );
		
		cvMatMul( T, x, xp );

		p[i].x = xp[i].x / cvmGet( xp, 2, 0 );
		p[i].y = xp[i].y / cvmGet( xp, 2, 0 );
	}

	cvReleaseMat( &x );
	cvReleaseMat( &xp );
}
	
/*}}}*/   /* end funchtion Normalization */	


/* 
 * ***  FUNCTION  **************************************************************
 *         Name:  isColinear
 *  Description:  check colinearity of a set of pts
 *       intput:  num & points 
 *       output:  ture if some points colinear, false if not
 *
 * *****************************************************************************
 */
void isColinear( int num, CvPoint2D64f *p ){ /*{{{*/
	int i, j, k;
	double value;
	bool isColinear = false;

	CvMat *pt1 = cvCreateMat(3,1,CV_64FC1);
	CvMat *pt2 = cvCreateMat(3,1,CV_64FC1);
	CvMat *pt3 = cvCreateMat(3,1,CV_64FC1);
	CvMat *line = cvCreateMat(3,1,CV_64FC1); 
	for ( i = 0; i < num-2; i++ ){
		cvmSet(pt1,0,0,p[i].x);
		cvmSet(pt1,1,0,p[i].y);
		cvmSet(pt1,2,0,1);
		for(j=i+1; j<num-1; j++){
			cvmSet(pt2,0,0,p[j].x);
			cvmSet(pt2,1,0,p[j].y);
			cvmSet(pt2,2,0,1);
			/* compute the line connecting pt1 & pt2 */
			cvCrossProduct(pt1, pt2, line);
			for(k=j+1; k<num; k++){
				cvmSet(pt3,0,0,p[k].x);
				cvmSet(pt3,1,0,p[k].y);
				cvmSet(pt3,2,0,1);
				/* check whether pt3 on the line */
				value = cvDotProduct(pt3, line);
				if(abs(value) < 10e-2){
					iscolinear = true;
					break;
				}/* end if */
			}/* end for */
			if(iscolinear == true) break;
		}/* end for */
		if(iscolinear == true) break;
	}/* end for */

	cvReleaseMat(&pt1);
	cvReleaseMat(&pt2);
	cvReleaseMat(&pt3);
	cvReleaseMat(&line);
	
	return iscolinear;
}   

/*}}}*/
















/***********************************************************************
 *	make points pairs unique
 * *********************************************************************/
