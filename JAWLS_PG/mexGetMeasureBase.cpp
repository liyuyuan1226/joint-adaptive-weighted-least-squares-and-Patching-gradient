#include<mex.h>
#include<math.h>

void GetBlock(int cur_i, int cur_j, int m, int n, int winR, double *Image, double *Block, double *Coor, double *CoorX, double *CoorY, int InvFlag);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *CR, *CG, *CB;
    //double sigmaC;
    int m, n, rSmooth, winSizeSmooth;
    int patchR;//patchR是围绕点求周围像素点加和的patch，而rsmooth是在原图上取得小窗口，与GetRelSmooth中的winR是一样的
    CR = mxGetPr(prhs[0]);
    //CG = mxGetPr(prhs[1]);
    //CB = mxGetPr(prhs[2]);
    //sigmaC = mxGetScalar(prhs[3]);
    rSmooth = (int)mxGetScalar(prhs[1]);
    //ColorWeightLUT = mxGetPr(prhs[5]);
	patchR = (int)mxGetScalar(prhs[2]);
	//自加	
    //double *ColorWeightLUT2;
    //ColorWeightLUT2 = mxGetPr(prhs[6]);
	//over
    m = mxGetM(prhs[0]); // row number
    n = mxGetN(prhs[0]); // column number
    
    winSizeSmooth = 2*rSmooth + 1; //window size of the smoothness term
    

    
    /************* compute weight for the smoothness term **********************/
    //plhs[0] = mxCreateSparse(m*n, m*n, winSizeSmooth*winSizeSmooth*m*n, mxREAL); 
    //double *WeightColorPr = mxGetPr(plhs[0]); // the output depth weight  sparse matrix for data term
    //mwIndex *WeightColorIr = mxGetIr(plhs[0]);
    //mwIndex *WeightColorJc = mxGetJc(plhs[0]);
	//自加�?	
	plhs[0] = mxCreateSparse(m*n, m*n, winSizeSmooth*winSizeSmooth*m*n, mxREAL); 
    double *RSPr = mxGetPr(plhs[0]); // the output depth weight  sparse matrix for data term
    mwIndex *RSIr = mxGetIr(plhs[0]);
    mwIndex *RSJc = mxGetJc(plhs[0]);
	//over
    
    /************ compute weight *************/
    int i, j, s, t, counter;
	int i_offset_Center, j_offset_Center, i_offset_Block, j_offset_Block;
	double tempCenter, tempBlock, RelSmooth;
    int InvFlag = -1; // if the pixle in the window is outside the image, set the corresponding coordinate in Coor as InvFlag.
    //double diffR, diffG, diffB;
    
    double *BlockR = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL));
    double *BlockG = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL));
    double *BlockB = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL));
    double *Coor = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL)); // to label the coordinate of the pixel in the sparse matrix
    double *CoorX = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL)); 
	double *CoorY = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL)); 
    
	RSJc[0] = 0;
    //WeightColorJc[0] = 0;
	//WeightColorJc2[0] = 0;//自加
    counter = 0;
    
    for(j=0;j<n;j++)
    {
        for(i=0;i<m;i++)
        {
            GetBlock(i, j, m, n, rSmooth, CR, BlockR, Coor, CoorX, CoorY, InvFlag);
            GetBlock(i, j, m, n, rSmooth, CG, BlockG, Coor, CoorX, CoorY, InvFlag);
            GetBlock(i, j, m, n, rSmooth, CB, BlockB, Coor, CoorX, CoorY, InvFlag);
			
			tempCenter = 0;
			for(p=-patchR;p<=patchR;p++)
						{
							i_offset_Center = i + p ;// -1 because the coordinate in MinCoorX is in real world coordinate
							if(i_offset_Center<0) i_offset_Center = 0;
							if(i_offset_Center>m-1) i_offset_Center = m-1;
				for(q=-patchR;q<=patchR;q++)
							{	
								j_offset_Center = j + q ;
								if(j_offset_Center<0) j_offset_Center = 0;
								if(j_offset_Center>n-1) j_offset_Center = n-1;
							tempCenter = tempCenter + CR[j_offset_Center*m + i_offset_Center];
							}
						}
			
            
            for(t=0;t<winSizeSmooth;t++)
            {
                for(s=0;s<winSizeSmooth;s++)
                {
                    if(Coor[t*winSizeSmooth + s]!=InvFlag)
                    { 
                        
						tempBlock = 0;
            
						for(p=-patchR;p<=patchR;p++)
						{
							/*i_offset_Center = i + p ;// -1 because the coordinate in MinCoorX is in real world coordinate
							if(i_offset_Center<0) i_offset_Center = 0;
							if(i_offset_Center>m-1) i_offset_Center = m-1;*/
                
							i_offset_Block = CoorX[t*winSizeSmooth + s] + p ;
							if(i_offset_Block<0) i_offset_Block = 0;
							if(i_offset_Block>m-1) i_offset_Block = m-1;
							
                
							for(q=-patchR;q<=patchR;q++)
							{	
								/*j_offset_Center = j + q ;
								if(j_offset_Center<0) j_offset_Center = 0;
								if(j_offset_Center>n-1) j_offset_Center = n-1;*/
                
								j_offset_Block = CoorY[t*winSizeSmooth + s] + q ;
								if(j_offset_Block<0) j_offset_Block = 0;
								if(j_offset_Block>n-1) j_offset_Block = n-1;
                
                    
								//tempCenter = tempCenter + CR[j_offset_Center*m + i_offset_Center];
								tempBlock = tempBlock + CR[j_offset_Block*m + i_offset_Block];
                       
							}
						}
            
						RelSmooth = tempCenter/tempBlock;
					
						RSPr[counter]	= RelSmooth;
						RSIr[counter] = (mwIndex)Coor[t*winSizeSmooth + s];	
								
						//diffR = CR[j*m + i] - BlockR[t*winSizeSmooth + s];
                        //diffG = CG[j*m + i] - BlockG[t*winSizeSmooth + s]; 
                        //diffB = CB[j*m + i] - BlockB[t*winSizeSmooth + s]; 
						
                        
                        //WeightColorPr[counter] = ColorWeightLUT[int(diffR*diffR + diffG*diffG + diffB*diffB)];
						
                        
                  
                        //WeightColorIr[counter] = (mwIndex)Coor[t*winSizeSmooth + s];
						

                        counter = counter + 1;
                    }
                }
            }
            //WeightColorJc[j*m + i + 1] = counter;
			RSJc[j*m + i + 1] = counter;
					
        }
        
    }
   
}


/////////////////////////////////////////////   functions ///////////////////////////////////
void GetBlock(int cur_i, int cur_j, int m, int n, int winR, double *Image, double *Block, double *Coor, double *CoorX, double *CoorY, int InvFlag)
{
    int i, j, sMin, sMax, tMin, tMax, winSize;
    winSize = 2*winR + 1;
    
    for(i=0;i<winSize;i++)
    {
        for(j=0;j<winSize;j++)
        {
            Block[j*winSize + i] = 0;
            Coor[j*winSize + i] = InvFlag;
			CoorX[j*winSize + i] = InvFlag;
			CoorY[j*winSize + i] = InvFlag;
        }
    }
    
    if(cur_i < winR) sMin = winR - cur_i; 
    else sMin = 0;
    
    if(cur_i + winR > m - 1) sMax = winR + m - 1 - cur_i;
    else sMax = winSize - 1;
    
    if(cur_j < winR) tMin = winR - cur_j;
    else tMin = 0;
    
    if(cur_j + winR > n - 1) tMax = winR + n - 1 - cur_j;
    else tMax = winSize - 1;
    
    for(i=sMin;i<=sMax;i++)
    {
        for(j=tMin;j<=tMax;j++)
        {
            Block[j*winSize + i] = Image[(cur_j - winR + j)*m + cur_i - winR + i];
            Coor[j*winSize + i] = (cur_j - winR + j)*m + cur_i - winR + i;//以小窗口为标准的坐标点对应原深度图的坐标点
			CoorX[j*winSize + i] = cur_i - winR + i;
			CoorY[j*winSize + i] = cur_j - winR + j;
        }
    }

}

double max(double x ,double y ,double z)//自加
{
	double max;
	if(x<y)
		max=y;
	else
		max=x;
	if(max<z)
		max=z;
	return max;
} 