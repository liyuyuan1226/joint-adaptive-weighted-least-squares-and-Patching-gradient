#include <mex.h>

void GetBlock(int cur_i, int cur_j, int m, int n, int winR, double *Image, double *Block, double *Coor, int InvFlag);

int GetBlockFromSparse(int winSize, double *WeightPr, int WeightPrCount, double *FlagMask, double *WeightBlock, int InvFlag);
/////////////////////////////////////////

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *DepthUpd;
    double *DepthWeightSmoothPr, *ColorWeightPr;
    int m, n, rSmooth, winSizeSmooth;
    
    DepthUpd = mxGetPr(prhs[0]);
    DepthWeightSmoothPr = mxGetPr(prhs[1]);
    ColorWeightPr = mxGetPr(prhs[2]);
    
    m = (int)mxGetScalar(prhs[3]);
    n = (int)mxGetScalar(prhs[4]);
    rSmooth = (int)mxGetScalar(prhs[5]);
	
	//自加：
	double *RScolor, *RS, *ColorWeightPr2, *RSColorMeasurePr, *RSDepthMeasurePr; //double *RSColorMeasure_RPr, *RSColorMeasure_GPr, *RSColorMeasure_BPr;
	RScolor = mxGetPr(prhs[6]);
	RS = mxGetPr(prhs[7]);
	ColorWeightPr2 = mxGetPr(prhs[8]);
	//RSColorMeasure_RPr = mxGetPr(prhs[9]);
	//RSColorMeasure_GPr = mxGetPr(prhs[10]);
	//RSColorMeasure_BPr = mxGetPr(prhs[11]);
	RSColorMeasurePr = mxGetPr(prhs[9])
	RSDepthMeasurePr = mxGetPr(prhs[10]);
	//over
	
	
	
    winSizeSmooth = 2*rSmooth + 1;
    
    
    ////////// Output ////////
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
   
    
    double *WeightSumSmooth = mxGetPr(plhs[0]);
    double *WeightedDepthUpd = mxGetPr(plhs[1]);

    double *DepthWeightSmoothBlock = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL));
    double *DepthUpdBlock = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL));
    double *ColorWeightBlock = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL));
    double *FlagMaskSmooth = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL));
	
	//自加：
	double *RScolorBlock = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL));
	double *RSBlock = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL));
	double *FlagMaskSmooth1 = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL));//这个变量没啥意义，就是为了占个坑，防止出错。
	double *ColorWeightBlock2 = mxGetPr(mxCreateDoubleMatrix(winSizeSmooth, winSizeSmooth, mxREAL));
	//over
    
    int DepthSmoothPrCount = 0, ColorPrCount = 0, counter = 0, InvFlag = -1;
    double WeightSum, WeightedDepthSum;
    int i, j, s, t;
	int ColorPrCount2 = 0;//自加
    int RSColorMeasurePrCount = 0;//自加
	int RSDepthMeasurePrCount = 0;//自加
    for(j=0;j<n;j++)
    {
        for(i=0;i<m;i++)
        {
            GetBlock(i, j, m, n, rSmooth, DepthUpd, DepthUpdBlock, FlagMaskSmooth, InvFlag);//针对于图像尺寸大小的用GetBlock,而对于把窗口都扩展出来的DepthWeightSmoothPr与ColorWeightPr是用GetBlockFromSparse
            //自加：
			GetBlock(i, j, m, n, rSmooth, RScolor, RScolorBlock, FlagMaskSmooth1, InvFlag);//FlagMaskSmooth改为FlagMaskSmooth1是因为不想因为多加的这个函数改变FlagMaskSmooth的值，虽然不知道即便用FlagMaskSmooth，是否真的会改变
			GetBlock(i, j, m, n, rSmooth, RS, RSBlock, FlagMaskSmooth1, InvFlag);
			//over
            DepthSmoothPrCount = GetBlockFromSparse(winSizeSmooth, DepthWeightSmoothPr, DepthSmoothPrCount, FlagMaskSmooth, DepthWeightSmoothBlock, InvFlag);
            ColorPrCount = GetBlockFromSparse(winSizeSmooth, ColorWeightPr, ColorPrCount, FlagMaskSmooth, ColorWeightBlock, InvFlag);
            ColorPrCount2 = GetBlockFromSparse(winSizeSmooth, ColorWeightPr2, ColorPrCount2, FlagMaskSmooth, ColorWeightBlock2, InvFlag);//自加
			RSColorMeasurePrCount = GetBlockFromSparse(winSizeSmooth, RSColorMeasurePr, RSColorMeasurePrCount, FlagMaskSmooth, RSColorMeasureBlock, InvFlag);//自加
            RSDepthMeasurePrCount = GetBlockFromSparse(winSizeSmooth, RSDepthMeasurePr, RSDepthMeasurePrCount, FlagMaskSmooth, RSDepthMeasureBlock, InvFlag);//自加
            /************ compute the weight sum and weighted depth sum for the smoothness term *******/
            
            WeightSum = 0; WeightedDepthSum = 0;
            for(t=0;t<winSizeSmooth;t++)
            {
                for(s=0;s<winSizeSmooth;s++)
                {
                    if(FlagMaskSmooth[t*winSizeSmooth + s] != InvFlag)
                    {
						//自加：
						if(((RScolorBlock[t*winSizeSmooth + s]<=0.96)&&(RSBlock[t*winSizeSmooth + s]>=0.96))
							||((RSColorMeasureBlock[t*winSizeSmooth + s]<0.9)&&(RSDepthMeasureBlock[t*winSizeSmooth + s]>0.9)))
						{
								WeightSum = WeightSum + DepthWeightSmoothBlock[t*winSizeSmooth + s];
								WeightedDepthSum = WeightedDepthSum + DepthUpdBlock[t*winSizeSmooth + s]
                                *DepthWeightSmoothBlock[t*winSizeSmooth + s];
						}
						//over
						else if(((RScolorBlock[t*winSizeSmooth + s]>0.9)&&(RSBlock[t*winSizeSmooth + s]<0.96))
							||((RSColorMeasureBlock[t*winSizeSmooth + s]>0.9)&&(RSDepthMeasureBlock[t*winSizeSmooth + s]<0.9)))
						{
							WeightSum = WeightSum + DepthWeightSmoothBlock[t*winSizeSmooth + s]*ColorWeightBlock2[t*winSizeSmooth + s];
							WeightedDepthSum = WeightedDepthSum + DepthUpdBlock[t*winSizeSmooth + s]
                            *DepthWeightSmoothBlock[t*winSizeSmooth + s]*ColorWeightBlock2[t*winSizeSmooth + s];
							
						}
						else
						{
							
                        WeightSum = WeightSum + DepthWeightSmoothBlock[t*winSizeSmooth + s]*ColorWeightBlock[t*winSizeSmooth + s];
                        WeightedDepthSum = WeightedDepthSum + DepthUpdBlock[t*winSizeSmooth + s]
                                *DepthWeightSmoothBlock[t*winSizeSmooth + s]*ColorWeightBlock[t*winSizeSmooth + s];
						}
									
								
                    }
                }
            }
            //mexPrintf("%f\n", WeightedDepthSum);
            WeightSumSmooth[j*m + i] = WeightSum;
            WeightedDepthUpd[j*m + i] = WeightedDepthSum;
            
        }
    } 
    
}


/////////////////////////////////////////////////////////

void GetBlock(int cur_i, int cur_j, int m, int n, int winR, double *Image, double *Block, double *Coor, int InvFlag)
{
    int i, j, sMin, sMax, tMin, tMax, winSize;
    winSize = 2*winR + 1;
    
    for(i=0;i<winSize;i++)
    {
        for(j=0;j<winSize;j++)
        {
            Block[j*winSize + i] = 0;
            Coor[j*winSize + i] = InvFlag;
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
            Coor[j*winSize + i] = (cur_j - winR + j)*m + cur_i - winR + i;
        }
    }

}


//////////////////////////


int GetBlockFromSparse(int winSize, double *WeightPr, int WeightPrCount, double *FlagMask, double *WeightBlock, int InvFlag)//FlagMask就是GetBlock里面的Coor,也就是坐标
{
    int i, j, count = 0;
    
    for(j=0;j<winSize;j++)
    {
        for(i=0;i<winSize;i++)
        {
            if(FlagMask[j*winSize + i]==InvFlag)
            {
                WeightBlock[j*winSize + i] = 0;
            }
            else
            {
                WeightBlock[j*winSize + i] = WeightPr[WeightPrCount + count];
                count = count + 1;
            }
        }
    }
    
    return WeightPrCount + count;
    
}
