//implement of  ApplyMulDPolN class
//
#include "ApplyMulDPolN.h"

template <typename T>
void ApplyMulDPolN::FreeDynArray(T *****arr, int bin1, int bin2, int bin3,int bin4)
{   
	for(int iw=0;iw<bin1;iw++)
	{  
		for(int iq=0;iq<bin2;iq++)
		{  
			for(int it=0;it<bin3;it++) 
			{
				for(int ip=0;ip<bin4;ip++) delete [] arr[iw][iq][it][ip];
			}
		}
	}
};

//create a buffer using a given value
template <typename T>
T ***** ApplyMulDPolN::IniDynamicArray(int bin1, int bin2, int bin3, int bin4,int bin5, T defautval)
{
	T ***** arr;
	arr=new T **** [bin1];
	for(int iw=0;iw<bin1;iw++)
	{  	  
		arr[iw]=new T *** [bin2];			 
		for(int iq=0;iq<bin2;iq++)
		{  
			arr[iw][iq]=new T ** [bin3];					 
			for(int it=0;it<bin3;it++)
			{	
				arr[iw][iq][it]=new T * [bin4];
				for(int ip=0;ip<bin4;ip++)  
				{
					arr[iw][iq][it][ip]=new T [bin5];
					for(int ie=0;ie<bin5;ie++)  arr[iw][iq][it][ip][ie]=defautval;
				}
			}
		}
	}
	return arr;
};

template <typename T>
void ApplyMulDPolN::DebugArray(T ***** arr, int bin1, int bin2, int bin3, int bin4,int bin5)
{	
	for(int iw=0;iw<bin1;iw++)
	{  	  
		cout<<"iw="<<iw<<endl;		 
		for(int iq=0;iq<bin2;iq++)
		{  
			cout<<"iq="<<iq<<endl;						 
			for(int it=0;it<bin3;it++)
			{	
				cout<<"it="<<it<<endl;	
				for(int ip=0;ip<bin4;ip++)  
				{
					cout<<"ip="<<ip<<endl;	
					for(int ie=0;ie<bin5;ie++)  
					{
						cout<<"ie="<<ie<<"  Array_"<<iw<<iq<<it<<ip<<ie
							<<"="<<arr[iw][iq][it][ip][ie]<<endl;
					}
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
ApplyMulDPolN::ApplyMulDPolN(const char* infile):
mPara(0),mFitPara(0),mFitParaFixedFlag(0),mFitParaIndexMap(0),mParaNum(0)
{
#ifdef ApplyMulDPolN_DEBUG 
	if(ApplyMulDPolN_DEBUG>0)
	{
		cout<<">>>>>ApplyMulDPolN::ApplyMulDPolN() started with debug level "<<ApplyMulDPolN_DEBUG<<"<<<<<"<<endl;
	}
#endif

	ReadParaFile(infile);
}

ApplyMulDPolN::~ApplyMulDPolN()
{
#ifdef ApplyMulDPolN_DEBUG
	if(ApplyMulDPolN_DEBUG>=3)
		cout<<"ApplyMulDPolN::~ApplyMulDPolN() "<<endl;
#endif
	FreeDynArray(mPara,MaxOrderA+1,MaxOrderB+1,MaxOrderC+1,MaxOrderD+1);
	delete [] mFitPara;
	delete [] mFitParaFixedFlag;
	for(int ii=0;ii<mParaNum;ii++) delete [] mFitParaIndexMap[ii];
}

void ApplyMulDPolN::ReadParaFile(const char* infile)
{
#ifdef ApplyMulDPolN_DEBUG
	if(ApplyMulDPolN_DEBUG>=3)
		cout<<"ApplyMulDPolN::ReadParaFile()"<<endl;
#endif

	ifstream fin;
	fin.open(infile);

	if(!fin.good())
	{
		cout<<"ApplyMulDPolN::ReadParaFile(): could not open file \""<<infile<<"\", I quit ...\n"<<endl;		
		exit(-1);
	}

	char buf[512];
	char label;  //only accept ABCDE
	int  ia,ib,ic,id,ie;
	int line=0;
	int index=0, pBufIndex=0; 

	//eat the first 8 lines
	for(int i=0;i<8;i++) {fin.getline(buf,511);line++;}

	//line9
	fin>>MaxOrderA>>MaxOrderB>>MaxOrderC>>MaxOrderD>>MaxOrderE;
	//eat the rest of the line
	fin.getline(buf,511);line++; 
#ifdef ApplyMulDPolN_DEBUG
	if(ApplyMulDPolN_DEBUG>=4) cout<<"Line="<<line<<"  buf=\""<<buf<<"\"  index="<<index<<endl;
#endif

	//now setup the buffer
	double defval=0.0; 
	ApplyMulDPolN::mPara = IniDynamicArray(MaxOrderA+1,MaxOrderB+1,MaxOrderC+1,MaxOrderD+1,MaxOrderE+1,defval);
	MaxOrder = max(max(max(MaxOrderA,MaxOrderB),max(MaxOrderC,MaxOrderD)),MaxOrderE);

	double *par = new double [MaxOrder+1];
	int *par_fixedflag = new int [MaxOrder+1];

	int pMaxParaNum=(MaxOrderA+1)*(MaxOrderB+1)*(MaxOrderC+1)*(MaxOrderD+1)*(MaxOrderE+1);	
	mFitPara = new double [pMaxParaNum];
	mFitParaErr = new double [pMaxParaNum];
	mFitParaFixedFlag = new int [pMaxParaNum];
	mFitParaIndexMap = new int * [pMaxParaNum];
	for(int ii=0;ii<pMaxParaNum;ii++) mFitParaIndexMap[ii]=new int [6];

	//now eat the 10th line
	fin.getline(buf,511);line++; 
#ifdef ApplyMulDPolN_DEBUG
	if(ApplyMulDPolN_DEBUG>=4) cout<<"Line="<<line<<"  buf=\""<<buf<<"\"  index="<<index<<endl;
#endif

	//now loop to the end
	while (!fin.eof())
	{
		fin>>label;

		//now build the map of parameter <--> array
		//and put all these values into the array
		if(label=='A') 
		{
			fin>>ib>>ic>>id>>ie;
			for(int i=0;i<=MaxOrderA;i++) fin>>par[i];
			for(int i=0;i<=MaxOrderA;i++) fin>>par_fixedflag[i];

			for(int i=0;i<=MaxOrderA;i++) 
			{
				ia=i;
				mFitParaIndexMap[index][0]=ia;
				mFitParaIndexMap[index][1]=ib;
				mFitParaIndexMap[index][2]=ic;
				mFitParaIndexMap[index][3]=id;
				mFitParaIndexMap[index][4]=ie;
				mFitParaIndexMap[index][5]=label;
				pBufIndex = (((ia*(MaxOrderB+1) + ib)*(MaxOrderC+1) + ic)*(MaxOrderD+1) +id)*(MaxOrderE+1) + ie;
				mFitPara[index]=par[i];
				mFitParaFixedFlag[index]=par_fixedflag[i];
				mPara[ia][ib][ic][id][ie]=par[i];

#ifdef ApplyMulDPolN_DEBUG
				if(ApplyMulDPolN_DEBUG>=3)
					cout<<"BufIndex="<<setw(4)<<pBufIndex<<"   K_"<<ia<<ib<<ic<<id<<ie
					<<"="<<setw(14)<<mFitPara[index]<<"  fixedflag="<<mFitParaFixedFlag[index]<<endl;
#endif

				index++;
				//sanity check
				if(index>pMaxParaNum) 
				{
					cout<<"***Error! Input file specify too many parameters!***\n"<<endl;
					break;
				}
			}
		}
		else if(label=='B') 
		{
			fin>>ia>>ic>>id>>ie;
			for(int i=0;i<=MaxOrderB;i++) fin>>par[i];
			for(int i=0;i<=MaxOrderB;i++) fin>>par_fixedflag[i];

			for(int i=0;i<=MaxOrderB;i++) 
			{
				ib=i;
				mFitParaIndexMap[index][0]=ia;
				mFitParaIndexMap[index][1]=ib;
				mFitParaIndexMap[index][2]=ic;
				mFitParaIndexMap[index][3]=id;
				mFitParaIndexMap[index][4]=ie;
				mFitParaIndexMap[index][5]=label;
				pBufIndex = (((ia*(MaxOrderB+1) + ib)*(MaxOrderC+1) + ic)*(MaxOrderD+1) +id)*(MaxOrderE+1) + ie;
				mFitPara[index]=par[i];
				mFitParaFixedFlag[index]=par_fixedflag[i];
				mPara[ia][ib][ic][id][ie]=par[i];
#ifdef ApplyMulDPolN_DEBUG
				if(ApplyMulDPolN_DEBUG>=3)
					cout<<"BufIndex="<<setw(4)<<pBufIndex<<"   K_"<<ia<<ib<<ic<<id<<ie
					<<"="<<setw(14)<<mFitPara[index]<<"  fixedflag="<<mFitParaFixedFlag[index]<<endl;
#endif
				index++;
				//sanity check
				if(index>pMaxParaNum) 
				{
					cout<<"***Error! Input file specify too many parameters!***\n"<<endl;
					break;
				}
			}
		}
		else if(label=='C') 
		{
			fin>>ia>>ib>>id>>ie;
			for(int i=0;i<=MaxOrderC;i++) fin>>par[i];
			for(int i=0;i<=MaxOrderC;i++) fin>>par_fixedflag[i];

			for(int i=0;i<=MaxOrderC;i++) 
			{
				ic=i;
				mFitParaIndexMap[index][0]=ia;
				mFitParaIndexMap[index][1]=ib;
				mFitParaIndexMap[index][2]=ic;
				mFitParaIndexMap[index][3]=id;
				mFitParaIndexMap[index][4]=ie;
				mFitParaIndexMap[index][5]=label;
				pBufIndex = (((ia*(MaxOrderB+1) + ib)*(MaxOrderC+1) + ic)*(MaxOrderD+1) +id)*(MaxOrderE+1) + ie;
				mFitPara[index]=par[i];
				mFitParaFixedFlag[index]=par_fixedflag[i];
				mPara[ia][ib][ic][id][ie]=par[i];
#ifdef ApplyMulDPolN_DEBUG
				if(ApplyMulDPolN_DEBUG>=3)
					cout<<"BufIndex="<<setw(4)<<pBufIndex<<"   K_"<<ia<<ib<<ic<<id<<ie
					<<"="<<setw(14)<<mFitPara[index]<<"  fixedflag="<<mFitParaFixedFlag[index]<<endl;
#endif
				index++;
				//sanity check
				if(index>pMaxParaNum) 
				{
					cout<<"***Error! Input file specify too many parameters!***\n"<<endl;
					break;
				}
			}
		}
		else if(label=='D') 
		{
			fin>>ia>>ib>>ic>>ie;
			for(int i=0;i<=MaxOrderD;i++) fin>>par[i];
			for(int i=0;i<=MaxOrderD;i++) fin>>par_fixedflag[i];

			for(int i=0;i<=MaxOrderD;i++) 
			{
				id=i;
				mFitParaIndexMap[index][0]=ia;
				mFitParaIndexMap[index][1]=ib;
				mFitParaIndexMap[index][2]=ic;
				mFitParaIndexMap[index][3]=id;
				mFitParaIndexMap[index][4]=ie;
				mFitParaIndexMap[index][5]=label;
				pBufIndex = (((ia*(MaxOrderB+1) + ib)*(MaxOrderC+1) + ic)*(MaxOrderD+1) +id)*(MaxOrderE+1) + ie;
				mFitPara[index]=par[i];
				mFitParaFixedFlag[index]=par_fixedflag[i];
				mPara[ia][ib][ic][id][ie]=par[i];
#ifdef ApplyMulDPolN_DEBUG
				if(ApplyMulDPolN_DEBUG>=3)
					cout<<"BufIndex="<<setw(4)<<pBufIndex<<"   K_"<<ia<<ib<<ic<<id<<ie
					<<"="<<setw(14)<<mFitPara[index]<<"  fixedflag="<<mFitParaFixedFlag[index]<<endl;
#endif
				index++;
				//sanity check
				if(index>pMaxParaNum) 
				{
					cout<<"***Error! Input file specify too many parameters!***\n"<<endl;
					break;
				}
			}
		}
		else if(label=='E') 
		{
			fin>>ia>>ib>>ic>>id;
			for(int i=0;i<=MaxOrderE;i++) fin>>par[i];
			for(int i=0;i<=MaxOrderE;i++) fin>>par_fixedflag[i];

			for(int i=0;i<=MaxOrderE;i++) 
			{
				ie=i;
				mFitParaIndexMap[index][0]=ia;
				mFitParaIndexMap[index][1]=ib;
				mFitParaIndexMap[index][2]=ic;
				mFitParaIndexMap[index][3]=id;
				mFitParaIndexMap[index][4]=ie;
				mFitParaIndexMap[index][5]=label;
				pBufIndex = (((ia*(MaxOrderB+1) + ib)*(MaxOrderC+1) + ic)*(MaxOrderD+1) +id)*(MaxOrderE+1) + ie;
				mFitPara[index]=par[i];
				mFitParaFixedFlag[index]=par_fixedflag[i];
				mPara[ia][ib][ic][id][ie]=par[i];
#ifdef ApplyMulDPolN_DEBUG
				if(ApplyMulDPolN_DEBUG>=3)
					cout<<"BufIndex="<<setw(4)<<pBufIndex<<"   K_"<<ia<<ib<<ic<<id<<ie
					<<"="<<setw(14)<<mFitPara[index]<<"  fixedflag="<<mFitParaFixedFlag[index]<<endl;
#endif
				index++;
				//sanity check
				if(index>pMaxParaNum) 
				{
					cout<<"***Error! Input file specify too many parameters!***\n"<<endl;
					break;
				}
			}
		}
		else if(label=='x') 
		{
			//already end of file, or an empty line
			//do nothing
			//break;
		}
		else
		{
			cout<<"Ignoring unknown variable label '"<<label<<"' in input file, Line "<<line<<endl;
		}
		//eat the rest of the line		
		fin.getline(buf,511); line++; 
		//reset the label to 'x' so I can check if it is an empty line
		label='x';

#ifdef ApplyMulDPolN_DEBUG
		if(ApplyMulDPolN_DEBUG>=4) cout<<"Line="<<line<<"  buf=\""<<buf<<"\"  index="<<index<<endl;
#endif

	}

	//now build other buffers
	mParaNum=index;
#ifdef ApplyMulDPolN_DEBUG
	if(ApplyMulDPolN_DEBUG>=1) cout<<"Number of fit parameters = "<<mParaNum<<endl;
#endif
	fin.close();
	delete par;
	delete par_fixedflag;


	//just for for debug
	//DebugArray(mPara,MaxOrderA+1,MaxOrderB+1,MaxOrderC+1,MaxOrderD+1,MaxOrderE+1);
	//double x[5]={0,1,2,3,4};
	//myApplyMulDPolN(x,mFitPara); exit(0);
}


void  ApplyMulDPolN::SetMatrix(double *par)
{
#ifdef ApplyMulDPolN_DEBUG
	if(ApplyMulDPolN_DEBUG>=3)
		cout<<"ApplyMulDPolN::SetMatrix()"<<endl;
#endif

	if(!mPara) return;
	//now need to set parameters values to the array element
	for(int i=0;i<mParaNum;i++)
	{
		mFitPara[i]=par[i];
		int ia=mFitParaIndexMap[i][0];
		int ib=mFitParaIndexMap[i][1];
		int ic=mFitParaIndexMap[i][2];
		int id=mFitParaIndexMap[i][3];
		int ie=mFitParaIndexMap[i][4];
		mPara[ia][ib][ic][id][ie]=par[i];
#ifdef ApplyMulDPolN_DEBUG 
		if(ApplyMulDPolN_DEBUG>=5)
		{
			cout<<"ApplyMulDPolN::SetParameters() set K_"<<ia<<ib<<ic<<id<<ie<<"="<<mPara[ia][ib][ic][id][ie]<<endl;
		}
#endif
	}
	//Debug5DArray(); 
}



//static
double ApplyMulDPolN::MulDPolN(double *x, double *****matrix, 
							   int bin1, int bin2, int bin3, int bin4, int bin5)
{
	//multi-dimention polN
	//  V = Sum_abcde {C_abcde * X^a T^b Y^C F^d P^e}
	//    = Sum_bcde {Sum_a{C_abcde * X^a} * T^b Y^C F^d P^e}

	double var=0;
	for(int ia=0;ia<=bin1;ia++)
	{
		double itemA = pow(x[0],double(ia));
		for(int ib=0;ib<=bin2;ib++)
		{
			double itemB = pow(x[1],double(ib));
			for(int ic=0;ic<=bin3;ic++)
			{
				double itemC = pow(x[2],double(ic));
				for(int id=0;id<=bin4;id++)
				{
					double itemD = pow(x[3],double(id));
					for(int ie=0;ie<=bin5;ie++)
					{
						double itemE = pow(x[4],double(ie));
#ifdef ApplyMulDPolN_DEBUG 
						if(ApplyMulDPolN_DEBUG>=6)
						{
							cout<<"ApplyMulDPolN::ApplyMulDPolN():  K_"<<ia<<ib<<ic<<id<<ie<<"="<<matrix[ia][ib][ic][id][ie]<<endl; 
						}
#endif
						var += matrix[ia][ib][ic][id][ie] * itemA * itemB *itemC *itemD *itemE;
					}
				}
			}
		}
	}
	return var;
}

//static
void   ApplyMulDPolN::SetMatrix(double *par, int parnum, double **parmap, double *****matrix, 
							   int bin1, int bin2, int bin3, int bin4, int bin5)
{
	//multi-dimention polN
	//  V = Sum_abcde {C_abcde * X^a T^b Y^C F^d P^e}
	//    = Sum_bcde {Sum_a{C_abcde * X^a} * T^b Y^C F^d P^e}

	if(!matrix) return;
	//now need to set parameters values to the array element
	for(int i=0;i<parnum;i++)
	{
		int ia=parmap[i][0];
		int ib=parmap[i][1];
		int ic=parmap[i][2];
		int id=parmap[i][3];
		int ie=parmap[i][4];
		matrix[ia][ib][ic][id][ie]=par[i];
#ifdef ApplyMulDPolN_DEBUG 
		if(ApplyMulDPolN_DEBUG>=5)
		{
			cout<<"ApplyMulDPolN::SetParameters() set K_"<<ia<<ib<<ic<<id<<ie<<"="<<matrix[ia][ib][ic][id][ie]<<endl;
		}
#endif
	}
	//Debug5DArray(); 
}



void   ApplyMulDPolN::SetMatrix(double *****matrix, 
							   int bin1, int bin2, int bin3, int bin4, int bin5)
{
	for(int ia=0;ia<=bin1;ia++)
	{
		for(int ib=0;ib<=bin2;ib++)
		{
			for(int ic=0;ic<=bin3;ic++)
			{
				for(int id=0;id<=bin4;id++)
				{
					for(int ie=0;ie<=bin5;ie++)
					{
						mPara[ia][ib][ic][id][ie] = matrix[ia][ib][ic][id][ie];
					}
				}
			}
		}
	}
}


double***** ApplyMulDPolN::GetMatrix(int &bin1, int &bin2, int &bin3, int &bin4, int &bin5) 
{
	bin1 = MaxOrderA;
	bin2 = MaxOrderB;
	bin3 = MaxOrderC;
	bin4 = MaxOrderD;
	bin5 = MaxOrderE;
	return mPara;
}


double ApplyMulDPolN::MulDPolN(double *x)
{
#ifdef ApplyMulDPolN_DEBUG
	if(ApplyMulDPolN_DEBUG>=4)
	{
		cout<<"ApplyMulDPolN::MulDPolN()\n";
	}
#endif
	//multi-dimention polN
	//  V = Sum_abcde {C_abcde * X^a T^b Y^C F^d P^e}
	//    = Sum_bcde {Sum_a{C_abcde * X^a} * T^b Y^C F^d P^e}

	double var=0;
	for(int ia=0;ia<=MaxOrderA;ia++)
	{
		double itemA = pow(x[0],double(ia));
		for(int ib=0;ib<=MaxOrderB;ib++)
		{
			double itemB = pow(x[1],double(ib));
			for(int ic=0;ic<=MaxOrderC;ic++)
			{
				double itemC = pow(x[2],double(ic));
				for(int id=0;id<=MaxOrderD;id++)
				{
					double itemD = pow(x[3],double(id));
					for(int ie=0;ie<=MaxOrderE;ie++)
					{
						double itemE = pow(x[4],double(ie));
#ifdef ApplyMulDPolN_DEBUG 
						if(ApplyMulDPolN_DEBUG>=6)
						{
							cout<<"ApplyMulDPolN::MulDPolN():  K_"<<ia<<ib<<ic<<id<<ie<<"="<<mPara[ia][ib][ic][id][ie]<<endl; 
						}
#endif
						var += mPara[ia][ib][ic][id][ie] * itemA * itemB *itemC *itemD *itemE;
					}
				}
			}
		}
	}
	return var;
}



double ApplyMulDPolN::Eval(double *x)
{
#ifdef ApplyMulDPolN_DEBUG
	if(ApplyMulDPolN_DEBUG>=4)
	{
		cout<<"ApplyMulDPolN::Eval()\n";
	}
#endif

	return MulDPolN(x);
}

