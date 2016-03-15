//This file contains lot of template to allocate and free N-Dimention arrays, where N=2-5.
//And also template function to print those array into a INC file of ascii file.
//Please note that template functions start with create but normal functions start with print.
#ifndef _ACCEPTANCE_TOOLS_
#define _ACCEPTANCE_TOOLS_

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <typeinfo>
#include <iomanip>
#include <iostream>
#include <fstream>
using namespace std;

//template PrintBuffer will print all types arrays into ascii files or terminals
//template CreateFile will print ascii files using all types of pointer
//template CreatIncFile will print include files using all types of pointer
//template CalculateRatio can calculate the ratio with all types of numerator and denominator
//template PrintACCTable can print ascii files with various types  of numerator and denominator and various levels

namespace ACCEPTANCE {
	/////////////////////////////////////////////////////////////////////////////

	int BinarySearch(const int    sortedArray[], int first, int last, const int    key);
	int BinarySearch(const double sortedArray[], int first, int last, const double key);

	void CreateQ2Table(double Q2Table[]);
	int  GetQ2Index(double Q2);

	/////////////////////////////////////////////////////////////////////////////
	//these function will print C type inc file, but support double type only
	//input: val[bin1][bin2][bin3][bin4]
	//input: char *name[] is string of "table_name, bin1name bin2name ...bin#name "
	//example: char *name[]={"Table_4D","Bin1","Bin2","Bin3","Bin4"}
	void PrintIncFile(double ****val, int bin1, int bin2, int bin3,int bin4, char **name=0, int ncol=8);
	void PrintIncFile(double ***val, int bin1, int bin2, int bin3, char **name=0, int ncol=8);
	void PrintIncFile(double **val, int bin1, int bin2, char **name=0, int ncol=8);	

	//these functions will print ascii table file, but support double type only
	void PrintFile(double *val, int bin1, char **name, int ncol=8, bool bPrintIndex=false);
	void PrintFile(double **val, int bin1, int bin2, char **name, int ncol=8, bool bPrintIndex=false);
	void PrintFile(double ***val, int bin1, int bin2, int bin3, char **name, int ncol=8, bool bPrintIndex=false);
	void PrintFile(double ****val, int bin1, int bin2, int bin3,int bin4, char **name, int ncol=8, bool bPrintIndex=false);

	//this function can print either inc or table file
	void PrintTable(double ****val, int bin1, int bin2, int bin3,int bin4, char **name, int ncol=8, 
		bool bPrintIncFile=true,bool bPrintIndex=true);

	//these functions will calculate the ratio using 2 int pointers
	double **** CalculateACC(int ****N_det, int ****N_true, int bin1, int bin2,int bin3,int bin4);
	double *** CalculateACC(int ***N_det, int ***N_true, int bin1, int bin2,int bin3);
	double ** CalculateACC(int **N_det, int **N_true, int bin1, int bin2);    
	double * CalculateACC(int *N_det, int *N_true, int bin1);

	void RemoveAll(std::string str, std::string key);

	/////////////////////////////////////////////////////////////////////////////

	//input: 5D array or 5D pointer
	template <class T_5DArr,typename T>
	T ***** CopyToTypeT(int bin1, int bin2, int bin3, int bin4,int bin5, T_5DArr defaultarr, T targettype)
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
						for(int ie=0;ie<bin5;ie++)  
						{
							//int index=iw*bin2*bin3*bin4*bin5+iq*bin3*bin4*bin5+it*bin4*bin5+ip*bin5+ie;
							//arr[iw][iq][it][ip][ie]=defaultarr[index];
							arr[iw][iq][it][ip][ie]=(T)defaultarr[iw][iq][it][ip][ie];
						}
					}
				}
			}
		}
		return arr;
	};
	/////////////////////////////////////////////////////////////////////////////

	//input: 4D array or 4D pointer
	template <class T_4DArr, typename T>
	T **** CopyToTypeT(int bin1, int bin2, int bin3, int bin4, T_4DArr defaultarr, T targettype)
	{
		T **** arr;
		arr=new T *** [bin1];
		for(int iw=0;iw<bin1;iw++)
		{  	  
			arr[iw]=new T ** [bin2];			 
			for(int iq=0;iq<bin2;iq++)
			{  
				arr[iw][iq]=new T * [bin3];					 
				for(int it=0;it<bin3;it++)
				{	
					arr[iw][iq][it]=new T [bin4];                    
					for(int ip=0;ip<bin4;ip++)  
					{                        
						//int index=iw*bin2*bin3*bin4+iq*bin3*bin4+it*bin4+ip;    
						//arr[iw][iq][it][ip]=defaultarr[index];
						arr[iw][iq][it][ip]=(T)defaultarr[iw][iq][it][ip];
					}
				}
			}
		}
		return arr;
	};

	/////////////////////////////////////////////////////////////////////////////

	//input: 3D array or 3D pointer
	template <class T_3DArr, typename T>
	T  *** CopyToTypeT( int bin1, int bin2, int bin3, T_3DArr defaultarr, T targettype)
	{
		T ***arr;
		arr=new T ** [bin1];
		for(int iw=0;iw<bin1;iw++)
		{  	  
			arr[iw]=new T * [bin2];			 
			for(int iq=0;iq<bin2;iq++)
			{  
				arr[iw][iq]=new T [bin3];				 
				for(int it=0;it<bin3;it++) 
				{
					//int index=iw*bin2*bin3+iq*bin3+it;	
					//arr[iw][iq][it]=defaultarr[index];	
					arr[iw][iq][it]=(T)defaultarr[iw][iq][it];
				}
			}
		}
		return arr;
	};

	/////////////////////////////////////////////////////////////////////////////

	//input: 2D array or 2D pointer
	template <class T_2DArr, typename T>
	T ** CopyToTypeT( int bin1, int bin2, T_2DArr defaultarr, T targettype)
	{
		T **arr;
		arr=new T * [bin1];
		for(int iw=0;iw<bin1;iw++)
		{  	  
			arr[iw]=new T [bin2];			 
			for(int iq=0;iq<bin2;iq++) 
			{                
				//int index=iw*bin2+iq;	
				//arr[iw][iq]=defaultarr[index];	
				arr[iw][iq]=(T)defaultarr[iw][iq];	
			}
		}
		return arr;
	};

	/////////////////////////////////////////////////////////////////////////////

	//input: 1D array or 1D pointer
	template <class T_1DArr, typename T>
	T * CopyToTypeT( int bin1, T_1DArr defaultarr, T targettype)
	{
		T *arr;
		arr=new T [bin1];
		for(int iw=0;iw<bin1;iw++) arr[iw]=(T)defaultarr[iw];
		return arr;
	};


	//////////////////////////////////////////////////////////////////////////
	//input: pointer to the buffer
	template <typename T>
	T ***** CloneArray(int bin1, int bin2, int bin3, int bin4,int bin5, T *pointer)
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
						for(int ie=0;ie<bin5;ie++)  
						{
							int index=iw*bin2*bin3*bin4*bin5+iq*bin3*bin4*bin5+it*bin4*bin5+ip*bin5+ie;
							arr[iw][iq][it][ip][ie]=pointer[index];
						}
					}
				}
			}
		}
		return arr;
	};
	/////////////////////////////////////////////////////////////////////////////
	
	template <typename T>
	T **** CloneArray(int bin1, int bin2, int bin3, int bin4, T *pointer)
	{
		T **** arr;
		arr=new T *** [bin1];
		for(int iw=0;iw<bin1;iw++)
		{  	  
			arr[iw]=new T ** [bin2];			 
			for(int iq=0;iq<bin2;iq++)
			{  
				arr[iw][iq]=new T * [bin3];					 
				for(int it=0;it<bin3;it++)
				{	
					arr[iw][iq][it]=new T [bin4];                    
					for(int ip=0;ip<bin4;ip++)  
					{                        
						int index=iw*bin2*bin3*bin4+iq*bin3*bin4+it*bin4+ip;    
						arr[iw][iq][it][ip]=pointer[index];
					}
				}
			}
		}
		return arr;
	};

	/////////////////////////////////////////////////////////////////////////////

	template <typename T>
	T *** CloneArray( int bin1, int bin2, int bin3, T *pointer)
	{
		T ***arr;
		arr=new T ** [bin1];
		for(int iw=0;iw<bin1;iw++)
		{  	  
			arr[iw]=new T * [bin2];			 
			for(int iq=0;iq<bin2;iq++)
			{  
				arr[iw][iq]=new T [bin3];				 
				for(int it=0;it<bin3;it++) 
				{
					int index=iw*bin2*bin3+iq*bin3+it;	
					arr[iw][iq][it]=pointer[index];
				}
			}
		}
		return arr;
	};

	/////////////////////////////////////////////////////////////////////////////

	template <typename T>
	T ** CloneArray( int bin1, int bin2, T *pointer)
	{
		T **arr;
		arr=new T * [bin1];
		for(int iw=0;iw<bin1;iw++)
		{  	  
			arr[iw]=new T [bin2];			 
			for(int iq=0;iq<bin2;iq++) 
			{                
				int index=iw*bin2+iq;	
				arr[iw][iq]=pointer[index];	
			}
		}
		return arr;
	};

	/////////////////////////////////////////////////////////////////////////////

	template <typename T>
	T * CloneArray( int bin1, T *pointer)
	{
		T *arr;
		arr=new T [bin1];
		for(int iw=0;iw<bin1;iw++) arr[iw]=pointer[iw];
		return arr;
	};


	//////////////////////////////////////////////////////////////////////////
	template <typename T>
	void FreeDynArray(T *****arr, int bin1, int bin2, int bin3,int bin4)
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
	T ***** IniDynamicArray(int bin1, int bin2, int bin3, int bin4,int bin5, T defautval)
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
	/////////////////////////////////////////////////////////////////////////////
	template <typename T>
	void FreeDynArray(T ****arr, int bin1, int bin2, int bin3)
	{   
		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)
			{  
				for(int it=0;it<bin3;it++)  delete [] arr[iw][iq][it];
			}
		}
	};

	//create a buffer using a given value
	template <typename T>
	T **** IniDynamicArray(int bin1, int bin2, int bin3, int bin4, T defautval)
	{
		T **** arr;
		arr=new T *** [bin1];
		for(int iw=0;iw<bin1;iw++)
		{  	  
			arr[iw]=new T ** [bin2];			 
			for(int iq=0;iq<bin2;iq++)
			{  
				arr[iw][iq]=new T * [bin3];					 
				for(int it=0;it<bin3;it++)
				{	
					arr[iw][iq][it]=new T [bin4];
					for(int ip=0;ip<bin4;ip++)  arr[iw][iq][it][ip]=defautval;
				}
			}
		}
		//printf("IniDynamicArray() address: %p\n",arr);
		return arr;
	};

	/////////////////////////////////////////////////////////////////////////////
	template <typename T>
	void FreeDynArray(T ***arr, int bin1, int bin2)
	{   
		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)  delete [] arr[iw][iq];		
		}
	};

	//create a buffer using a given value
	template <typename T>
	T *** IniDynamicArray( int bin1, int bin2, int bin3, T defautval)
	{
		T ***arr;
		arr=new T ** [bin1];
		for(int iw=0;iw<bin1;iw++)
		{  	  
			arr[iw]=new T * [bin2];			 
			for(int iq=0;iq<bin2;iq++)
			{  
				arr[iw][iq]=new T [bin3];					 
				for(int it=0;it<bin3;it++) arr[iw][iq][it]=defautval;	 			
			}
		}
		return arr;
	};


	/////////////////////////////////////////////////////////////////////////////

	template <typename T>
	void FreeDynArray(T **arr, int bin1)
	{   
		for(int iw=0;iw<bin1;iw++) delete [] arr[iw];	
	};

	//create a buffer using a given value
	template <typename T>
	T ** IniDynamicArray( int bin1, int bin2, T defautval)
	{
		T **arr;
		arr=new T * [bin1];
		for(int iw=0;iw<bin1;iw++)
		{  	  
			arr[iw]=new T [bin2];			 
			for(int iq=0;iq<bin2;iq++) arr[iw][iq]=defautval;	
		}
		return arr;
	};

	/////////////////////////////////////////////////////////////////////////////

	template <typename T>
	void FreeDynArray(T *arr)
	{   
		delete [] arr;	
	};

	//create a buffer using a given value
	template <typename T>
	T * IniDynamicArray( int bin1, T defautval)
	{
		T *arr;
		arr=new T [bin1];
		for(int iw=0;iw<bin1;iw++) arr[iw]=defautval;
		return arr;
	};

//////////////////////////////////////////////////////////////////////////////
		//print include file with 5-D Array or 5-D Array pointer
	template <class T_5DArr>
	void CreateIncFile(T_5DArr val, int bin1, int bin2, int bin3,int bin4,int bin5, char **name=0, int ncol=8, int precision=6)
	{
		//write out this table
		char *DefaultName[]={"Table_4D","Bin1","Bin2","Bin3","Bin4","Bin5"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(val[0][0][0][0][0])==typeid(int) || typeid(val[0][0][0][0][0])==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(val[0][0][0][0][0])==typeid(float) ||  typeid(val[0][0][0][0][0])==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 5-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"]["<<name[3]<<"]["<<name[4]<<"]["<<name[5]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";
		facc<<"//static const int "<<name[3]<<"="<<bin3<<"; // bin3 \n";
		facc<<"//static const int "<<name[4]<<"="<<bin4<<"; // bin4 \n";
		facc<<"//static const int "<<name[5]<<"="<<bin5<<"; // bin5 \n";

		int NCol = (bin5>ncol)?ncol:bin5;


		facc<<"static " <<GetTypeName(val[0][0][0][0][0])<<"  "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]["<<bin4<<"]["<<bin5<<"]={"<<std::endl;
		for(int iw=0;iw<bin1;iw++)
		{  
			facc<<"\t{ //start of bin1 "<<name[1]<<"="<<iw<<std::endl;
			for(int iq=0;iq<bin2;iq++)
			{  
				facc<<"\t\t{ //start of bin2 "<<name[2]<<"="<<iq<<std::endl;
				for(int it=0;it<bin3;it++)
				{
					facc<<"\t\t\t{ //start of bin3 "<<name[3]<<"="<<it<<std::endl;
					for(int ip=0;ip<bin4;ip++)
					{
						facc<<"\t\t\t\t{ ";
						if(bin5/NCol>10) facc<<" //start of bin4 "<<name[4]<<"="<<ip;
						facc<<std::endl;
						for(int i5=0;i5<bin5;i5++)
						{
							if (!(i5%NCol)) facc<<"\t\t\t\t"; 
							facc<<std::setw(span)<<val[iw][iq][it][ip][i5];
							if(i5<bin5-1) facc<<", ";
							if ( !((i5+1)%NCol) || i5==bin5-1) facc<<"\n"; 
						}	

						if(ip<bin4-1) facc<<"\t\t\t\t},";
						else facc<<"\t\t\t\t}";
						if(bin5/NCol>10) facc<<" //  end of bin4 "<<name[4]<<"="<<ip;
						facc<<std::endl;
					}	

					if(it<bin3-1) facc<<"\t\t\t},";
					else facc<<"\t\t\t} //  end of bin3 "<<name[3]<<"="<<it<<std::endl;
				}  
				if(iq<bin2-1) facc<<"\t\t}, //  end of bin2 "<<name[2]<<"="<<iq<<std::endl;
				else facc<<"\t\t} //  end of bin2 "<<name[2]<<"="<<iq<<std::endl;
			} 
			if(iw<bin1-1) facc<<"\t}, //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
			else facc<<"\t} //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
		}
		facc<<"};// end of const double "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]["<<bin4<<"]["<<bin5<<"]"<<std::endl;
		facc.close();
	}

	/////////////////////////////////////////////////////////////////////////////    
	//print include file with 4-D Array or 4-D array pointer
	template <class T_4DArr>
	void CreateIncFile(T_4DArr val, int bin1, int bin2, int bin3,int bin4, char **name=0, int ncol=8,int precision=6)
	{
		//write out this table
		char *DefaultName[]={"Table_4D","Bin1","Bin2","Bin3","Bin4"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(val[0][0][0][0])==typeid(int) || typeid(val[0][0][0][0])==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(val[0][0][0][0])==typeid(float) ||  typeid(val[0][0][0][0])==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 4-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"]["<<name[3]<<"]["<<name[4]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";
		facc<<"//static const int "<<name[3]<<"="<<bin3<<"; // bin3 \n";
		facc<<"//static const int "<<name[4]<<"="<<bin4<<"; // bin4 \n";

		int NCol = (bin4>ncol)?ncol:bin4;

		facc<<"static " <<GetTypeName(val[0][0][0][0])<<"  "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]["<<bin4<<"]={"<<std::endl;
		for(int iw=0;iw<bin1;iw++)
		{  
			facc<<"\t{ //start of bin1 "<<name[1]<<"="<<iw<<std::endl;
			for(int iq=0;iq<bin2;iq++)
			{  
				facc<<"\t\t{ //start of bin2 "<<name[2]<<"="<<iq<<std::endl;
				for(int it=0;it<bin3;it++)
				{
					facc<<"\t\t\t{ ";
					if(bin4/NCol>10) facc<<" //start of bin3 "<<name[3]<<"="<<it;
					facc<<std::endl;
					for(int ip=0;ip<bin4;ip++)
					{
						if (!(ip%NCol)) facc<<"\t\t\t\t"; 
						facc<<std::setw(span)<<val[iw][iq][it][ip];
						if(ip<bin4-1) facc<<", ";
						if ( !((ip+1)%NCol) || ip==bin4-1) facc<<"\n"; 
					}	

					if(it<bin3-1) facc<<"\t\t\t},";
					else facc<<"\t\t\t}";
					if(bin4/NCol>10) facc<<" //  end of bin3 "<<name[3]<<"="<<it;
					facc<<std::endl;
				}  
				if(iq<bin2-1) facc<<"\t\t}, //  end of bin2 "<<name[2]<<"="<<iq<<std::endl;
				else facc<<"\t\t} //  end of bin2 "<<name[2]<<"="<<iq<<std::endl;
			} 
			if(iw<bin1-1) facc<<"\t}, //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
			else facc<<"\t} //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
		}
		facc<<"};// end of const double "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]["<<bin4<<"]"<<std::endl;
		facc.close();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////    
	//print include file with 3-D Array or 3-D Array pointer
	template <class T_3DArr>
	void CreateIncFile(T_3DArr val, int bin1, int bin2, int bin3, char **name=0, int ncol=8, int precision=6)
	{
		//write out this table
		char *DefaultName[]={"Table_3D","Bin1","Bin2","Bin3"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(val[0][0][0])==typeid(int) || typeid(val[0][0][0])==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(val[0][0][0])==typeid(float) ||  typeid(val[0][0][0])==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 3-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"]["<<name[3]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";
		facc<<"//static const int "<<name[3]<<"="<<bin3<<"; // bin3 \n";

		int NCol = (bin3>ncol)?ncol:bin3;
		facc<<"static " <<GetTypeName(val[0][0][0])<<"  "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]={"<<std::endl;
		for(int iw=0;iw<bin1;iw++)
		{  
			facc<<"\t{ //start of bin1 "<<name[1]<<"="<<iw<<std::endl;
			for(int iq=0;iq<bin2;iq++)
			{  
				facc<<"\t\t{ ";
				if(bin3/NCol>10) facc<<" //start of bin2 "<<name[2]<<"="<<iq;
				facc<<std::endl;
				for(int it=0;it<bin3;it++)
				{			 
					if (!(it%NCol)) facc<<"\t\t\t"; 
					facc<<std::setw(span)<<val[iw][iq][it];
					if(it<bin3-1) facc<<", ";
					if ( !((it+1)%NCol) || it==bin3-1) facc<<"\n"; 				
				}  
				if(iq<bin2-1) facc<<"\t\t},";
				else facc<<"\t\t}";
				if(bin3/NCol>10) facc<<" //  end of bin2 "<<name[2]<<"="<<iq;
				facc<<std::endl;
			} 
			if(iw<bin1-1) facc<<"\t}, //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
			else facc<<"\t} //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
		}
		facc<<"};// end of const double "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]"<<std::endl;
		facc.close();
	}

	/////////////////////////////////////////////////////////////////////////////////////////////    
	//print include file with 2-D Array or 2-D Array pointer
	template <class T_2DArr>
	void CreateIncFile(T_2DArr val, int bin1, int bin2, char **name=0, int ncol=8, int precision=6)
	{  
		char *DefaultName[]={"Table_2D","Bin1","Bin2"};
		if(name==0) name=DefaultName;
		//write out this table
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(val[0][0])==typeid(int) || typeid(val[0][0])==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(val[0][0])==typeid(float) ||  typeid(val[0][0])==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 2-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";


		int NCol = (bin2>ncol)?ncol:bin2;
		facc<<"static " <<GetTypeName(val[0][0])<<"  "<<table<<"["<<bin1<<"]["<<bin2<<"]={"<<std::endl;
		for(int iw=0;iw<bin1;iw++)
		{  
			facc<<"\t{ //start of bin1 "<<name[1]<<"="<<iw<<std::endl;
			for(int iq=0;iq<bin2;iq++)
			{  
				if (!(iq%NCol)) facc<<"\t\t"; 
				facc<<std::setw(span)<<val[iw][iq];
				if(iq<bin2-1) facc<<", ";
				if ( !((iq+1)%NCol) || iq==bin2-1) facc<<"\n"; 				
			} 
			if(iw<bin1-1) facc<<"\t}, //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
			else facc<<"\t} //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
		}
		facc<<"};// end of const double "<<table<<"["<<bin1<<"]["<<bin2<<"]"<<std::endl;
		facc.close();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////      
	//print include file with 1-D Array or a pointer
	template <class T_1DArr>
	void CreateIncFile(T_1DArr val, int bin1, char **name=0, int ncol=8, int precision=6)
	{  
		char *DefaultName[]={"Table_1D","Bin1"};
		if(name==0) name=DefaultName;
		//write out this table
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(val[0])==typeid(int) || typeid(val[0])==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(val[0])==typeid(float) ||  typeid(val[0])==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);

		facc<<"//This is the table for 1-d array "<<table<<"["<<name[1]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";


		int NCol = (bin1>ncol)?ncol:bin1;
		facc<<"static " <<GetTypeName(val[0])<<"  "<<table<<"["<<bin1<<"]={"<<std::endl;
		for(int iw=0;iw<bin1;iw++)
		{  
			if (!(iw%NCol)) facc<<"\t"; 
			facc<<std::setw(span)<<val[iw];
			if(iw<bin1-1) facc<<", ";
			if ( !((iw+1)%NCol) || iw==bin1-1) facc<<"\n"; 				

		}
		facc<<"};// end of const double "<<table<<"["<<bin1<<"]"<<std::endl;
		facc.close();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////    
	//print include file with 5-D Array pointer
	template <typename T>
	void CreateIncFile_ArrPointer(T *****val, int bin1, int bin2, int bin3,int bin4,int bin5, char **name=0, int ncol=8, int precision=6)
	{
		//write out this table
		char *DefaultName[]={"Table_4D","Bin1","Bin2","Bin3","Bin4","Bin5"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(T)==typeid(int) || typeid(T)==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(T)==typeid(float) ||  typeid(T)==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 5-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"]["<<name[3]<<"]["<<name[4]<<"]["<<name[5]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";
		facc<<"//static const int "<<name[3]<<"="<<bin3<<"; // bin3 \n";
		facc<<"//static const int "<<name[4]<<"="<<bin4<<"; // bin4 \n";
		facc<<"//static const int "<<name[5]<<"="<<bin5<<"; // bin5 \n";

		int NCol = (bin5>ncol)?ncol:bin5;


		facc<<"static " <<GetTypeName((T) 0.0)<<"  "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]["<<bin4<<"]["<<bin5<<"]={"<<std::endl;
		for(int iw=0;iw<bin1;iw++)
		{  
			facc<<"\t{ //start of bin1 "<<name[1]<<"="<<iw<<std::endl;
			for(int iq=0;iq<bin2;iq++)
			{  
				facc<<"\t\t{ //start of bin2 "<<name[2]<<"="<<iq<<std::endl;
				for(int it=0;it<bin3;it++)
				{
					facc<<"\t\t\t{ //start of bin3 "<<name[3]<<"="<<it<<std::endl;
					for(int ip=0;ip<bin4;ip++)
					{
						facc<<"\t\t\t\t{ ";
						if(bin5/NCol>10) facc<<" //start of bin4 "<<name[4]<<"="<<ip;
						facc<<std::endl;
						for(int i5=0;i5<bin5;i5++)
						{
							if (!(i5%NCol)) facc<<"\t\t\t\t"; 
							facc<<std::setw(span)<<val[iw][iq][it][ip][i5];
							if(i5<bin5-1) facc<<", ";
							if ( !((i5+1)%NCol) || i5==bin5-1) facc<<"\n"; 
						}	

						if(ip<bin4-1) facc<<"\t\t\t\t},";
						else facc<<"\t\t\t\t}";
						if(bin5/NCol>10) facc<<" //  end of bin4 "<<name[4]<<"="<<ip;
						facc<<std::endl;
					}	

					if(it<bin3-1) facc<<"\t\t\t},";
					else facc<<"\t\t\t} //  end of bin3 "<<name[3]<<"="<<it<<std::endl;
				}  
				if(iq<bin2-1) facc<<"\t\t}, //  end of bin2 "<<name[2]<<"="<<iq<<std::endl;
				else facc<<"\t\t} //  end of bin2 "<<name[2]<<"="<<iq<<std::endl;
			} 
			if(iw<bin1-1) facc<<"\t}, //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
			else facc<<"\t} //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
		}
		facc<<"};// end of const double "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]["<<bin4<<"]["<<bin5<<"]"<<std::endl;
		facc.close();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////    
	template <typename T>
	void CreateIncFile_ArrPointer(T ****val, int bin1, int bin2, int bin3,int bin4, char **name=0, int ncol=8,int precision=6)
	{
		//write out this table
		char *DefaultName[]={"Table_4D","Bin1","Bin2","Bin3","Bin4"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(T)==typeid(int) || typeid(T)==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(T)==typeid(float) ||  typeid(T)==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 4-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"]["<<name[3]<<"]["<<name[4]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";
		facc<<"//static const int "<<name[3]<<"="<<bin3<<"; // bin3 \n";
		facc<<"//static const int "<<name[4]<<"="<<bin4<<"; // bin4 \n";

		int NCol = (bin4>ncol)?ncol:bin4;

		facc<<"static " <<GetTypeName((T) 0.0)<<"  "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]["<<bin4<<"]={"<<std::endl;
		for(int iw=0;iw<bin1;iw++)
		{  
			facc<<"\t{ //start of bin1 "<<name[1]<<"="<<iw<<std::endl;
			for(int iq=0;iq<bin2;iq++)
			{  
				facc<<"\t\t{ //start of bin2 "<<name[2]<<"="<<iq<<std::endl;
				for(int it=0;it<bin3;it++)
				{
					facc<<"\t\t\t{ ";
					if(bin4/NCol>10) facc<<" //start of bin3 "<<name[3]<<"="<<it;
					facc<<std::endl;
					for(int ip=0;ip<bin4;ip++)
					{
						if (!(ip%NCol)) facc<<"\t\t\t\t"; 
						facc<<std::setw(span)<<val[iw][iq][it][ip];
						if(ip<bin4-1) facc<<", ";
						if ( !((ip+1)%NCol) || ip==bin4-1) facc<<"\n"; 
					}	

					if(it<bin3-1) facc<<"\t\t\t},";
					else facc<<"\t\t\t}";
					if(bin4/NCol>10) facc<<" //  end of bin3 "<<name[3]<<"="<<it;
					facc<<std::endl;
				}  
				if(iq<bin2-1) facc<<"\t\t}, //  end of bin2 "<<name[2]<<"="<<iq<<std::endl;
				else facc<<"\t\t} //  end of bin2 "<<name[2]<<"="<<iq<<std::endl;
			} 
			if(iw<bin1-1) facc<<"\t}, //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
			else facc<<"\t} //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
		}
		facc<<"};// end of const double "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]["<<bin4<<"]"<<std::endl;
		facc.close();
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////    
	//print include file with 3-D Array pointer
	template <typename T>
	void CreateIncFile_ArrPointer(T ***val, int bin1, int bin2, int bin3, char **name=0, int ncol=8, int precision=6)
	{
		//write out this table
		char *DefaultName[]={"Table_3D","Bin1","Bin2","Bin3"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(T)==typeid(int) || typeid(T)==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(T)==typeid(float) ||  typeid(T)==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 3-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"]["<<name[3]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";
		facc<<"//static const int "<<name[3]<<"="<<bin3<<"; // bin3 \n";

		int NCol = (bin3>ncol)?ncol:bin3;
		facc<<"static " <<GetTypeName((T) 0.0)<<"  "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]={"<<std::endl;
		for(int iw=0;iw<bin1;iw++)
		{  
			facc<<"\t{ //start of bin1 "<<name[1]<<"="<<iw<<std::endl;
			for(int iq=0;iq<bin2;iq++)
			{  
				facc<<"\t\t{ ";
				if(bin3/NCol>10) facc<<" //start of bin2 "<<name[2]<<"="<<iq;
				facc<<std::endl;
				for(int it=0;it<bin3;it++)
				{			 
					if (!(it%NCol)) facc<<"\t\t\t"; 
					facc<<std::setw(span)<<val[iw][iq][it];
					if(it<bin3-1) facc<<", ";
					if ( !((it+1)%NCol) || it==bin3-1) facc<<"\n"; 				
				}  
				if(iq<bin2-1) facc<<"\t\t},";
				else facc<<"\t\t}";
				if(bin3/NCol>10) facc<<" //  end of bin2 "<<name[2]<<"="<<iq;
				facc<<std::endl;
			} 
			if(iw<bin1-1) facc<<"\t}, //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
			else facc<<"\t} //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
		}
		facc<<"};// end of const double "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]"<<std::endl;
		facc.close();
	}

	/////////////////////////////////////////////////////////////////////////////////////////////    
	//print include file with 2-D Array pointer
	template <typename T>
	void CreateIncFile_ArrPointer(T **val, int bin1, int bin2, char **name=0, int ncol=8, int precision=6)
	{  
		char *DefaultName[]={"Table_2D","Bin1","Bin2"};
		if(name==0) name=DefaultName;
		//write out this table
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(T)==typeid(int) || typeid(T)==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(T)==typeid(float) ||  typeid(T)==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 2-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";


		int NCol = (bin2>ncol)?ncol:bin2;
		facc<<"static " <<GetTypeName((T) 0.0)<<"  "<<table<<"["<<bin1<<"]["<<bin2<<"]={"<<std::endl;
		for(int iw=0;iw<bin1;iw++)
		{  
			facc<<"\t{ //start of bin1 "<<name[1]<<"="<<iw<<std::endl;
			for(int iq=0;iq<bin2;iq++)
			{  
				if (!(iq%NCol)) facc<<"\t\t"; 
				facc<<std::setw(span)<<val[iw][iq];
				if(iq<bin2-1) facc<<", ";
				if ( !((iq+1)%NCol) || iq==bin2-1) facc<<"\n"; 				
			} 
			if(iw<bin1-1) facc<<"\t}, //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
			else facc<<"\t} //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
		}
		facc<<"};// end of const double "<<table<<"["<<bin1<<"]["<<bin2<<"]"<<std::endl;
		facc.close();
	}

	//////////////////////////////////////////////////////////////////////////////////////////////      
	//print include file with 1-D Array pointer
	template <typename T>
	void CreateIncFile_ArrPointer(T *val, int bin1, char **name=0, int ncol=8, int precision=6)
	{  
		char *DefaultName[]={"Table_1D","Bin1"};
		if(name==0) name=DefaultName;
		//write out this table
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(T)==typeid(int) || typeid(T)==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(T)==typeid(float) ||  typeid(T)==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);

		facc<<"//This is the table for 1-d array "<<table<<"["<<name[1]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";


		int NCol = (bin1>ncol)?ncol:bin1;
		facc<<"static " <<GetTypeName((T) 0.0)<<"  "<<table<<"["<<bin1<<"]={"<<std::endl;
		for(int iw=0;iw<bin1;iw++)
		{  
			if (!(iw%NCol)) facc<<"\t"; 
			facc<<std::setw(span)<<val[iw];
			if(iw<bin1-1) facc<<", ";
			if ( !((iw+1)%NCol) || iw==bin1-1) facc<<"\n"; 				

		}
		facc<<"};// end of const double "<<table<<"["<<bin1<<"]"<<std::endl;
		facc.close();
	}
	//////////////////////////////////////////////////////////////////////////////////////////////  
	//print ascii file with 5-D Array or 5-D Array pointer
	template <class T_5DArr>
	void CreateFile(T_5DArr val, int bin1, int bin2, int bin3,int bin4, int bin5, char **name=0, int ncol=8, bool bPrintIndex=false, int precision=6)
	{
		//write out this table
		char *DefaultName[]={"Table_4D","Bin1","Bin2","Bin3","Bin4","Bin5"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.txt",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(val[0][0][0][0][0])==typeid(int) || typeid(val[0][0][0][0][0])==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(val[0][0][0][0][0])==typeid(float) ||  typeid(val[0][0][0][0][0])==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 5-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"]["<<name[3]<<"]["<<name[4]<<"]["<<name[5]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";
		facc<<"//static const int "<<name[3]<<"="<<bin3<<"; // bin3 \n";
		facc<<"//static const int "<<name[4]<<"="<<bin4<<"; // bin4 \n";
		facc<<"//static const int "<<name[5]<<"="<<bin5<<"; // bin5 \n";

		int NCol = (bin4>ncol)?ncol:bin4;

		if(bPrintIndex) facc<<"//SuperIndex=001002003004005, = bin1*pow(10,12)+bin2*pow(10,9)+bin3*pow(10,6)+bin4*pow(10,3)+bin5 "<<std::endl; 
		char stridx[20];
		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)
			{  
				for(int it=0;it<bin3;it++)
				{
					for(int ip=0;ip<bin4;ip++)
					{
						for(int ie=0;ie<bin5;ie++)
						{
							if(bPrintIndex && !(ie%NCol)) 
							{ 
								sprintf(stridx,"%03d%03d%03d%03d%03d",iw,iq,it,ip,ie);
								facc<<stridx<<"\t";
							}
							facc<<std::setw(span)<<val[iw][iq][it][ip][ie];
							if(ie<bin5-1) facc<<", ";
							if ( !((ie+1)%NCol) || ie==bin5-1) facc<<"\n"; 
						}	
					}
				}  
			} 
		}
		facc.close();
	}    

	//////////////////////////////////////////////////////////////////////////////////////////////  
	//print ascii file with 4-D Array or 4-D Array pointer
	template <class T_4DArr>
	void CreateFile(T_4DArr val, int bin1, int bin2, int bin3,int bin4, char **name=0, int ncol=8, bool bPrintIndex=false, int precision=6)
	{
		//write out this table
		char *DefaultName[]={"Table_4D","Bin1","Bin2","Bin3","Bin4"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.txt",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(val[0][0][0][0])==typeid(int) || typeid(val[0][0][0][0])==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(val[0][0][0][0])==typeid(float) ||  typeid(val[0][0][0][0])==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 4-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"]["<<name[3]<<"]["<<name[4]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";
		facc<<"//static const int "<<name[3]<<"="<<bin3<<"; // bin3 \n";
		facc<<"//static const int "<<name[4]<<"="<<bin4<<"; // bin4 \n";

		int NCol = (bin4>ncol)?ncol:bin4;

		if(bPrintIndex) facc<<"//SuperIndex=001002003004, = bin1*pow(10,9)+bin2*pow(10,6)+bin3*pow(10,3)+bin4 "<<std::endl; 
		char stridx[20];
		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)
			{  
				for(int it=0;it<bin3;it++)
				{
					for(int ip=0;ip<bin4;ip++)
					{
						if(bPrintIndex && !(ip%NCol)) 
						{ 
							sprintf(stridx,"%03d%03d%03d%03d",iw,iq,it,ip);
							facc<<stridx<<"\t";
						}
						facc<<std::setw(span)<<val[iw][iq][it][ip];
						if(ip<bin4-1) facc<<", ";
						if ( !((ip+1)%NCol) || ip==bin4-1) facc<<"\n"; 
					}	
				}  
			} 
		}
		facc.close();
	}    

	/////////////////////////////////////////////////////////////////////////////      
	//print ascii file with 3-D Array or 3-D Array pointer
	template <class T_3DArr>
	void CreateFile(T_3DArr val, int bin1, int bin2, int bin3, char **name=0, int ncol=8, bool bPrintIndex=false, int precision=6)
	{
		//write out this table
		char *DefaultName[]={"Table_3D","Bin1","Bin2","Bin3"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.txt",table);
		std::ofstream facc;
		facc.open(filename); 
		facc.precision(precision);
		int span=(typeid(val[0][0][0])==typeid(int) || typeid(val[0][0][0])==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(val[0][0][0])==typeid(float) ||  typeid(val[0][0][0])==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 3-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"]["<<name[3]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";
		facc<<"//static const int "<<name[3]<<"="<<bin3<<"; // bin3 \n";

		int NCol = (bin3>ncol)?ncol:bin3;

		if(bPrintIndex) facc<<"//SuperIndex=001002003, = bin1*pow(10,6)+bin2*pow(10,3)+bin3 "<<std::endl; 
		char stridx[20];
		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)
			{  
				for(int it=0;it<bin3;it++)
				{
					if(bPrintIndex && !(it%NCol)) 
					{ 
						sprintf(stridx,"%03d%03d%03d",iw,iq,it);
						facc<<stridx<<"\t";
					}
					facc<<std::setw(span)<<val[iw][iq][it];
					if(it<bin3-1) facc<<", ";
					if ( !((it+1)%NCol) || it==bin3-1) facc<<"\n"; 						
				}  
			} 
		}
		facc.close();
	}    


	/////////////////////////////////////////////////////////////////////////////  
	//print ascii file with 2-D Array or 2-D Array pointer
	template <class T_2DArr>
	void CreateFile(T_2DArr val, int bin1, int bin2, char **name=0, int ncol=8, bool bPrintIndex=false, int precision=6)
	{
		//write out this table
		char *DefaultName[]={"Table_2D","Bin1","Bin2"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.txt",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(val[0][0])==typeid(int) || typeid(val[0][0])==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(val[0][0])==typeid(float) ||  typeid(val[0][0])==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 2-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";

		int NCol = (bin2>ncol)?ncol:bin2;

		if(bPrintIndex) facc<<"//SuperIndex=001002, = bin1*pow(10,3)+bin2 "<<std::endl; 
		char stridx[20];
		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)
			{  
				if(bPrintIndex && !(iq%NCol)) 
				{ 
					sprintf(stridx,"%03d%03d",iw,iq);
					facc<<stridx<<"\t";
				}
				facc<<std::setw(span)<<val[iw][iq];
				if(iq<bin2-1) facc<<", ";
				if ( !((iq+1)%NCol) || iq==bin2-1) facc<<"\n"; 					 
			} 
		}
		facc.close();
	}    
	/////////////////////////////////////////////////////////////////////////////    
	//print ascii file with 1-D Array or 1-D Array pointer
	template <class T_1DArr>
	void CreateFile(T_1DArr val, int bin1, char **name=0, int ncol=8, bool bPrintIndex=false, int precision=6)
	{
		//write out this table
		char *DefaultName[]={"Table_1D","Bin1"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.txt",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(precision);
		int span=(typeid(val[0])==typeid(int) || typeid(val[0])==typeid(long)) ? 6 : 6+precision;

		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		if( typeid(val[0])==typeid(float) ||  typeid(val[0])==typeid(double))
			facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 1-d array "<<table<<"["<<name[1]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";

		int NCol = (bin1>ncol)?ncol:bin1;

		if(bPrintIndex) facc<<"//SuperIndex=### = bin1 "<<std::endl; 
		char stridx[20];
		for(int iw=0;iw<bin1;iw++)
		{  
			if(bPrintIndex && !(iw%NCol)) 
			{ 
				sprintf(stridx,"%03",iw);
				facc<<stridx<<"\t";
			}
			facc<<std::setw(span)<<val[iw];
			if(iw<bin1-1) facc<<", ";
			if ( !((iw+1)%NCol) || iw==bin1-1) facc<<"\n"; 			 
		}
		facc.close();
	}    
	/////////////////////////////////////////////////////////////////////////////
	//input, pointer to 2 buffer, can be 1D to 5D arrays, used bin#>0 to find out the dimention
	//return a pointer to the ratio array
	template <typename T1, typename T2>
	double * CalculateRatio_T(T1 *numerator, T2 *denominator, int bin1, int bin2=0,int bin3=0,int bin4=0,int bin5=0)
	{
		double *result=0;
		double dDef=0.0;
		int dimention=1;
		if(bin5>0)      dimention=5;
		else if(bin4>0) dimention=4;
		else if(bin3>0) dimention=3;
		else if(bin2>0) dimention=2;
		if(bin1<1) 
		{
			cout<<"Error! template CalculateRatio_T(*) argument bin1 invalid!"<<endl;
			return NULL;
		}

		int index=0;
		if(dimention==5)
		{
			double *****Acc=IniDynamicArray( bin1, bin2, bin3, bin4,bin5, dDef);
			for(int iw=0;iw<bin1;iw++)
			{  
				for(int iq=0;iq<bin2;iq++)
				{  
					for(int it=0;it<bin3;it++)
					{
						for(int ip=0;ip<bin4;ip++) 
						{
							for(int i5=0;i5<bin4;i5++) 
							{
								index=iw+iq*bin2+it*bin3+ip*bin4+i5*bin5;
								Acc[iw][iq][it][ip][i5]=(denominator[index]==0)?0.0:(double)numerator[index]/(double)denominator[index];
							}
						}
					}  
				} 
			}
			result=&Acc[0][0][0][0][0];
		}
		else if(dimention==4)
		{
			double ****Acc=IniDynamicArray( bin1, bin2, bin3, bin4, dDef);
			for(int iw=0;iw<bin1;iw++)
			{  
				for(int iq=0;iq<bin2;iq++)
				{  
					for(int it=0;it<bin3;it++)
					{
						for(int ip=0;ip<bin4;ip++) 
						{                                
							index=iw+iq*bin2+it*bin3+ip*bin4;
							Acc[iw][iq][it][ip]=(denominator[index]==0)?0.0:(double)numerator[index]/(double)denominator[index];
						}
					}  
				} 
			}
			result=&Acc[0][0][0][0];
		}
		else  if(dimention==3)
		{
			double ***Acc=IniDynamicArray( bin1, bin2, bin3, dDef);
			for(int iw=0;iw<bin1;iw++)
			{  
				for(int iq=0;iq<bin2;iq++)
				{  
					for(int it=0;it<bin3;it++)
					{
						index=iw+iq*bin2+it*bin3;
						Acc[iw][iq][it]=(denominator[index]==0)?0.0:(double)numerator[index]/(double)denominator[index];
					}  
				} 
			}
			result=&Acc[0][0][0];
		}
		else  if(dimention==2)
		{
			double **Acc=IniDynamicArray( bin1, bin2, dDef);
			for(int iw=0;iw<bin1;iw++)
			{  
				for(int iq=0;iq<bin2;iq++)
				{  
					index=iw+iq*bin2;
					Acc[iw][iq]=(denominator[index]==0)?0.0:(double)numerator[index]/(double)denominator[index];				   
				} 
			}
			result=&Acc[0][0];
		}
		else
		{
			double *Acc=IniDynamicArray( bin1, dDef);
			for(int iw=0;iw<bin1;iw++)
			{  
				Acc[iw]=(denominator[iw]==0)?0.0:numerator[iw]/denominator[iw];	   			   
			}
			result=&Acc[0];
		}
		return result;
	}

	/////////////////////////////////////////////////////////////////////////////
	//return the ratio array
	template <typename T1, typename T2>
	double ***** CalculateRatio(T1 *****numerator, T2 *****denominator, int bin1, int bin2,int bin3,int bin4,int bin5)
	{
		double dDef=0.0;
		double *****Acc=IniDynamicArray( bin1, bin2, bin3, bin4,bin5, dDef);

		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)
			{  
				for(int it=0;it<bin3;it++)
				{
					for(int ip=0;ip<bin4;ip++) 
					{
						for(int i5=0;i5<bin4;i5++) 
						{
							Acc[iw][iq][it][ip][i5]=(denominator[iw][iq][it][ip][i5]==0) ? 
								0.0 : (double)numerator[iw][iq][it][ip][i5]/(double)denominator[iw][iq][it][ip][i5];
						}
					}
				}  
			} 
		}
		return Acc;
	}

	template <typename T1, typename T2>
	double **** CalculateRatio(T1 ****numerator, T2 ****denominator, int bin1, int bin2,int bin3,int bin4)
	{
		double dDef=0.0;
		double ****Acc=IniDynamicArray( bin1, bin2, bin3, bin4, dDef);

		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)
			{  
				for(int it=0;it<bin3;it++)
				{
					for(int ip=0;ip<bin4;ip++) 
					{                                
						Acc[iw][iq][it][ip]=(denominator[iw][iq][it][ip]==0) ? 
							0.0 : (double)numerator[iw][iq][it][ip]/(double)denominator[iw][iq][it][ip];
					}
				}  
			} 
		}
		return Acc;
	}

	template <typename T1, typename T2>
	double *** CalculateRatio(T1 ***numerator, T2 ***denominator, int bin1, int bin2,int bin3)
	{
		double dDef=0.0;
		double ***Acc=IniDynamicArray( bin1, bin2, bin3, dDef);

		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)
			{  
				for(int it=0;it<bin3;it++)
				{
					Acc[iw][iq][it]=(denominator[iw][iq][it]==0) ? 
						0.0 : (double)numerator[iw][iq][it]/(double)denominator[iw][iq][it];
				}  
			} 
		}
		return Acc;
	}

	template <typename T1, typename T2>
	double ** CalculateRatio(T1 **numerator, T2 **denominator, int bin1, int bin2)
	{
		double dDef=0.0;
		double **Acc=IniDynamicArray( bin1, bin2, dDef);

		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)
			{  
				Acc[iw][iq]=(denominator[iw][iq]==0) ? 0.0 : (double)numerator[iw][iq]/(double)denominator[iw][iq];				   
			} 
		}
		return Acc;
	}


	template <typename T1, typename T2>
	double * CalculateRatio(T1 *numerator, T2 *denominator, int bin1)
	{
		double dDef=0.0;
		double *Acc=IniDynamicArray( bin1, dDef);

		for(int iw=0;iw<bin1;iw++)
		{  
			Acc[iw]=(denominator[iw]==0) ? 0.0 : (double)numerator[iw]/(double)denominator[iw];	   			   
		}
		return Acc;
	}

	//////////////////////////////////////////////////////////////////////////////
	//due to typeid(T).name() return different string in WINDOWS and LINUX, I need to identify it here

	/*
	windows vc8			linux gcc v3.4.3	
	int is: int           	int is: i       	
	i is: int           	  i is: i       	
	pi is: int *         	 pi is: Pi      	
	*pi is: int           	*pi is: i       	
	f is: float           	f is: f         	
	d is: double          	d is: d         	
	l is: long            	l is: l         	
	ui is: unsigned int   	ui is: j        	
	c is: char            	c is: c         	
	sz is: char [5]       	sz is: A5_c     	
	pi_2d is: int * *       pi_2d is: PPi  	
	**pi_2d is: int         **pi_2d is: i   	
	pi_3d is: int * * *     pi_3d is: PPPi 	
	***pi_3d is: int        ***pi_3d is: i  	
	pi_4d is: int * * * *   pi_4d is: PPPPi	
	****pi_4d is: int       ****pi_4d is: i 	

	*/


	template <typename T>
	const char *GetTypeName(T var)
	{   
		//I need to put it as static otherwise can not return correctly
		static std::string sysType; 

		if (typeid(var)==typeid(int)) sysType="int";
		else if (typeid(var)==typeid(float)) sysType="float";
		else if (typeid(var)==typeid(double)) sysType="double";
		else if (typeid(var)==typeid(long)) sysType="long";
		else if (typeid(var)==typeid(unsigned int)) sysType="unsigned int";
		else if (typeid(var)==typeid(short)) sysType="short";
		else if (typeid(var)==typeid(char)) sysType="char";
		else  sysType="unknow_type";

		return sysType.c_str();
	}
	//////////////////////////////////////////////////////////////////////////////

	//print ascii file with 4-D Array     
	template <class T_4DArr>
	void PrintBuffer(std::ostream &pOut, T_4DArr val, int bin1, int bin2, int bin3, int bin4, int ncol)
	{
		int precision=pOut.precision();
		int span=(typeid(val[0][0][0][0])==typeid(int) || typeid(val[0][0][0][0])==typeid(long)) ? 6 : 6+precision;
		int NCol = (bin4>ncol)?ncol:bin4;
		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)
			{  
				for(int it=0;it<bin3;it++)
				{
					for(int ip=0;ip<bin4;ip++)
					{
						pOut<<std::setw(span)<<val[iw][iq][it][ip];
						if(ip<bin4-1) pOut<<", ";
						if ( !((ip+1)%NCol) || ip==bin4-1) pOut<<"\n"; 	
					}
				} 
				pOut<<"\n"; 
			} 
			pOut<<"\n";
		}
	}    

	//print ascii file with 3-D Array 
	template <class T_3DArr>
	void PrintBuffer(std::ostream &pOut, T_3DArr val, int bin1, int bin2, int bin3, int ncol)
	{
		int precision=pOut.precision();
		int span=(typeid(val[0][0][0])==typeid(int) || typeid(val[0][0][0])==typeid(long)) ? 6 : 6+precision;
		int NCol = (bin3>ncol)?ncol:bin3;
		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)
			{  
				for(int it=0;it<bin3;it++)
				{
					pOut<<std::setw(span)<<val[iw][iq][it];
					if(it<bin3-1) pOut<<", ";
					if ( !((it+1)%NCol) || it==bin3-1) pOut<<"\n"; 						
				};
			} 
			pOut<<"\n";
		}
	}    


	//print ascii file with 2-D Array   
	template <class T_2DArr>
	void PrintBuffer(std::ostream &pOut, T_2DArr val, int bin1, int bin2, int ncol)
	{		
		int precision=pOut.precision();
		int span=(typeid(val[0][0])==typeid(int) || typeid(val[0][0])==typeid(long)) ? 6 : 6+precision;
		int NCol = (bin2>ncol)?ncol:bin2;
		for(int iw=0;iw<bin1;iw++)
		{  
			for(int iq=0;iq<bin2;iq++)
			{  
				pOut<<std::setw(span)<<val[iw][iq];
				if(iq<bin2-1) pOut<<", ";
				if ( !((iq+1)%NCol) || iq==bin2-1) pOut<<"\n"; 					  
			} 
			pOut<<"\n";
		}
	}    

	//print ascii file with 1-D Array 
	template <class T_1DArr>
	void PrintBuffer(std::ostream &pOut, T_1DArr val, int bin1, int ncol)
	{
		int precision=pOut.precision();
		int span=(typeid(val[0])==typeid(int) || typeid(val[0])==typeid(long)) ? 6 : 6+precision;
		int NCol = (bin1>ncol)?ncol:bin1;
		for(int iw=0;iw<bin1;iw++)
		{  
			pOut<<std::setw(span)<<val[iw];
			if(iw<bin1-1) pOut<<", ";
			if ( !((iw+1)%NCol) || iw==bin1-1) pOut<<"\n"; 									
		}
	}    

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	template <typename  T1, typename T2>
	void PrintACCTable(const double bin1[], int n1, const double bin2[], int n2, const double bin3[], int n3,
		const double bin4[], int n4, int level, const char *filename, char *headblock, T1 ****N_det, T2 ****N_true)
	{
		//I also generate the regular table for other people
		std::ofstream fout;
		fout.open(filename);
		if(!headblock)
		{
			fout<<"//This is a 4-D acceptance table:"<<std::endl;  
			fout<<"//static const int bin1="<<n1<<std::endl;
			fout<<"//static const int bin2="<<n2<<std::endl;
			fout<<"//static const int bin3="<<n3<<std::endl;
			fout<<"//static const int bin4="<<n4<<std::endl;

			fout<<"This table will not be printed out in the following ways in order to save space: "<<std::endl
				<<" LEVEL >=3, the full table will be printed, "<<std::endl
				<<" LEVEL ==2, skip bins which thrown# is zero "<<std::endl
				<<" LEVEL ==1, skip bins which detected# or thrown# is zero  "<<std::endl
				<<" LEVEL ==0, don't print bin variables, use super index and skip bins which detected# or thrown# is zero  "<<std::endl
				<<"            where SuperIndex=001002003004, = bin1*pow(10,9)+bin2*pow(10,6)+bin3*pow(10,3)+bin4 "<<std::endl; 

			fout<<"This table is generated with LEVEL "<<level<<std::endl;

			if(level) fout<<"//bin1_l    bin2_l   bin3_l   bin4_l  Thrown#     Det#   Acc=Det#/Thrown#"<<std::endl;
			else fout<<"//SuperIndex  Thrown#   Det#   Acc=Det#/Thrown#"<<std::endl;
		}
		else fout<<headblock;

		int span=7;
		int iw=n1,iq=n2,it=n3,ip=n4;
		char str[255];
		for( iw=0;iw<n1;iw++)
		{  
			for( iq=0;iq<n2;iq++)
			{  
				for( it=0;it<n3;it++)
				{
					for( ip=0;ip<n4;ip++)
					{					
						double acc=(N_true[iw][iq][it][ip]==0)? 0.0:(double)N_det[iw][iq][it][ip]/(double)N_true[iw][iq][it][ip];
						if(level>=3) 
						{
							fout<<setw(span)<<bin1[iw]<<"  "<<setw(span)<<bin2[iq]<<"  "<<setw(span)<< bin3[it]<<"  "<<setw(span)<<bin4[ip]
							<<"  "<<setw(span)<< N_true[iw][iq][it][ip]<<"  "<<setw(span)<< N_det[iw][iq][it][ip]
							<<"  "<<setw(span)<<acc<<endl;
						}
						else if(level==2) 
						{
							if(N_true[iw][iq][it][ip]!=0)
							{
								fout<<setw(span)<<bin1[iw]<<"  "<<setw(span)<<bin2[iq]<<"  "<<setw(span)<< bin3[it]<<"  "<<setw(span)<<bin4[ip]
								<<"  "<<setw(span)<< N_true[iw][iq][it][ip]<<"  "<<setw(span)<< N_det[iw][iq][it][ip]
								<<"  "<<setw(span)<<acc<<endl;
							}			 
						}else if(level==1) 
						{
							if(N_det[iw][iq][it][ip]!=0  &&  N_true[iw][iq][it][ip]!=0)
							{
								fout<<setw(span)<<bin1[iw]<<"  "<<setw(span)<<bin2[iq]<<"  "<<setw(span)<< bin3[it]<<"  "<<setw(span)<<bin4[ip]
								<<"  "<<setw(span)<< N_true[iw][iq][it][ip]<<"  "<<setw(span)<< N_det[iw][iq][it][ip]
								<<"  "<<setw(span)<<acc<<endl;
							}			 
						}
						else 
						{
							if(N_det[iw][iq][it][ip]!=0  &&  N_true[iw][iq][it][ip]!=0)
							{
								sprintf(str,"%03d%03d%03d%03d",iw,iq,it,ip);
								fout<<str<<"  "<<setw(span)<< N_true[iw][iq][it][ip]<<"  "<<setw(span)<< N_det[iw][iq][it][ip]
								<<"  "<<setw(span)<<acc<<endl;
							}			 
						}
					}	  
				}  
			} 
		}
		//print last boundary
		if(level>0) 
		{
			sprintf(str,"%5.3f  %5.3f  %5.3f  %5.3f  %7.0f  %7.0f  %7.4f",
				bin1[iw], bin2[iq], bin3[it],  bin4[ip],  0.0,  0.0,  0.0 );
			fout<<str<<std::endl;
		}
		//////////////////////////////////////////////////////////////////////////
		fout<<std::endl;
		fout.close();
	}    
	/////////////////////////////////////////////////////////////////////////////////////


	//input: 2 arrays of 4 elements which used to define what binindex will be removed
	template <class T_4DArr, typename T>
    T**** RemoveBins(T_4DArr arr, T datatype, int bin1, int bin2, int bin3,int bin4, 
		int removestart[4], int removeend[4], char **name=0,int precision=6)
    {   
        T defaultval=T(0.0);
		
		int newbin[]={bin1,bin2,bin3,bin4}; 

		for(int i=0;i<4;i++) 
		{
			//check if the user given a wrong end index
			if( removestart[i]>=0 && removestart[i]<=removeend[i] && removeend[i]<newbin[i])
			{
				newbin[i]=newbin[i]-(removeend[i]-removestart[i]+1);
			}
			if (newbin[i]==0)
			{
				cout<<"Error! RemoveBins() can not shrink dimention. Check your removestart and removeend\n";
				return NULL;
			}
		}
        T ****newarr=IniDynamicArray( newbin[0], newbin[1], newbin[2], newbin[3],defaultval);
#ifdef ACCTOOLS_DEBUG
		defaultval=T(-1.0);
        T ****testarr=IniDynamicArray(bin1,bin2,bin3,bin4,defaultval);
#endif

		int newi1=-1;
		for(int i1=0;i1<bin1;i1++)
		{  
			if(i1>=removestart[0] && i1<=removeend[0]) continue;
			newi1++;
			int newi2=-1;
			for(int i2=0;i2<bin2;i2++)
			{  
				if(i2>=removestart[1] && i2<=removeend[1]) continue;
				newi2++;
				int newi3=-1;
				for(int i3=0;i3<bin3;i3++)
				{                    
					if(i3>=removestart[2] && i3<=removeend[2]) continue;
					newi3++;
					int newi4=-1;
					for(int i4=0;i4<bin4;i4++) 
					{
						if(i4>=removestart[3] && i4<=removeend[3]) continue;
						newi4++;
						newarr[newi1][newi2][newi3][newi4]=T(arr[i1][i2][i3][i4]);
#ifdef ACCTOOLS_DEBUG
						testarr[i1][i2][i3][i4]=T(arr[i1][i2][i3][i4]);
#endif
					}
				}
			}
		}

		//crate new inc file ant txt file
		CreateIncFile(newarr,newbin[0],newbin[1], newbin[2], newbin[3],name,newbin[3],precision);
		CreateFile(newarr,newbin[0],newbin[1], newbin[2], newbin[3],name,newbin[3],true,precision);
		
#ifdef ACCTOOLS_DEBUG
		name[0]="testarr";
		CreateFile(testarr,bin1,bin2,bin3,bin4,name,bin4,true,precision);
#endif
		return newarr;
	}

	
	template <class T_4DArr, typename T>
    T**** Remove1Bin(T_4DArr arr, T datatype, int bin1, int bin2, int bin3,int bin4, 
		int bin2remove,int removestart, int removeend, char **name=0,int precision=6)
    {   
		int start[]={-1,-1,-1,-1};
		int end[]={-1,-1,-1,-1};
		start[bin2remove]=removestart;
		end[bin2remove]=removeend;
		return RemoveBins(arr,datatype,bin1,bin2,bin3,bin4,start,end,name,precision);
	}
	///////////////////////////////////////////////////////////////////////////////////////////

};

#endif //_ACCEPTANCE_TOOLS_

