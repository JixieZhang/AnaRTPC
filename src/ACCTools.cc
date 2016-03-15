#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <string>

#include "ACCTools.h"

namespace ACCEPTANCE {

    //remove all occurance of a substring from the orignal string
    void RemoveAll(std::string str, std::string key)
    {   
        size_t found;
        found=str.rfind(key);
        while (found!=string::npos)
        {
            str.erase(found,key.length());
        }
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	void CreateQ2Table(double Q2Table[])
	{
		//double Q2Table[50]
		double Q2=0,Q2const=(1.0 + pow(10.0,(-1.0/13.0)))/2.0;

		FILE* myfile = fopen("Q2BinTable.txt","w");
		if(myfile==NULL) return;

		fprintf(myfile,"Index   Q2\n");
		Q2Table[0]=Q2;
		fprintf(myfile,"%2d   %7.4f\n",0,Q2);
		int i=1;
		do{
			if(i>=50) break;
			Q2=Q2const*pow(10.0,((i-27.0)/13.0)); 
			//if(Q2<0.01) continue;
			Q2Table[i]=Q2;   
			fprintf(myfile,"%2d   %7.4f\n",i,Q2);
			i++;
		}while(Q2<5.0);

		fprintf(myfile,"static const double kQ2Bin[%d]={",i);
		for (int j=0;j<i;j++)
		{	  
			if(!(j%10)) fprintf(myfile,"\n\t");
			fprintf(myfile,"%7.4f",Q2Table[j]);
			if(j<i-1) fprintf(myfile,", ");
			else fprintf(myfile,"\n");
		}
		fprintf(myfile,"};\n");
		fclose(myfile);
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	int BinarySearch(const int sortedArray[], int first, int last, const int key) 
	{
		// function:
		//   Searches sortedArray[first]..sortedArray[last] for key.  
		// returns: index of the matching element if it finds key, 
		//          otherwise  -index (after which it could be inserted).
		//			which satisfy sortedArray[index]< key < sortedArray[index+1].
		// parameters:
		//   sortedArray in  array of sorted (ascending) values.
		//   first, last in  lower and upper subscript bounds
		//   key         in  value to search for.
		// returns:
		//   index of key, or -insertion_position if key is not 
		//                 in the array.
		if(key<sortedArray[0] || key>sortedArray[last-1]) return -last+1;	   
		while (first <= last) {
			int mid = (first + last) / 2;  // compute mid point.
			if (key > sortedArray[mid]) 
				first = mid + 1;  // repeat search in top half.
			else if (key < sortedArray[mid]) 
				last = mid - 1; // repeat search in bottom half.
			else
				return mid;     // found it. return position /////
		}
		return -first+1;    // failed to find key
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	int BinarySearch(const double sortedArray[], int first, int last, const double key) 
	{
		// function:
		//   Searches sortedArray[first]..sortedArray[last] for key.  
		// returns: index of the matching element if it finds key, 
		//          otherwise  -index (after which it could be inserted).
		//			which satisfy sortedArray[index]< key < sortedArray[index+1].
		// parameters:
		//   sortedArray in  array of sorted (ascending) values.
		//   first, last in  lower and upper subscript bounds
		//   key         in  value to search for.
		// returns:
		//   index of key, or -insertion_position if key is not 
		//                 in the array. 
		if(key<sortedArray[0] || key>sortedArray[last-1]) return -last+1;	 	   
		while (first <= last) 
		{
			int mid = (first + last) / 2;  // compute mid point.
			if (key > sortedArray[mid]) 
				first = mid + 1;  // repeat search in top half.
			else if (key < sortedArray[mid]) 
				last = mid - 1;  // repeat search in bottom half.
			else
				return mid;      // found it. return position /////
		}
		return -first+1;		// failed to find key
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////
	void PrintIncFile(double ****val, int bin1, int bin2, int bin3,int bin4, char **name, int ncol)
	{
		//write out this table
		char *DefaultName[]={"Table_4D","Bin1","Bin2","Bin3","Bin4"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(6);
		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 4-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"]["<<name[3]<<"]["<<name[4]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";
		facc<<"//static const int "<<name[3]<<"="<<bin3<<"; // bin3 \n";
		facc<<"//static const int "<<name[4]<<"="<<bin4<<"; // bin4 \n";

		int NCol = (bin4>ncol)?ncol:bin4;

		facc<<"const double "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]["<<bin4<<"]={"<<std::endl;
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
						facc<<std::setw(8)<<val[iw][iq][it][ip];
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
	void PrintIncFile(double ***val, int bin1, int bin2, int bin3, char **name, int ncol)
	{
		//write out this table
		char *DefaultName[]={"Table_3D","Bin1","Bin2","Bin3"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(6);
		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		facc.setf(std::ios_base::fixed,std::ios_base::floatfield);  

		facc<<"//This is the table for 3-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"]["<<name[3]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";
		facc<<"//static const int "<<name[3]<<"="<<bin3<<"; // bin3 \n";

		int NCol = (bin3>ncol)?ncol:bin3;
		facc<<"const double "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]={"<<std::endl;
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
					facc<<std::setw(8)<<val[iw][iq][it];
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
	void PrintIncFile(double **val, int bin1, int bin2, char **name, int ncol)
	{  
		char *DefaultName[]={"Table_2D","Bin1","Bin2"};
		if(name==0) name=DefaultName;
		//write out this table
		char *table=name[0];
		char filename[255];
		sprintf(filename,"%s.inc",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(6);
		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 2-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";


		int NCol = (bin2>ncol)?ncol:bin2;
		facc<<"const double "<<table<<"["<<bin1<<"]["<<bin2<<"]={"<<std::endl;
		for(int iw=0;iw<bin1;iw++)
		{  
			facc<<"\t{ //start of bin1 "<<name[1]<<"="<<iw<<std::endl;
			for(int iq=0;iq<bin2;iq++)
			{  
				if (!(iq%NCol)) facc<<"\t\t"; 
				facc<<std::setw(8)<<val[iw][iq];
				if(iq<bin2-1) facc<<", ";
				if ( !((iq+1)%NCol) || iq==bin2-1) facc<<"\n"; 				
			} 
			if(iw<bin1-1) facc<<"\t}, //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
			else facc<<"\t} //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
		}
		facc<<"};// end of const double "<<table<<"["<<bin1<<"]["<<bin2<<"]"<<std::endl;
		facc.close();
	}

	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    //print either inc or ascii file for given 4-D array
	void PrintTable(double ****val, int bin1, int bin2, int bin3,int bin4, char **name, int ncol, 
        bool bPrintIncFile, bool bPrintIndex)
	{
		//write out this table
		char *DefaultName[]={"Table_4D","Bin1","Bin2","Bin3","Bin4"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
          
		if(bPrintIncFile) sprintf(filename,"%s.inc",table);
        else sprintf(filename,"%s.txt",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(6);
		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 4-d array "<<table<<"["<<name[1]<<"]["<<name[2]<<"]["<<name[3]<<"]["<<name[4]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";
		facc<<"//static const int "<<name[2]<<"="<<bin2<<"; // bin2 \n";
		facc<<"//static const int "<<name[3]<<"="<<bin3<<"; // bin3 \n";
		facc<<"//static const int "<<name[4]<<"="<<bin4<<"; // bin4 \n";
        
        char stridx[20];
        if(!bPrintIncFile && bPrintIndex) facc<<"//SuperIndex=001002003004, = bin1*pow(10,9)+bin2*pow(10,6)+bin3*pow(10,3)+bin4 "<<std::endl; 

		int NCol = (bin4>ncol)?ncol:bin4;
      
        
		if(bPrintIncFile) facc<<"const double "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]["<<bin4<<"]={"<<std::endl;
		for(int iw=0;iw<bin1;iw++)
		{  
			if(bPrintIncFile) facc<<"\t{ //start of bin1 "<<name[1]<<"="<<iw<<std::endl;
			for(int iq=0;iq<bin2;iq++)
			{  
				if(bPrintIncFile) facc<<"\t\t{ //start of bin2 "<<name[2]<<"="<<iq<<std::endl;
				for(int it=0;it<bin3;it++)
				{
					if(bPrintIncFile)
                    {
                        facc<<"\t\t\t{ ";
                        if(bin4/NCol>10) facc<<" //start of bin3 "<<name[3]<<"="<<it;
                        facc<<std::endl;
                    }
					for(int ip=0;ip<bin4;ip++)
					{
                        if(!(ip%NCol)) 
                        {
                            if (bPrintIncFile) facc<<"\t\t\t\t"; 
                            else if(bPrintIndex)
                            {
                                sprintf(stridx,"%03d%03d%03d%03d",iw,iq,it,ip);
                                facc<<stridx<<"\t";
                            }
                        }
						facc<<std::setw(8)<<val[iw][iq][it][ip];
						if(ip<bin4-1) facc<<", ";
						if ( !((ip+1)%NCol) || ip==bin4-1) facc<<"\n"; 
					}	

                    if(bPrintIncFile) 
                    {
                        if(it<bin3-1) facc<<"\t\t\t},";
                        else facc<<"\t\t\t}";
                        if(bin4/NCol>10) facc<<" //  end of bin3 "<<name[3]<<"="<<it;
                        facc<<std::endl;
                    }
				}  
				if(bPrintIncFile) 
                {
                    if(iq<bin2-1) facc<<"\t\t}, //  end of bin2 "<<name[2]<<"="<<iq<<std::endl;
				    else facc<<"\t\t} //  end of bin2 "<<name[2]<<"="<<iq<<std::endl;
                }
			} 
			if(bPrintIncFile) 
            {
                if(iw<bin1-1) facc<<"\t}, //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
			    else facc<<"\t} //  end of bin1 "<<name[1]<<"="<<iw<<std::endl;
            }
		}
		if(bPrintIncFile) facc<<"};// end of const double "<<table<<"["<<bin1<<"]["<<bin2<<"]["<<bin3<<"]["<<bin4<<"]"<<std::endl;
		facc.close();
	}    

/////////////////////////////////////////////////////////////////////////////////////////////////////////

	void PrintFile(double ****val, int bin1, int bin2, int bin3,int bin4, char **name, int ncol, bool bPrintIndex)
	{
		//write out this table
		char *DefaultName[]={"Table_4D","Bin1","Bin2","Bin3","Bin4"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
        sprintf(filename,"%s.txt",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(6);
		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		//facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

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
						facc<<std::setw(8)<<val[iw][iq][it][ip];
						if(ip<bin4-1) facc<<", ";
						if ( !((ip+1)%NCol) || ip==bin4-1) facc<<"\n"; 
					}	
				}  
			} 
		}
		facc.close();
	}    


	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	void PrintFile(double ***val, int bin1, int bin2, int bin3, char **name, int ncol, bool bPrintIndex)
	{
		//write out this table
		char *DefaultName[]={"Table_3D","Bin1","Bin2","Bin3"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
        sprintf(filename,"%s.txt",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(6);
		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		//facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

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
					facc<<std::setw(8)<<val[iw][iq][it];
					if(it<bin3-1) facc<<", ";
					if ( !((it+1)%NCol) || it==bin3-1) facc<<"\n"; 						
				}  
			} 
		}
		facc.close();
	}    

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	void PrintFile(double **val, int bin1, int bin2, char **name, int ncol, bool bPrintIndex)
	{
		//write out this table
		char *DefaultName[]={"Table_2D","Bin1","Bin2"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
        sprintf(filename,"%s.txt",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(6);
		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		//facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

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
				facc<<std::setw(8)<<val[iw][iq];
				if(iq<bin2-1) facc<<", ";
				if ( !((iq+1)%NCol) || iq==bin2-1) facc<<"\n"; 	
			} 
		}
		facc.close();
	}    


	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	void PrintFile(double *val, int bin1, char **name, int ncol, bool bPrintIndex)
	{
		//write out this table
		char *DefaultName[]={"Table_1D","Bin1"};
		if(name==0) name=DefaultName;
		char *table=name[0];
		char filename[255];
        sprintf(filename,"%s.txt",table);
		std::ofstream facc;
		facc.open(filename);  
		facc.precision(6);
		// floatfield set to fixed, this = %.precisionf, take precision digit after the point
		//facc.setf(std::ios_base::fixed,std::ios_base::floatfield);   

		facc<<"//This is the table for 1-d array "<<table<<"["<<name[1]<<"].\n";
		facc<<"//static const int "<<name[1]<<"="<<bin1<<"; // bin1 \n";

		int NCol = (bin1>ncol)?ncol:bin1;
              
        if(bPrintIndex) facc<<"//SuperIndex=001 = bin1 "<<std::endl; 
        char stridx[20];
		for(int iw=0;iw<bin1;iw++)
		{  
                if(bPrintIndex && !(iw%NCol)) 
                { 
                    sprintf(stridx,"%03d",iw);
                    facc<<stridx<<"\t";
                }
				facc<<std::setw(8)<<val[iw];
				if(iw<bin1-1) facc<<", ";
				if ( !((iw+1)%NCol) || iw==bin1-1) facc<<"\n"; 			
		}
		facc.close();
	}    

    /////////////////////////////////////////////////////////////////////////////
 
    //return 4-d array
    double **** CalculateACC(int ****N_det, int ****N_true, int bin1, int bin2,int bin3,int bin4)
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
                        Acc[iw][iq][it][ip]=(N_true[iw][iq][it][ip]==0)?0.0:(double)N_det[iw][iq][it][ip]/(double)N_true[iw][iq][it][ip];
				}  
			} 
		}
        return Acc;
    }
    
    /////////////////////////////////////////////////////////////////////////////
    //return 3-d array
    double *** CalculateACC(int ***N_det, int ***N_true, int bin1, int bin2,int bin3)
    {
        double dDef=0.0;
        double ***Acc=IniDynamicArray( bin1, bin2, bin3, dDef);
        for(int iw=0;iw<bin1;iw++)
	    {  
		    for(int iq=0;iq<bin2;iq++)
		    {  
			    for(int it=0;it<bin3;it++)
			    {
                    Acc[iw][iq][it]=(N_true[iw][iq][it]==0)?0.0:(double)N_det[iw][iq][it]/(double)N_true[iw][iq][it];
			    }  
		    } 
	    }
        return Acc;
    }

    //return 2-d array
    double ** CalculateACC(int **N_det, int **N_true, int bin1, int bin2)
    {
        double dDef=0.0;
        double **Acc=IniDynamicArray( bin1, bin2, dDef);
        for(int iw=0;iw<bin1;iw++)
	    {  
		    for(int iq=0;iq<bin2;iq++)
		    {  
                Acc[iw][iq]=(N_true[iw][iq]==0)?0.0:(double)N_det[iw][iq]/(double)N_true[iw][iq];				   
		    } 
	    }
        return Acc;

    }

    //return 1-d array
    double * CalculateACC(int *N_det, int *N_true, int bin1)
    {
        double dDef=0.0;
        double *Acc=IniDynamicArray( bin1, dDef);
        for(int iw=0;iw<bin1;iw++)
	    {  
            Acc[iw]=(N_true[iw]==0)?0.0:(double)N_det[iw]/(double)N_true[iw];	   			   
	    }
        return Acc;
    }

    /////////////////////////////////////////////////////////////////////////////


} //end of namespace ACCEPTANCE



