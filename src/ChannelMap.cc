//////////////////////////////////////////////////////////////////////
// RTPC12 ChannelMap
// Define how the channal id is assigned in location
// Jixie Zhang, v1.0
//////////////////////////////////////////////////////////////////////

#include "ChannelMap.hh"
#include <stdio.h>
#include <math.h>

//#define PAD_LOC_DEBUG 2
//#define CHANNELMAP_DUBUG 0

///////////////////////////Constant Array Definitions///////////////////////////////
//Each connector has 16 channels, which are connected to 4x4 pads
//Each Pad is 4mm(z) x 2.7mm(phi), pads are 4 mil separated in z and 1.5mil in phi
//Basic Assumption:
//1) The connector col index increases along z, row index increases along phi 
//2) Pads are separated in phi by 1.5 mil? need to verify 
//3) Zoffsets({ 1.8, 1.2, 0.6, 0}) need to verify
/**********************************************************************************/
//Note that we have to verify 1) with DAQ to make sure it is written in this way!
/**********************************************************************************/

//kConnChanMap[phi_idx(row)][z_idx(col)] is used to find chan id for given (row,col)
static const int kConnChanMap[4][4] = { 
  { 3, 5, 9, 13}, { 1, 7, 11, 15}, { 2, 6, 8, 12}, { 0, 4, 10, 14}
};
 
//The connector col index increases along z, row index increases along phi 
//Each row (in phi) of connectors (6 of them) will be connected to one adapter board
//kConnColMap[16] will tell the col index for a given chan id
static const int kConnColMap[16] = { 
  0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 
};  //can use (chan%4) instead
  
//kConnColMap[16] will tell the row index for a given chan id
static const int kConnRowMap[16] = { 
  3, 1, 2, 0 , 3, 0, 2, 1, 2, 0, 3, 1, 2, 0, 3, 1
};

/* Rows of pads may be shifted along z... 
 * starting with bottom row and proceeding in increasing phi:
 * By Jixie: I do not get the final number for these yet
  * * *
  repeat
  4th row   no shift 
  3rd row   shifted by 0.6 mm in z 
  2nd row   shifted by 1.2 mm in z 
  1st row   shifted by 1.8 mm in z (downstream)
  repeat
  * * *
So the pad centers of the first pads on each row are at these positions:
*/
static const double kZOffset[4]= { 1.8, 1.2, 0.6, 0.0};
//static const double kZOffset[4]= { 0, 0, 0, 0};

////////////////////////////Constant Array Definitions end /////////////////////////////


static const double PI=acos(0.0)*2.0;
/*=========================================================================*/
ChannelMap::ChannelMap(double pad_w, double pad_l, double pad_r_in, double phi_conn_edge):
  fDeltaA(pad_w), fDeltaZ(pad_l), fPadRi(pad_r_in), fPhiStart(phi_conn_edge)
{
  //suggested RTPC12 parameters
  //double fDeltaZ;//=4.1016;     //a single pad length in z (4mm+4mil gap)
  //double fDeltaA;//=2.7381;     //phi dimension of a pad (the arc length) (2.7mm+1.5mil gap)
  //double fPadRi;   //=79.0;     //inner radius where the pad plane is located (79mm), 10mil thick
  
  //For the above number, fDeltaPhi=1.9950535479061624 degree 
  //the dead area is about 3.5mm in length

  fDeltaPhi=fDeltaA/fPadRi;  //the angular coverage of a single pad (rad)
  
  //define the left edge of z for these 4 rows of the connector connect to pad_id==0 
  for(int i=0;i<4;i++) fZStart[i] = -fDeltaZ*NUM_COLS/2 + kZOffset[i];
  
 #ifdef CHANNELMAP_DUBUG  
    printf("Creating Channel map: %d Rows x %d Columns = %d pads in total\n",NUM_ROWS,NUM_COLS,NUM_PADS);
    printf("\t Width(%.4fmm)  Length(%.4fmm) Inner_S(%.3fmm) fDeltaPhi(%.3fdeg)\n",
      fDeltaA, fDeltaZ, fPadRi, fDeltaPhi*180./PI);
    printf("\t fPhiStart(%.3fdeg) fZStart[4]={%.3f, %.3f, %.3f, %.3f}\n\n",
      fPhiStart*180./PI, fZStart[0], fZStart[1], fZStart[2], fZStart[3]);
 #endif 
 
	for(int i=0;i<NUM_PADS;i++)
	{
		fZPhiMap[i][0]=fZPhiMap[i][1]=-99999.;
	}

	for(int rr=0;rr<NUM_ROWS;rr++)
	{
    for(int cc=0;cc<NUM_COLS;cc++) fChanMap[rr][cc]=-1; 
  }

	GenerateMap();
}

/*=========================================================================*/
ChannelMap::~ChannelMap()
{
	;
}

/*=========================================================================*/
// Set up the pad locations array by creating the locations from knowledge 
//of the pad layout and the readout channel order.
void ChannelMap::GenerateMap()
{
	int pad_id, conn_row, conn_col, pad_seq, pad_col, pad_row;
	double pad_z, pad_phi;

#if defined PAD_LOC_DEBUG && (PAD_LOC_DEBUG>=2)
	printf("RTPC12 Channel Map:\n");
	printf("  ID     z_mm  phi_deg  row  col\n");
#endif

	for(pad_id=0; pad_id<NUM_PADS; pad_id++)
	{
    conn_row = pad_id/384;          //adapter idx, each adapter has 384 chan
    conn_col = (pad_id%384)/16;     //col idx, each adapter connected to 24 connectors
    pad_seq = pad_id%16;            //sequence number within connector: 0-15
    
    pad_row = 4*conn_row + kConnRowMap[pad_seq];
    pad_col = 4*conn_col + kConnColMap[pad_seq];

		pad_z = fZStart[pad_row%4] + fDeltaZ*(pad_col+0.5);
		pad_phi = fPhiStart + fDeltaPhi*(pad_row+0.5);
		if(pad_phi < 0.0)   pad_phi +=2*PI;
		if(pad_phi >= 2*PI) pad_phi -=2*PI;

		//fill the buffer
		fZPhiMap[pad_id][0] = pad_z;
		fZPhiMap[pad_id][1] = pad_phi;
    fChanMap[pad_row][pad_col] = pad_id;
    
    //just to debug routine GetPadID
#if defined PAD_LOC_DEBUG
		int idx = GetPadID(pad_z,pad_phi);
		if(idx!=pad_id || fZPhiMap[idx][0]!=pad_z || fZPhiMap[idx][1]<pad_phi)
		{
			printf("*Error in GetPadID()!! this cell index=%d z=%7.2f phi=%7.2f row=%d col=%d\n",
				   idx, pad_z, pad_phi*180./PI,pad_row,pad_col);
		}
#endif

#if defined PAD_LOC_DEBUG && (PAD_LOC_DEBUG>=2)
		printf("%4d  %7.2f  %7.2f  %3d  %3d\n", pad_id,
			pad_z, pad_phi*180./PI, pad_row, pad_col);
#endif
	}
}

/*====================================================================== */
//return (global) pad row and col using given z_mm and phi_rad
//return -1 if out of range
int ChannelMap::GetPadRowNCol(double z_mm, double phi_rad, int &row, int &col)
{
	//make sure the phi is in range
  if(phi_rad>=2.*PI) phi_rad-=2.*PI;
	if(phi_rad<0.) phi_rad+=2.*PI;

  row = int((phi_rad-fPhiStart)/fDeltaPhi);	 
	col = int((z_mm-fZStart[row%4])/fDeltaZ);

	if(row<0 || row>=NUM_ROWS || col<0 || col>=NUM_COLS)
	{
		//User will give all cracy z and phi to ask for an index
    //necessary to check the result here 
#ifdef CHANNELMAP_DUBUG
		printf("**Error!!! GetRowNCol(z=%3.2fmm phi=%.2fdeg) return row=%3d col=%3d\n",
			z_mm,phi_rad*180./PI,row,col);
#endif      
		if(row<0) row=0;
		if(row>=NUM_ROWS) row=NUM_ROWS-1;
		if(col<0) col=0;
		if(col>=NUM_COLS) col=NUM_COLS-1;
		return -1;
	}
#ifdef CHANNELMAP_DUBUG
	if(CHANNELMAP_DUBUG>=3)
	{
		printf("GetRowNCol(z=%3.2fmm phi=%.4f) ==> row=%3d col=%3d\n",
			z_mm,phi_rad*180./PI,row,col);
	}
#endif
	return 0;
}

/* ============================================================================= */
int ChannelMap::GetPadRowNCol(int pChan, int &row, int &col)
{
	if(pChan<0 || pChan>=NUM_PADS) return -1;
	return GetPadRowNCol(fZPhiMap[pChan][0],fZPhiMap[pChan][1],row,col);
}

/* ============================================================================= */
//return the vertex Z position, in unit of mm
double ChannelMap::GetZ(int pChan)
{
	if(pChan<0 || pChan>=NUM_PADS) return -99999.;
	return  fZPhiMap[pChan][0];
}

/* ============================================================================= */
//return the Phi angle in rad
double  ChannelMap::GetPhi(int pChan)
{
	if(pChan<0 || pChan>=NUM_PADS) return -99999.;
	return  fZPhiMap[pChan][1];
}

/* ============================================================================= */
//return the Phi angle in rad and z in mm
void  ChannelMap::GetZPhi(int pChan, double &z_mm, double &phi_rad)
{
	if(pChan<0 || pChan>=NUM_PADS) return;
	z_mm = fZPhiMap[pChan][0];
	phi_rad = fZPhiMap[pChan][1];
}

/* ============================================================================= */
//Get pad_id using global ROW and COL
//In order to speed up, do not call this routine, use this buffer 'fChanMap[pad_row][pad_col]'
int ChannelMap::CalculatePadIDByRowNCol(int pad_row, int pad_col)
{	
	if(pad_col<0 || pad_col>=NUM_COLS || pad_row<0 || pad_row>=NUM_ROWS)
	{
		printf("***Error!! CalculatePadIDByRowNCol(row=%3d, col=%3d): Input out of RTPC range\n",
     pad_row, pad_col);
		return -1;
	}
  
  //row and col index of the connector
  int conn_row = pad_row/4;
  int conn_col = pad_col/4;
  
  //the row and col in a connector
  int row_seq = pad_row%4;
  int col_seq = pad_col%4;
  
	int pad_id = conn_row*384 + conn_col*16 + kConnChanMap[row_seq][col_seq];
	if(pad_id<0 || pad_id>=NUM_PADS)
	{ 
		printf("**Error!!! CalculatePadIDByRowNCol(row=%3d, col=%3d) return: pad_id=%5d\n",
			pad_row, pad_col, pad_id);        
		pad_id=-1;
	}
#ifdef CHANNELMAP_DUBUG
	if(CHANNELMAP_DUBUG>=2)
	{
		printf("CalculatePadIDByRowNCol(row=%3d, col=%3d) return: pad_id=%5d\n",
			pad_row, pad_col, pad_id);    
	}
#endif
	return pad_id;
}

/* ============================================================================= */
//Get pad_id using global ROW and COL
int ChannelMap::GetPadIDByRowNCol(int pad_row, int pad_col)
{	
	if(pad_col<0 || pad_col>=NUM_COLS || pad_row<0 || pad_row>=NUM_ROWS)
	{
		printf("***Error!! GetPadIDByRowNCol(row=%3d, col=%3d): Input out of RTPC range\n",
      pad_row, pad_col);
		return -1;
	}
  return fChanMap[pad_row][pad_col];
}

/* ============================================================================= */
//return pad_id for given (z_mm,phi_rad), negative if invalid
int	ChannelMap::GetPadID(double z_mm, double phi_rad)
{
  //global row and col of this position
	int pad_row=-200, pad_col=-200;
	int ret=GetPadRowNCol(z_mm,phi_rad,pad_row,pad_col);
	if(ret<0) return -1;  //error!!! this given (z,phi) is out of RTPC range
  
	//return fChanMap[pad_row][pad_col];
	return CalculatePadIDByRowNCol(pad_row,pad_col);  //this is slow
}
