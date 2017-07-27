//////////////////////////////////////////////////////////////////////
// RTPC12 ChannelMap
// Define how the channal id is assigned in location
// Jixie Zhang, v1.0
//////////////////////////////////////////////////////////////////////

#ifndef _CHANNELMAP_
#define _CHANNELMAP_

#define NUM_PADS 17280
#define NUM_ROWS 180
#define NUM_COLS 96


class ChannelMap
{
public:
  //providing pad width(phi), length(z), inner_s, phi at center of the 
  //bottom edge of the connector connected to channel 0
  //ChannelMap* fChanMap = new ChannelMap(2.7381,4.1016,79.00,0);
  //ChannelMap* fChanMap = new ChannelMap(2.8016,4.1016,80.42,0);
	ChannelMap(double pad_w=2.8016, double pad_l=4.1016, double pad_r_in=80.42, 
             double phi_conn_edge=0.);
	virtual ~ChannelMap();

	double GetZ(int chan);    //return the vertex Z position in mm
	double GetPhi(int chan);  //return the Phi angle in rad
	void   GetZPhi(int chan, double &z_mm, double &phi_rad);

	int    GetPadRowNCol(int chan,int &row,int &col);
	int    GetPadRowNCol(double z_mm, double phi_rad, int &row, int &col);
  
	int	   GetPadID(double z_mm,double phi_rad);
	int    GetPadIDByRowNCol(int row, int col);

  void   SetPhiStart(double val){fPhiStart=val;};
  void   SetZStart(double *val){for(int i=0;i<4;i++) fZStart[i]=val[i];};
  
	void   GenerateMap(); 
  
private:
	int    CalculatePadIDByRowNCol(int row, int col);

private:
  // the buffer to store the map: zphi map and pad_id map
	double fZPhiMap[NUM_PADS][2]; //in the order of z_mm:phi_rad
	int    fChanMap[NUM_ROWS][NUM_COLS]; 
  
//some parameter of the Detector  
private:
  double fDeltaA;//=2.7381;     //phi dimension of a pad (the arc length) (2.7mm+2mil gap)
  double fDeltaZ;//=4.1016;     //a single pad length in z (4mm+4mil gap)
  double fPadRi; //=79.0;        //inner radius where the pad plane is located (79mm), 10mil thick
  double fDeltaPhi;//=fDeltaA/fPadRi;  //the angular coverage of a single pad (rad)
 
  //For the above number, fDeltaPhi=1.9950535479061624 degree 
  //the dead area is about 3.5mm in length
  
  //define phi at the bottom edge of the the connector connect to pad_id==0 
  double fPhiStart;// = 0.5*fDeltaPhi; //in rad
  //define left edge z of these 4 rows for the connector connect to piad_id==0 
  double fZStart[4];
  
};
typedef ChannelMap RTPC12ChannelMap;

#endif //_CHANNELMAP_
