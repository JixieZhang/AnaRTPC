// this is just a header file to turn on or off 
// the debug information
// #define ECHO(i) printf(">>>debug position %g\n",(float)(i));

#include "GlobalDebuger.hh"
#include <stdlib.h>
#include <iostream>
#include <string>
#include <stdio.h>
using namespace std;

int  Global_Debug_Level=1;
int  iKeepGoing=0;	  // 0 false; 1 true
int  iPause=1;	          // 0 false; 1 true
int  iNoPauseThisEvent=0; // 0 false; 1 true
int  Global_Skip_Counter=0;

void PrintDebugMenu()
{  
  cout<<endl;
  cout<<"==================================================="<<endl;
  cout<<"GlobalDebug: The help menu of global debug system is:"<<endl;
  cout<<"\t h|H:  print help menu and current status"<<endl
    <<"\t l|L #:  change the debug level, if the level==-999, exit immediately"<<endl
    <<"\t j|J #:  process # of events without stopping or showing this message"<<endl
    <<"\t k|K:    keep running without showing this message till the end"<<endl
    <<"\t q|Q:    quit now"<<endl;
  cout<<"---------------------------------------------------"<<endl;
  cout<<"Current debug information: Global_Debug_Level="<<Global_Debug_Level
    <<" iKeepGoing="<<((iKeepGoing)? "true" : "false")<<endl;
  cout<<"==================================================="<<endl;
  cout<<endl;
}	

int Stop4Debug(int iNewEvent)
{
  iNoPauseThisEvent=0;
  if(Global_Skip_Counter>0) {Global_Skip_Counter--;return 0;}
  if(!iKeepGoing) 
  { 
    char foo;
    if(iNewEvent) cout<<"=======================Event End  ============================"<<endl;
    else cout<<"=======================DEBUG BLOCK End  ============================"<<endl;

    cout<<"GlobalDebug: Press k to skip asking, L# change debug level, q to quit, any other key to continue:";
    scanf("%c",&foo);
    if(foo=='k' || foo=='K')  iKeepGoing=1;
    else if(foo=='h' || foo=='H')  PrintDebugMenu();
    else if(foo=='l' || foo=='L')  
    {
      scanf("%d",&Global_Debug_Level);
      if(Global_Debug_Level==-999) exit(Global_Debug_Level);
      cout<<"==================================================="<<endl;
      cout<<"GlobalDebug: Changing global debug level to "<<Global_Debug_Level<<endl;
      cout<<"==================================================="<<endl;
    }
    else if(foo=='j' || foo=='J')  
    {
      scanf("%d",&Global_Skip_Counter); 
      if(Global_Skip_Counter==-999) exit(Global_Skip_Counter);
      cout<<"==================================================="<<endl;
      cout<<"GlobalDebug: skip the next "<<Global_Skip_Counter<<" events"<<endl;
      cout<<"==================================================="<<endl;
      Global_Skip_Counter-=2;  //j# here includes 2 scanf events, I have to subtract 2 
    }
    else if(foo=='q' || foo=='Q')  
    {
      abort();
    }

    if(iNewEvent) cout<<"=======================Event Start ============================"<<endl;
    else cout<<"=======================DEBUG BLOCK Start============================"<<endl;
  }
  return iKeepGoing;
}


int Pause4Debug()
{
  if(iPause && !iNoPauseThisEvent) 
  { 
    char foo;

    cout<<"GlobalDebug: Paused. Press k to disable pause, p to stop asking for this event, q to quit, any other keys to continue:";
    scanf("%c",&foo);
    if(foo=='k' || foo=='K')  iPause=0;
    else if(foo=='p' || foo=='P')  iNoPauseThisEvent=1;
    else if(foo=='h' || foo=='H')  PrintDebugMenu();
    else if(foo=='l' || foo=='L')  
    {
      scanf("%d",&Global_Debug_Level);
      if(Global_Debug_Level==-999) exit(Global_Debug_Level);
      cout<<"==================================================="<<endl;
      cout<<"GlobalDebug: Changing global debug level to "<<Global_Debug_Level<<endl;
      cout<<"==================================================="<<endl;
    }
    else if(foo=='q' || foo=='Q')  
    {
      abort();
    }
  }
  return iPause;
}

