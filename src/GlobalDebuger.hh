// this is just a header file to turn on or off 
// the debug information

#ifndef _GLOBALDEBUGER_H_
#define _GLOBALDEBUGER_H_ 1

#include <iostream>

//#define STOP4DEBUG Stop4Debug();

extern int  Global_Debug_Level;
extern int  iKeepGoing;
extern int  iPause;
extern int  Global_Skip_Counter;

//in gcc 4.4.6, this routine return an error of
//"expected unqualified-id before numeric constant"
template <typename T>
void ECHO(T i) 
{ 	
	std::cout<<"GlobalDebuger:  debug position "<<i<<std::endl; 
}

template <typename T> 
void ECHO(const char *name, T value)
{
    std::cout<<"GlobalDebuger:  "<<name<<" = "<<value<<std::endl;
}

template <typename T>
void SetGlobalDebugLevel(const char *caller, T val)
{
    Global_Debug_Level=int(val);
    std::cout<<">>> "<<caller<<" set Global_Debug_Level to "<<Global_Debug_Level<<" <<<"<<std::endl;
}

int Pause4Debug();
int Stop4Debug(int iNewEvent=0);
void PrintDebugMenu();

#endif //_GLOBALDEBUGER_H_
