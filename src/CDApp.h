//
//  CDApp.h
//  candide
//
//  Created by Damian Stewart on 03.11.13.
//  Copyright (c) 2013 Damian Stewart. All rights reserved.
//

#ifndef __candide__CandideApp__
#define __candide__CandideApp__

#include <iostream>

#include "CDFaceWindow.h"

class CDApp
{
public:
	CDApp( int argc, const char* argv[] );
	
	static CDApp* getInstance();
	int argc;
	const char** argv;
		
	int run( int (*align_callback)(int argc, const char** argv, CDFaceWindow* fw) =0 );
	
private:
	
};

#endif /* defined(__candide__CandideApp__) */
