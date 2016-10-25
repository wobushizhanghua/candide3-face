//
//  CDApp.cpp
//  candide
//
//  Created by Damian Stewart on 03.11.13.
//  Copyright (c) 2013 Damian Stewart. All rights reserved.
//

#include "CDApp.h"
#include "CDWindow.h"

#include <Fl/Fl.h>
#include <assert.h>


static CDApp* instance = NULL;

CDApp* CDApp::getInstance()
{
	assert(instance);
	return instance;
}

CDApp::CDApp( int argc, const char* argv[] )
{
	instance = this;
	this->argc = argc;
	this->argv = argv;
}


static void exitCallback( Fl_Widget* widget )
{
	// disable escape key
	if (Fl::event()==FL_SHORTCUT && Fl::event_key()==FL_Escape) {
		// do nothing
	} else {
		exit(0);
	}
}

int CDApp::run(int (*align_callback)(int, const char**, CDFaceWindow* fw))
{
	CDWindow window(450,500,"Candide");
	
	if (align_callback)
	{
		align_callback( argc, argv, (window.faceWindow) );
	}
	
	window.callback(&exitCallback);
	
	window.end();
	window.show(0,NULL);
	
	return Fl::run();
}


