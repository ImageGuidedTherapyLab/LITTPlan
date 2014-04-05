LITTPlan
========

LITT plan module for Slicer 4

## Build Instructions

	http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Developers/Build_Instructions
	https://www.slicer.org/slicerWiki/index.php/Documentation/4.1/Developers/Build_Module

## Mac Build

	$ uname -a
	Darwin Einstein.local 12.5.0 Darwin Kernel Version 12.5.0: Mon Jul 29 16:33:49 PDT 2013; root:xnu-2050.48.11~1/RELEASE_X86_64 x86_64
	$ sw_vers | grep 'ProductVersion:'
	ProductVersion: 10.8.5

	install Xcode from app store

## download dependencies

	$ curl -O http://packages.kitware.com/download/item/3735/qt-mac-opensource-4.8.4.dmg
	$ ls /Users/fuentes/Downloads/SlicerBuildDownloads
	XQuartz-2.7.5.dmg                                               commandline_tools_os_x_mountain_lion_for_xcode__march_2014.dmg
	cmake-2.8.12.2-Darwin64-universal.dmg                           qt-mac-opensource-4.8.4.dmg


## configure
	$svn co http://svn.slicer.org/Slicer4/trunk Slicer4
	$ svn info Slicer4
	URL: http://svn.slicer.org/Slicer4/trunk
	Repository Root: http://svn.slicer.org/Slicer4
	Repository UUID: 3bd1e089-480b-0410-8dfb-8563597acbee
	Revision: 23045
	Node Kind: directory
	Schedule: normal
	Last Changed Author: alexy
	Last Changed Rev: 23045
	Last Changed Date: 2014-04-04 07:59:06 -0500 (Fri, 04 Apr 2014)

        mkdir Slicer4-SuperBuild-Debug/; cd Slicer4-SuperBuild-Debug/
	cmake -DCMAKE_OSX_DEPLOYMENT_TARGET=10.8 -DCMAKE_OSX_SYSROOT=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk -DSlicer_USE_PYTHONQT_WITH_TCL:BOOL=OFF ../Slicer4

