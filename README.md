LITTPlan
========

LITT plan module for Slicer 4

## Build Instructions

 * http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Developers/Build_Instructions
 * https://www.slicer.org/slicerWiki/index.php/Documentation/4.1/Developers/Build_Module

## Mac Version

	$ uname -a
	Darwin Einstein.local 12.5.0 Darwin Kernel Version 12.5.0: Mon Jul 29 16:33:49 PDT 2013; root:xnu-2050.48.11~1/RELEASE_X86_64 x86_64
	$ sw_vers | grep 'ProductVersion:'
	ProductVersion: 10.8.5

## download dependencies

	$ ll SlicerBuildDownloads/
	total 5298848
	-rw-r--r--@ 1 staff  staff   191356132 Sep 13  2012 qt-mac-opensource-4.8.3.dmg
	-rw-r--r--@ 1 staff  staff    67440114 Nov 10 18:40 XQuartz-2.7.5.dmg
	-rw-r--r--@ 1 staff  staff    42478875 Jan 16 13:48 cmake-2.8.12.2-Darwin64-universal.dmg
	-rw-r--r--@ 1 staff  staff   124724212 Mar 10 14:03 commandline_tools_os_x_mountain_lion_for_xcode__march_2014.dmg
	-rw-r--r--@ 1 staff  staff  2286987512 Mar 10 14:21 xcode_5.1.dmg
	$ which qmake
	/usr/bin/qmake
	$ qmake --version
	QMake version 2.01a
	Using Qt version 4.8.3 in /Library/Frameworks
	$ which cmake
	/usr/bin/cmake
	$ cmake --version
	cmake version 2.8.12.2
	$ xcode-select --print-path
	/Applications/Xcode.app/Contents/Developer
	$ xcode-select --version
	xcode-select version 2311.

## configure
	$svn co http://svn.slicer.org/Slicer4/trunk Slicer4
	$ svn info Slicer4
	URL: http://svn.slicer.org/Slicer4/trunk
	Repository Root: http://svn.slicer.org/Slicer4
	Repository UUID: 3bd1e089-480b-0410-8dfb-8563597acbee
	Revision: 23047
	Node Kind: directory
	Schedule: normal
	Last Changed Author: finetjul
	Last Changed Rev: 23047
	Last Changed Date: 2014-04-04 17:34:32 -0500 (Fri, 04 Apr 2014)

        mkdir Slicer4-SuperBuild-Debug/; cd Slicer4-SuperBuild-Debug/
	cmake -DCMAKE_OSX_DEPLOYMENT_TARGET=10.8 -DCMAKE_OSX_SYSROOT=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk -DSlicer_USE_PYTHONQT_WITH_TCL:BOOL=OFF  -DCMAKE_INSTALL_PREFIX=/Users/fuentes/MyProjects/Slicer4-Install ../Slicer4

	$ vim VTK/CMakeLists.txt 
	22,26c22,26
	< ##  IF(OSX_SDK_VERSION)
	< ##    IF(${OSX_SDK_VERSION} VERSION_GREATER "10.4")
	< ##      SET(VTK_OBJCXX_FLAGS_DEFAULT "-fobjc-gc")
	< ##    ENDIF(${OSX_SDK_VERSION} VERSION_GREATER "10.4")
	< ##  ENDIF(OSX_SDK_VERSION)
	---
	>   IF(OSX_SDK_VERSION)
	>     IF(${OSX_SDK_VERSION} VERSION_GREATER "10.4")
	>       SET(VTK_OBJCXX_FLAGS_DEFAULT "-fobjc-gc")
	>     ENDIF(${OSX_SDK_VERSION} VERSION_GREATER "10.4")
	>   ENDIF(OSX_SDK_VERSION)
