LITTPlan
========

LITT plan module for Slicer 4

## Build Instructions

http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Developers/Build_Instructions

svn co http://svn.slicer.org/Slicer4/trunk


## Mac Build

install Xcode from app store

## download dependencies

$ ls /Users/fuentes/Downloads/SlicerBuildDownloads
XQuartz-2.7.5.dmg                                               commandline_tools_os_x_mountain_lion_for_xcode__march_2014.dmg  qt-mac-opensource-4.7.4.dmg
cmake-2.8.12.2-Darwin64-universal.dmg                           qt-mac-opensource-4.7.4-debug-libs.dmg
$ uname -a
Darwin Einstein.local 12.5.0 Darwin Kernel Version 12.5.0: Mon Jul 29 16:33:49 PDT 2013; root:xnu-2050.48.11~1/RELEASE_X86_64 x86_64


## configure
cmake -DCMAKE_OSX_DEPLOYMENT_TARGET=10.8 -DCMAKE_OSX_SYSROOT=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk -DSlicer_USE_PYTHONQT_WITH_TCL:BOOL=OFF ../Slicer4

