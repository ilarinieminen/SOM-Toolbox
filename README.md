SOM-Toolbox
===========

A Matlab toolbox for Self-Organizing Maps (SOM) and others.

SOM Toolbox 2.0, a software library for Matlab 5 implementing the
Self-Organizing Map algorithm is Copyright (C) 1999 by Esa Alhoniemi,
Johan Himberg, Jukka Parviainen and Juha Vesanto.

SOM Toolbox 2.1 is a revision to the SOM Toolbox made in 
December 2012. See CHANGELOG for changes.

Setup
-----

1. Download
2. addpath(genpath('SOM-Toolbox'));

This will add SOM-Toolbox and all subfolders into matlabpath.

Directories included in the toolbox
-----------------------------------

`som/`      - SOM functionality included in SOM Toolbox 2.0 (revised for this version)

`gtm/`      - GTM functionality eploiting Netlab GTM functions.

`contrib/`  - External contributions to SOM Toolbox

`demo/`     - Demo scripts

`data/`     - Data for demo scripts

`dijkstra/` - Dijkstra algorithm from pmtk3 package (used by som_quality)


COPYRIGHT
---------

This package is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or any later version.

Note: only part of the files distributed in the package belong to the
SOM Toolbox. The package also contains contributed files, which may
have their own copyright notices. If not, the GNU General Public
License holds for them, too, but so that the author(s) of the file
have the Copyright.

This package is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this package (file COPYING); if not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

LICENSING THE LIBRARY FOR PROPRIETARY PROGRAMS
----------------------------------------------

As stated in the GNU General Public License (see the license in COPYING)
it is not possible to include this software library in a commercial
proprietary program without written permission from the owners of the
copyright.
