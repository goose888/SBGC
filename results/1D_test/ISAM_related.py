# -*- Python source file -*-
"""
File Name : ISAM_related.py

Description : Test code related to the ISAM model calculations in isamcalc_lib

Created on : Mon Jun 18 23:33:19 2018

Last Modified : Mon Jun 18 23:34:17 2018

Author : Shijie Shu
"""

import isamcalc_lib as isam
z, dz, zsoih = isam.get_isam_soildp(10)
mod = isam.get_depth_modifier(10, z)
mod
zsoih
