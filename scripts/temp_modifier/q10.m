function [tfact] = q10(q10const, tsoil)
tref = 25;
tfact = q10const^((tsoil-tref)/10);
