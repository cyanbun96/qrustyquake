/*
Copyright (C) 1996-1997 Id Software, Inc.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#define	MAX_HANDLES 32

// File IO
int Sys_FileOpenRead (char *path, int *hndl);
int Sys_FileOpenWrite (char *path);
void Sys_FileClose (int handle);
void Sys_FileSeek (int handle, int position);
int Sys_FileRead (int handle, void *dest, int count);
int Sys_FileWrite (int handle, void *data, int count);
int Sys_FileTime (char *path);
void Sys_mkdir (char *path);

// System IO
void Sys_DebugLog(char *file, char *fmt, ...);
void Sys_Error (char *error, ...);
void Sys_Printf (char *fmt, ...);
void Sys_Quit ();
double Sys_FloatTime ();
