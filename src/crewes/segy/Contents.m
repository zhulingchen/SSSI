% CREWES SEGY toolbox
% Tools to perform SEG-Y I/O
%
%
%Basic tools
% ALTREADSEGY ... Read SEG-Y files
% ALTWRITESEGY ... Write SEG-Y files
% READSEGY ... Read SEG-Y files (old version, bugs, not recommended)
% WRITESEGY ... Write SEG-Y files (old version, bugs, not recommended)
%
%SEGY Toolbox ... a collection of higher level tools, 
%       type: >> help segy_toolbox 
%       for more information
%
% SEGY_* A set of functions for manipulating rev 1 (2002) SEG-Y data.
% 
% Written mostly by Chad Hogan, bug him for details.
% cmhogan@ucalgary.ca
%
% SEGY_OpenFile ... Opens a SEG-Y formatted disk file, MUST BE USED!
% ^^^^^^ USE THIS FIRST ^^^^^^
%
% SEGY_FindCMPs ... Finds all CMP locations in a SEG-Y file
% SEGY_FindIndex ... Finds all [trace header value] in a SEG-Y file
% SEGY_FindShots ... Finds all shot locations in a SEG-Y file
% SEGY_ModifyTextHeaderLine .. Changes a text header line in-place
% SEGY_ReadBinaryHeader ... Reads the binary header from a SEG-Y file
% SEGY_ReadCMPGather ... Reads a CMP gather from a SEG-Y file
% SEGY_ReadShotGather ... Reads a shot gather from a SEG-Y file
% SEGY_ReadTextHeader ... Reads a text header from a SEG-Y file
% SEGY_ReadTrace ... Reads a trace from a SEG-Y file
% SEGY_ReleaseFile ... Closes a SEG-Y file after you're done with it
% SEGY_ReplaceShotGather ... Replaces a shot gather in a SEG-Y file
% SEGY_ReplaceTrace ... Replace a single trace in a SEG-Y file
%
% FUNCTIONS DESIGNED TO CREATE A NEW SEG-Y FORMATTED FILE
%
% SEGY_GetBinaryHeader ... Returns a binary header template
% SEGY_GetTextHeader ... Returns a text header template
% SEGY_GetTrace ... Returns a complete trace template, header and data.
%
% SEGY_WriteTextHeader ... Write a text header into a SEG-Y file
% SEGY_WriteBinaryHeader ... Writes a binary header into a SEG-Y file
% SEGY_WriteTrace ... Appends a trace (header and data) onto a SEG-Y file
% SEGY_WriteGathers ... example code
% SEGY_WriteStack ... example code
%
% UTILITY FUNCTIONS YOU PROBABLY WON'T USE MUCH
%
% SEGY_TraceSeek ... Seeks the file pointer to a trace in a SEG-Y file
%
% $Id: Contents.m,v 1.4 2008/03/05 00:48:00 cmhogan Exp $