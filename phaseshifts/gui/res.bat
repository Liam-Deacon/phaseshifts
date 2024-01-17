@echo OFF
REM Update resource file
set RES="%DROPBOX%\Programming\Python\LEED-PyV\phaseshifts\gui"
echo Updating resource file: %RES%\res\res.qrc...
set COMPILER=pyrcc4
cd %RES%
%COMPILER% Res\res.qrc -o res_rc.py
echo Updated 'res_rc.py'
pause
