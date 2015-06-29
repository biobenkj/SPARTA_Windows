REM This command adds python and R to the PATH. Assumes R-3.2.0 is installed based on the version distributed with SPARTA for Windows
setx path "C:\python27;C:\python27\Scripts\;C:\Program Files\R\R-3.2.0\bin\i386;C:\Program Files (x86)\GnuWin32\bin"
REM Setting the local library environment variable to be able to install edgeR
setx R_LIBS_USER %cd%\R_local\library_local