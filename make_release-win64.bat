@echo off
setlocal enabledelayedexpansion

set "PROJECT_NAME=qrustyquake"
set "EXE_NAME=%PROJECT_NAME%.exe"
set "VERSION=0.8.1"
set "ZIP_NAME=%PROJECT_NAME%_%VERSION%-win64.zip"
set "GAME_NAME=id1"
set "BUILD_DIR=%CD%\build\Release"
set "EXPORT_DIR=%CD%\dist"
set "DATA_DIR=C:\Quake\%GAME_NAME%"
set "RESHACKER=%USERPROFILE%\Desktop\ResHacker\ResourceHacker.exe"
set "SDL3=%CD%\lib\SDL\build\Release\SDL3.dll"
set "SDL3_MIXER=%CD%\lib\SDL_mixer\build\Release\SDL3_mixer.dll"
set "DLLS=%SDL3% %SDL3_MIXER%"
set "PAKS=%DATA_DIR%\PAK0.PAK %DATA_DIR%\PAK1.PAK"

if exist "%EXPORT_DIR%" (
echo Cleaning up existing files...
rd /s /q "%EXPORT_DIR%"
)
md "%EXPORT_DIR%" "%EXPORT_DIR%\temp_package" "%EXPORT_DIR%\temp_package\%GAME_NAME%"

%RESHACKER% -open %BUILD_DIR%\%EXE_NAME% -save %EXPORT_DIR%\temp_package\%EXE_NAME% ^
-resource %BUILD_DIR%\icon.ico -mask ICONGROUP,MAINICON, -action addoverwrite -log CONSOLE

echo Gathering libraries and game data...
for %%f in (%DLLS%) do (
    if exist "%%f" (
		echo Fetching: %%f
        copy "%%f" "%EXPORT_DIR%\temp_package\"
    ) else (
        echo WARNING: Missing library %%f!
    )
)

for %%f in (%PAKS%) do (
	if exist "%%f" (
		copy "%%f" "%EXPORT_DIR%\temp_package\%GAME_NAME%"
	)
)

echo Creating ZIP archive...
cd "%EXPORT_DIR%\temp_package"
tar -a -c -f "..\%ZIP_NAME%" *
cd ..\..
rd /s /q "%EXPORT_DIR%\temp_package"
echo.
echo.
echo ------------------------------------------------
echo  Your build is ready: '%EXPORT_DIR%\%ZIP_NAME%'
echo ------------------------------------------------
echo Press any key to exit.
pause >nul
explorer.exe %EXPORT_DIR%
exit /b 0