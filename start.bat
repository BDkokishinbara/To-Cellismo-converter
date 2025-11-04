@echo off
chcp 65001 >nul
REM To-Cellismo Converter - 起動スクリプト (Windows)

echo =================================
echo 🧬 To-Cellismo Converter
echo =================================
echo.
echo アプリケーションを起動しています...
echo.

REM Change to script directory
cd /d "%~dp0"

REM Check if virtual environment exists
if exist "venv\Scripts\activate.bat" (
    echo ✓ 仮想環境を検出しました
    call venv\Scripts\activate.bat
)

REM Check if required packages are installed
python -c "import flask" 2>nul
if errorlevel 1 (
    echo 必要なパッケージをインストールしています...
    pip install -r requirements.txt
)

REM Start the Flask app
echo ✓ サーバーを起動しています...
echo.
start /B python app.py

REM Wait for server to start
timeout /t 3 /nobreak >nul

REM Open browser
echo ✓ ブラウザを開いています...
start http://127.0.0.1:5000

echo.
echo =================================
echo ✅ アプリケーションが起動しました！
echo =================================
echo.
echo URL: http://127.0.0.1:5000
echo.
echo 終了するには、このウィンドウを閉じてください
echo.
pause
