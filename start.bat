@echo off
chcp 65001 >nul
setlocal enabledelayedexpansion

REM To-Cellismo Converter - 起動スクリプト (Windows)
REM このファイルをダブルクリックするだけでアプリが立ち上がります。
REM 必要な Python パッケージはプロジェクト内 .venv\ にだけインストールされるため、
REM ユーザーのシステム Python は一切汚しません。

echo =================================
echo 🧬 To-Cellismo Converter
echo =================================
echo.

REM スクリプトのあるフォルダへ移動
cd /d "%~dp0"

set "VENV_DIR=.cellismo_venv"

REM Python 実行コマンドを検出
REM 安定版 (py -3.13/3.12/3.11) を優先、無ければ py / python を使う
set "PY="
for %%v in (3.13 3.12 3.11) do (
    py -%%v --version >nul 2>nul
    if !errorlevel!==0 (
        set "PY=py -%%v"
        goto :py_found
    )
)

where py >nul 2>nul
if %errorlevel%==0 (
    set "PY=py -3"
    goto :py_found
)

where python >nul 2>nul
if %errorlevel%==0 (
    set "PY=python"
    goto :py_found
)

echo ❌ Python が見つかりません。
echo    https://www.python.org/downloads/ から Python 3.10〜3.13 をインストールしてください。
pause
exit /b 1

:py_found
echo ✓ 使用する Python: %PY%

REM 仮想環境を準備（無ければ作成）
if not exist "%VENV_DIR%\Scripts\python.exe" (
    echo 📦 専用の仮想環境を作成しています（初回のみ・数十秒かかります）...
    %PY% -m venv "%VENV_DIR%"
    if errorlevel 1 (
        echo ❌ 仮想環境の作成に失敗しました。
        pause
        exit /b 1
    )
    echo ✓ 仮想環境を作成しました: %VENV_DIR%\
)

set "VENV_PY=%VENV_DIR%\Scripts\python.exe"

REM 依存パッケージが入っているか確認
"%VENV_PY%" -c "import flask, pandas, anndata, mudata" >nul 2>nul
if errorlevel 1 (
    echo 📥 必要なパッケージを仮想環境にインストールしています（初回のみ・数分かかります）...
    "%VENV_PY%" -m pip install --upgrade pip >nul
    "%VENV_PY%" -m pip install -r requirements.txt
    if errorlevel 1 (
        echo.
        echo ❌ パッケージのインストールに失敗しました。
        echo    不完全な仮想環境フォルダを削除します: %VENV_DIR%\
        rmdir /s /q "%VENV_DIR%"
        echo.
        echo 💡 もう一度 start.bat をダブルクリックして再試行してください。
        echo    それでも失敗する場合は、お手元の Python バージョンを確認してください。
        echo    Python 3.10〜3.13 を推奨します（最新の 3.14 だと一部パッケージが対応していない場合あり）。
        pause
        exit /b 1
    )
    echo ✓ インストール完了
)

REM サーバー起動
echo.
echo ✓ サーバーを起動しています...
start /B "" "%VENV_PY%" app.py

REM 起動を待ってブラウザを開く
timeout /t 3 /nobreak >nul
echo ✓ ブラウザを開いています...
start "" http://127.0.0.1:5050

echo.
echo =================================
echo ✅ アプリケーションが起動しました！
echo =================================
echo.
echo URL: http://127.0.0.1:5050
echo.
echo 終了するには、このウィンドウを閉じてください
echo.
pause
