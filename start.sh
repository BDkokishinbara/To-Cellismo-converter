#!/bin/bash
# To-Cellismo Converter - 起動スクリプト

echo "================================="
echo "🧬 To-Cellismo Converter"
echo "================================="
echo ""
echo "アプリケーションを起動しています..."
echo ""

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$SCRIPT_DIR"

# Check if virtual environment exists
if [ -d "venv" ]; then
    echo "✓ 仮想環境を検出しました"
    source venv/bin/activate
fi

# Check if required packages are installed
if ! python -c "import flask" &> /dev/null; then
    echo "必要なパッケージをインストールしています..."
    pip install -r requirements.txt
fi

# Start the Flask app in the background
echo "✓ サーバーを起動しています..."
echo ""
python app.py &
APP_PID=$!

# Wait for server to start
sleep 3

# Open browser
echo "✓ ブラウザを開いています..."
if command -v xdg-open > /dev/null; then
    xdg-open http://127.0.0.1:5000
elif command -v gnome-open > /dev/null; then
    gnome-open http://127.0.0.1:5000
elif command -v open > /dev/null; then
    open http://127.0.0.1:5000
fi

echo ""
echo "================================="
echo "✅ アプリケーションが起動しました！"
echo "================================="
echo ""
echo "URL: http://127.0.0.1:5000"
echo ""
echo "終了するには Ctrl+C を押してください"
echo ""

# Wait for the Flask app
wait $APP_PID
