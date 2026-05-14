#!/bin/bash
# To-Cellismo Converter - 起動スクリプト
# このスクリプトをダブルクリックするだけでアプリが立ち上がります。
# 必要な Python パッケージはプロジェクト内 .venv/ にだけインストールされるため、
# ユーザーのシステム Python は一切汚しません。

set -e

echo "================================="
echo "🧬 To-Cellismo Converter"
echo "================================="
echo ""

# スクリプトのあるディレクトリへ移動
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$SCRIPT_DIR"

VENV_DIR=".cellismo_venv"

# Python 実行コマンドを検出
# 安定版 (3.11/3.12/3.13) を優先し、無ければ python3 / python を使う
PY=""
for candidate in python3.13 python3.12 python3.11 python3 python; do
    if command -v "$candidate" > /dev/null 2>&1; then
        PY="$candidate"
        break
    fi
done

if [ -z "$PY" ]; then
    echo "❌ Python が見つかりません。"
    echo "   https://www.python.org/downloads/ から Python 3.10〜3.13 をインストールしてください。"
    read -n 1 -s -r -p "何かキーを押すと終了します..."
    exit 1
fi

echo "✓ 使用する Python: $($PY --version) ($PY)"

# 仮想環境を準備（無ければ作成）
if [ ! -d "$VENV_DIR" ]; then
    echo "📦 専用の仮想環境を作成しています（初回のみ・数十秒かかります）..."
    "$PY" -m venv "$VENV_DIR"
    echo "✓ 仮想環境を作成しました: $VENV_DIR/"
fi

# 仮想環境内の python を直接使う（activate しなくて良い）
VENV_PY="$VENV_DIR/bin/python"

# 依存パッケージが入っているか確認
if ! "$VENV_PY" -c "import flask, pandas, anndata, mudata" > /dev/null 2>&1; then
    echo "📥 必要なパッケージを仮想環境にインストールしています（初回のみ・数分かかります）..."

    # インストール失敗時に壊れた仮想環境を消すための後始末
    cleanup_on_failure() {
        echo ""
        echo "❌ パッケージのインストールに失敗しました。"
        echo "   不完全な仮想環境フォルダを削除します: $VENV_DIR/"
        rm -rf "$VENV_DIR"
        echo ""
        echo "💡 もう一度 start.command (または start.sh) をダブルクリックして再試行してください。"
        echo "   それでも失敗する場合は、お手元の Python バージョンを確認してください："
        echo "   $($PY --version)"
        echo "   Python 3.10〜3.13 を推奨します（最新の 3.14 だと一部パッケージが対応していない場合あり）。"
        read -n 1 -s -r -p "何かキーを押すと終了します..."
        exit 1
    }
    trap cleanup_on_failure ERR

    "$VENV_PY" -m pip install --upgrade pip > /dev/null
    "$VENV_PY" -m pip install -r requirements.txt

    trap - ERR
    echo "✓ インストール完了"
fi

# サーバー起動
echo ""
echo "✓ サーバーを起動しています..."
"$VENV_PY" app.py &
APP_PID=$!

# 起動を待ってブラウザを開く
sleep 3
echo "✓ ブラウザを開いています..."
if command -v open > /dev/null 2>&1; then
    open http://127.0.0.1:5050
elif command -v xdg-open > /dev/null 2>&1; then
    xdg-open http://127.0.0.1:5050
elif command -v gnome-open > /dev/null 2>&1; then
    gnome-open http://127.0.0.1:5050
fi

echo ""
echo "================================="
echo "✅ アプリケーションが起動しました！"
echo "================================="
echo ""
echo "URL: http://127.0.0.1:5050"
echo ""
echo "終了するには Ctrl+C を押してください"
echo ""

# Ctrl+C で子プロセスを終了
trap 'kill $APP_PID 2>/dev/null; exit 0' INT TERM
wait $APP_PID
