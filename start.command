#!/bin/bash
# To-Cellismo Converter - macOS 用ダブルクリック起動ファイル
#
# macOS では .sh ファイルはエディタ（Cursor / VSCode など）で開いてしまうため、
# 「ダブルクリックで実行」用に .command 拡張子のファイルを用意しています。
# 中身は start.sh を呼び出すだけです。

# スクリプトのあるフォルダへ移動
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$SCRIPT_DIR"

# 本体スクリプトを実行
exec bash ./start.sh
