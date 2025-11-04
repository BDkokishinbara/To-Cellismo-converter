# Render デプロイ手順書

このガイドに従って、アプリケーションをRenderに無料でデプロイできます。

## 📋 前提条件

- GitHubアカウント
- Gitがインストールされていること

## 🚀 ステップ1: GitHubにコードをプッシュ

### 1-1. GitHubで新しいリポジトリを作成

1. https://github.com にアクセスしてログイン
2. 右上の「+」ボタン → 「New repository」をクリック
3. リポジトリ情報を入力:
   - **Repository name**: `csv-h5mu-converter` (任意の名前でOK)
   - **Description**: `Single-cell file converter (CSV/RDS/MEX to h5mu)`
   - **Public** または **Private** を選択
   - **Initialize this repository with a README** のチェックは**外す**
4. 「Create repository」をクリック

### 1-2. ローカルからGitHubにプッシュ

リポジトリ作成後に表示されるURLをコピーして、以下のコマンドを実行:

```bash
# プロジェクトディレクトリに移動
cd /home/bdkoki/デスクトップ/csv_h5mu

# GitHubリポジトリをリモートとして追加（URLは自分のものに変更）
git remote add origin https://github.com/YOUR_USERNAME/csv-h5mu-converter.git

# mainブランチにプッシュ
git push -u origin main
```

**注意**: `YOUR_USERNAME` を自分のGitHubユーザー名に変更してください。

GitHubの認証が必要な場合:
- ユーザー名: GitHubのユーザー名
- パスワード: Personal Access Token（パスワードではなくトークンが必要）

**Personal Access Tokenの作成方法**:
1. GitHub → Settings → Developer settings → Personal access tokens → Tokens (classic)
2. "Generate new token" → "Generate new token (classic)"
3. Noteに「Render Deploy」など入力
4. `repo` にチェック
5. "Generate token" をクリック
6. 表示されたトークンをコピー（一度しか表示されません）

---

## 🌐 ステップ2: Renderでデプロイ

### 2-1. Renderアカウントを作成

1. https://render.com にアクセス
2. 「Get Started for Free」または「Sign Up」をクリック
3. **GitHubでサインアップ**を選択（推奨）
4. GitHubでRenderを認証

### 2-2. Web Serviceを作成

1. Renderダッシュボードで「New +」ボタンをクリック
2. 「Web Service」を選択

### 2-3. GitHubリポジトリを接続

1. 「Connect a repository」画面で、先ほど作成したリポジトリを探す
2. 見つからない場合は「Configure account」をクリックして、Renderにリポジトリへのアクセスを許可
3. リポジトリを選択して「Connect」をクリック

### 2-4. サービス設定

以下の情報を入力:

#### Basic Settings
- **Name**: `singlecell-converter` (任意、URLになります)
- **Region**: お好みのリージョン（例: Oregon (US West)、Singapore）
- **Branch**: `main`
- **Root Directory**: (空白のまま)
- **Runtime**: `Python 3`

#### Build & Deploy
- **Build Command**:
  ```
  pip install -r requirements.txt
  ```

- **Start Command**:
  ```
  gunicorn app:app --bind 0.0.0.0:$PORT --workers 2 --timeout 300
  ```

#### Instance Type
- **Free** を選択（無料プラン）

### 2-5. 環境変数を設定（オプション）

「Environment」セクションで「Add Environment Variable」をクリック:

```
SECRET_KEY=ランダムな長い文字列（例: abcd1234efgh5678ijkl9012mnop3456）
FLASK_ENV=production
MAX_CONTENT_LENGTH=524288000
```

**SECRET_KEYの生成方法**:
```bash
python -c "import secrets; print(secrets.token_hex(32))"
```

### 2-6. デプロイを開始

1. 「Create Web Service」ボタンをクリック
2. デプロイが自動的に開始されます

### 2-7. デプロイの進行状況を確認

- 「Logs」タブでビルドとデプロイのログを確認できます
- 初回デプロイは5〜10分かかることがあります
- 「Deploy successful」のメッセージが表示されればOK

---

## ✅ ステップ3: アプリケーションにアクセス

### 3-1. URLを確認

1. Renderのダッシュボードでサービスを選択
2. 上部に表示されるURL（例: `https://singlecell-converter.onrender.com`）をクリック

### 3-2. 動作確認

1. ブラウザでアプリケーションが開きます
2. `sample_data.csv` をアップロードしてテスト:
   - 「CSV」を選択
   - 「マトリックスを転置」にチェック
   - 「h5muに変換」をクリック
3. ファイルがダウンロードされれば成功！

---

## 📝 重要な注意事項

### 無料プランの制限

- **スリープ**: 15分間アクセスがないとアプリがスリープ状態になります
- **起動時間**: スリープから復帰する際、初回アクセスに30秒〜1分かかることがあります
- **メモリ**: 512MB RAM
- **実行時間**: リクエストタイムアウト（大きなファイルの変換に注意）

### 推奨事項

- テスト用途には無料プランで十分
- 本番環境で常時稼働させる場合は有料プラン（$7/月〜）を検討

---

## 🔄 コードの更新方法

コードを変更した後、Renderに再デプロイする方法:

```bash
# 変更をコミット
git add .
git commit -m "Update: 変更内容の説明"

# GitHubにプッシュ
git push origin main
```

Renderは自動的に新しいコミットを検出して再デプロイします。

---

## 🐛 トラブルシューティング

### デプロイが失敗する

**症状**: ビルドやデプロイが失敗する

**対処法**:
1. Logsタブでエラーメッセージを確認
2. `requirements.txt` のパッケージバージョンを確認
3. Start Commandが正しいか確認:
   ```
   gunicorn app:app --bind 0.0.0.0:$PORT --workers 2 --timeout 300
   ```

### アプリケーションが起動しない

**症状**: デプロイは成功するが、URLにアクセスできない

**対処法**:
1. Logsで起動ログを確認
2. 環境変数が正しく設定されているか確認
3. Renderのステータスページを確認: https://status.render.com/

### ファイルアップロードが失敗する

**症状**: 大きなファイルのアップロードが失敗

**対処法**:
1. 無料プランのメモリ制限（512MB）を確認
2. より小さいファイルでテスト
3. 必要に応じて有料プランにアップグレード

### アプリがスリープから復帰しない

**症状**: 久しぶりにアクセスすると応答がない

**対処法**:
1. 1〜2分待つ（スリープからの復帰には時間がかかります）
2. ページをリロード
3. 常時稼働させたい場合は有料プランを検討

---

## 💰 料金について

### 無料プラン
- 750時間/月の無料実行時間
- 複数のサービスで共有
- スリープ機能あり

### 有料プラン（Starter: $7/月）
- スリープなし
- より多いリソース
- カスタムドメイン対応

詳細: https://render.com/pricing

---

## 🎉 完了！

これで、あなたのシングルセル解析ファイル変換アプリケーションがインターネット上で公開されました！

URLを友人や同僚と共有して使ってもらいましょう。

---

## 📚 次のステップ

- カスタムドメインの設定
- SSL証明書の設定（Renderが自動で設定）
- モニタリングとログの確認
- パフォーマンスの最適化

詳細なドキュメント: [DEPLOYMENT.md](DEPLOYMENT.md)
