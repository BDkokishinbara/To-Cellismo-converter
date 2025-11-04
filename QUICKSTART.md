# クイックスタートガイド

最も簡単にデプロイする方法を3つ紹介します。

## 🐳 方法1: Docker（最も簡単）

**必要なもの**: Docker と Docker Compose

```bash
# 1. Dockerをインストール（まだの場合）
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh

# 2. Docker Composeでアプリを起動
docker-compose up -d

# 3. ブラウザでアクセス
# http://localhost:5000
```

**停止する場合:**
```bash
docker-compose down
```

---

## ☁️ 方法2: Render（無料、最も手軽）

**必要なもの**: GitHubアカウント

### 手順

1. **GitHubにコードをプッシュ**
   ```bash
   git init
   git add .
   git commit -m "Initial commit"
   git remote add origin https://github.com/your-username/csv-h5mu.git
   git push -u origin main
   ```

2. **Renderにサインアップ**
   - https://render.com にアクセス
   - GitHubでログイン

3. **Web Serviceを作成**
   - "New" → "Web Service"
   - GitHubリポジトリを選択
   - 以下を入力:
     - **Name**: `singlecell-converter`
     - **Environment**: `Python 3`
     - **Build Command**: `pip install -r requirements.txt`
     - **Start Command**: `gunicorn app:app --bind 0.0.0.0:$PORT --workers 2 --timeout 300`

4. **デプロイをクリック**
   - 数分待つ
   - URLが発行されます（例: https://singlecell-converter.onrender.com）

### 注意点
- 無料プランは15分間使用がないとスリープします
- 初回アクセス時は起動に数秒かかります

---

## 🚂 方法3: Railway（簡単、無料枠あり）

**必要なもの**: GitHubアカウント

### 手順

1. **GitHubにコードをプッシュ**（上記参照）

2. **Railwayにアクセス**
   - https://railway.app
   - GitHubでログイン

3. **デプロイ**
   - "New Project" → "Deploy from GitHub repo"
   - リポジトリを選択
   - 自動的にデプロイが開始

4. **ドメインを生成**
   - Settings → "Generate Domain"
   - 生成されたURLでアクセス

### 料金
- 無料: $5相当/月
- その後は従量課金

---

## 🖥️ 方法4: VPS/自分のサーバー

**必要なもの**: Ubuntu 20.04以上のサーバー、SSHアクセス

### ワンコマンドデプロイ

```bash
# サーバーにSSH接続
ssh user@your-server-ip

# リポジトリをクローン
git clone https://github.com/your-username/csv-h5mu.git
cd csv-h5mu

# 自動デプロイスクリプトを実行
chmod +x deploy.sh
./deploy.sh
```

スクリプトが以下を自動で行います:
- 必要なパッケージのインストール
- Python環境のセットアップ
- Nginxの設定
- Systemdサービスの作成
- ファイアウォールの設定

完了後、サーバーのIPアドレスまたはドメインでアクセス可能になります。

---

## 🔒 SSL証明書の追加（推奨）

### VPS/自分のサーバーの場合

```bash
# Certbotをインストール
sudo apt install -y certbot python3-certbot-nginx

# SSL証明書を取得（your-domain.comを実際のドメインに変更）
sudo certbot --nginx -d your-domain.com

# 自動更新の確認
sudo certbot renew --dry-run
```

### Render/Railway
- 自動的にSSL証明書が発行されます（設定不要）

---

## 📊 デプロイ後の確認

### アプリケーションが動作しているか確認

#### Docker
```bash
docker ps
docker logs singlecell-converter
```

#### VPS/Systemd
```bash
sudo systemctl status singlecell-converter
sudo journalctl -u singlecell-converter -f
```

### テストファイルでの動作確認

1. ブラウザでアプリにアクセス
2. `sample_data.csv`をアップロード
3. 「CSV」を選択
4. 「h5muに変換」をクリック
5. ダウンロードが開始されればOK

---

## ⚙️ 環境変数の設定（オプション）

### Docker
`docker-compose.yml`を編集:
```yaml
environment:
  - SECRET_KEY=your-random-secret-key
  - MAX_CONTENT_LENGTH=524288000
```

### Render/Railway
ダッシュボードから環境変数を追加:
- `SECRET_KEY`: ランダムな文字列
- `FLASK_ENV`: `production`
- `MAX_CONTENT_LENGTH`: `524288000`

### VPS
`.env`ファイルを編集:
```bash
nano .env
```

---

## 🆘 トラブルシューティング

### アプリケーションが起動しない

```bash
# ログを確認
docker logs singlecell-converter  # Docker
sudo journalctl -u singlecell-converter  # VPS
```

### アップロードが失敗する

- ファイルサイズ制限を確認（デフォルト: 500MB）
- `MAX_CONTENT_LENGTH`環境変数を増やす

### メモリ不足エラー

```bash
# Swapを追加（VPS）
sudo fallocate -l 2G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
```

---

## 📚 さらに詳しく

詳細なデプロイオプション、設定、トラブルシューティングについては:
- [DEPLOYMENT.md](DEPLOYMENT.md) - 完全なデプロイガイド
- [README.md](README.md) - アプリケーション概要とローカル開発

---

## 🎉 デプロイ成功！

アプリケーションが正常に動作していれば、以下が可能になります:
- CSVファイルをh5mu形式に変換
- 10x Genomics MEXファイルの変換
- RDSファイルの変換（Rがインストールされている場合）

問題が発生した場合は、[Issues](https://github.com/your-username/csv-h5mu/issues)で報告してください。
