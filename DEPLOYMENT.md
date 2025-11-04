# デプロイガイド

このドキュメントでは、シングルセル解析ファイル変換Webアプリケーションを様々な環境にデプロイする方法を説明します。

## 目次

1. [Docker でのデプロイ](#docker-でのデプロイ)
2. [Render でのデプロイ](#render-でのデプロイ)
3. [Railway でのデプロイ](#railway-でのデプロイ)
4. [VPS/クラウドサーバーでのデプロイ](#vpsクラウドサーバーでのデプロイ)
5. [環境変数の設定](#環境変数の設定)

---

## Docker でのデプロイ

### 前提条件
- Docker がインストールされていること
- Docker Compose がインストールされていること（オプション）

### 方法1: Docker Compose（推奨）

```bash
# リポジトリのクローンまたはダウンロード
cd csv_h5mu

# Docker Composeでビルドして起動
docker-compose up -d

# ログを確認
docker-compose logs -f

# 停止
docker-compose down
```

アプリケーションは `http://localhost:5000` でアクセス可能になります。

### 方法2: Docker コマンド

```bash
# イメージのビルド
docker build -t singlecell-converter .

# コンテナの起動
docker run -d \
  --name singlecell-converter \
  -p 5000:5000 \
  -v $(pwd)/uploads:/app/uploads \
  -v $(pwd)/outputs:/app/outputs \
  singlecell-converter

# ログを確認
docker logs -f singlecell-converter

# 停止
docker stop singlecell-converter
docker rm singlecell-converter
```

### カスタムポートで起動

```bash
docker run -d \
  --name singlecell-converter \
  -p 8080:5000 \
  -e PORT=5000 \
  singlecell-converter
```

`http://localhost:8080` でアクセス可能になります。

---

## Render でのデプロイ

[Render](https://render.com/) は無料枠があり、簡単にデプロイできるクラウドプラットフォームです。

### 手順

1. **GitHubにプッシュ**
   ```bash
   git init
   git add .
   git commit -m "Initial commit"
   git remote add origin <your-github-repo-url>
   git push -u origin main
   ```

2. **Render アカウント作成**
   - [Render](https://render.com/) にアクセスしてアカウントを作成

3. **新しいWeb Serviceを作成**
   - Dashboard → "New" → "Web Service"
   - GitHubリポジトリを接続
   - 以下の設定を入力：

   ```
   Name: singlecell-converter
   Region: お好みのリージョン
   Branch: main
   Root Directory: (空白)
   Environment: Python 3
   Build Command: pip install -r requirements.txt
   Start Command: gunicorn app:app --bind 0.0.0.0:$PORT --workers 2 --timeout 300
   ```

4. **環境変数を設定**（オプション）
   - Environment → "Add Environment Variable"
   ```
   SECRET_KEY=<ランダムな文字列>
   FLASK_ENV=production
   MAX_CONTENT_LENGTH=524288000
   ```

5. **デプロイ**
   - "Create Web Service" をクリック
   - 数分後にデプロイ完了

6. **アクセス**
   - Renderが提供するURL（例: `https://singlecell-converter.onrender.com`）でアクセス

### 注意事項

- **無料プラン**の場合、15分間アクセスがないとスリープ状態になります
- 大きなファイルのアップロードには時間がかかる場合があります
- 無料プランのディスク容量は一時的なもので、再起動時にリセットされます

---

## Railway でのデプロイ

[Railway](https://railway.app/) もGitHubから簡単にデプロイできるプラットフォームです。

### 手順

1. **GitHubにプッシュ**（上記参照）

2. **Railway アカウント作成**
   - [Railway](https://railway.app/) にアクセスしてGitHubでログイン

3. **新しいプロジェクトを作成**
   - "New Project" → "Deploy from GitHub repo"
   - リポジトリを選択

4. **自動検出**
   - Railwayが自動的にPythonアプリケーションを検出
   - `Procfile` を使用して起動コマンドを実行

5. **環境変数を設定**（オプション）
   - Variables タブで環境変数を追加
   ```
   SECRET_KEY=<ランダムな文字列>
   FLASK_ENV=production
   ```

6. **ドメインを生成**
   - Settings → "Generate Domain"
   - 生成されたURLでアクセス

### 料金

- 無料枠: $5相当のクレジット/月
- その後は従量課金

---

## VPS/クラウドサーバーでのデプロイ

### 前提条件
- Ubuntu 20.04/22.04 またはDebian 11/12
- SSHアクセス可能
- sudo権限

### 手順

#### 1. サーバーのセットアップ

```bash
# システムのアップデート
sudo apt update && sudo apt upgrade -y

# 必要なパッケージのインストール
sudo apt install -y python3 python3-pip python3-venv git nginx
```

#### 2. アプリケーションのデプロイ

```bash
# アプリケーション用ユーザーの作成（オプション）
sudo useradd -m -s /bin/bash webapp
sudo su - webapp

# リポジトリのクローン
git clone <your-repo-url> /home/webapp/csv_h5mu
cd /home/webapp/csv_h5mu

# 仮想環境の作成
python3 -m venv venv
source venv/bin/activate

# 依存パッケージのインストール
pip install -r requirements.txt

# 環境変数の設定
cp .env.example .env
nano .env  # 編集
```

#### 3. Systemdサービスの作成

```bash
# 元のユーザーに戻る
exit

# サービスファイルの作成
sudo nano /etc/systemd/system/singlecell-converter.service
```

以下の内容を入力：

```ini
[Unit]
Description=Single-cell File Converter Web Application
After=network.target

[Service]
Type=simple
User=webapp
WorkingDirectory=/home/webapp/csv_h5mu
Environment="PATH=/home/webapp/csv_h5mu/venv/bin"
ExecStart=/home/webapp/csv_h5mu/venv/bin/gunicorn \
    --bind 127.0.0.1:5000 \
    --workers 2 \
    --timeout 300 \
    --access-logfile /var/log/singlecell-converter/access.log \
    --error-logfile /var/log/singlecell-converter/error.log \
    app:app
Restart=always

[Install]
WantedBy=multi-user.target
```

```bash
# ログディレクトリの作成
sudo mkdir -p /var/log/singlecell-converter
sudo chown webapp:webapp /var/log/singlecell-converter

# サービスの有効化と起動
sudo systemctl daemon-reload
sudo systemctl enable singlecell-converter
sudo systemctl start singlecell-converter

# ステータス確認
sudo systemctl status singlecell-converter
```

#### 4. Nginxの設定

```bash
sudo nano /etc/nginx/sites-available/singlecell-converter
```

以下の内容を入力：

```nginx
server {
    listen 80;
    server_name your-domain.com;  # ドメイン名に変更

    client_max_body_size 500M;

    location / {
        proxy_pass http://127.0.0.1:5000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;

        # タイムアウト設定（大きなファイルの変換用）
        proxy_read_timeout 600s;
        proxy_connect_timeout 600s;
        proxy_send_timeout 600s;
    }
}
```

```bash
# シンボリックリンクの作成
sudo ln -s /etc/nginx/sites-available/singlecell-converter /etc/nginx/sites-enabled/

# デフォルトサイトの無効化（オプション）
sudo rm /etc/nginx/sites-enabled/default

# 設定のテスト
sudo nginx -t

# Nginxの再起動
sudo systemctl restart nginx
```

#### 5. SSL証明書の設定（Let's Encrypt）

```bash
# Certbotのインストール
sudo apt install -y certbot python3-certbot-nginx

# SSL証明書の取得
sudo certbot --nginx -d your-domain.com

# 自動更新のテスト
sudo certbot renew --dry-run
```

#### 6. ファイアウォールの設定

```bash
# UFWのインストールと設定
sudo apt install -y ufw
sudo ufw allow 22/tcp
sudo ufw allow 80/tcp
sudo ufw allow 443/tcp
sudo ufw enable
```

---

## 環境変数の設定

アプリケーションは以下の環境変数をサポートしています：

| 変数名 | 説明 | デフォルト値 |
|--------|------|-------------|
| `FLASK_ENV` | 実行環境（development/production） | `development` |
| `SECRET_KEY` | Flaskのシークレットキー | `dev-secret-key-change-in-production` |
| `PORT` | アプリケーションのポート | `5000` |
| `HOST` | バインドするホスト | `0.0.0.0` |
| `MAX_CONTENT_LENGTH` | 最大アップロードサイズ（バイト） | `524288000` (500MB) |
| `UPLOAD_FOLDER` | アップロードディレクトリ | `uploads` |
| `OUTPUT_FOLDER` | 出力ディレクトリ | `outputs` |

### 安全なSECRET_KEYの生成

```python
python -c "import secrets; print(secrets.token_hex(32))"
```

---

## トラブルシューティング

### メモリ不足エラー

大きなファイルを変換する際にメモリ不足が発生する場合：

```bash
# サーバーのメモリを確認
free -h

# Swapの追加（2GB）
sudo fallocate -l 2G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab
```

### アップロードサイズ制限

Nginxを使用している場合、`client_max_body_size`を設定してください：

```nginx
server {
    # ...
    client_max_body_size 500M;
    # ...
}
```

### ログの確認

```bash
# Systemdサービスのログ
sudo journalctl -u singlecell-converter -f

# アプリケーションログ
sudo tail -f /var/log/singlecell-converter/error.log

# Nginxログ
sudo tail -f /var/log/nginx/error.log
```

---

## パフォーマンスチューニング

### Gunicornワーカー数の調整

推奨: `(2 x CPU cores) + 1`

```bash
# 4コアの場合
--workers 9
```

### タイムアウトの調整

大きなファイルの変換には時間がかかるため、タイムアウトを調整：

```bash
--timeout 600  # 10分
```

---

## セキュリティ推奨事項

1. **SECRET_KEYを必ず変更**
2. **HTTPSを使用**（Let's Encryptで無料）
3. **ファイアウォールを設定**
4. **定期的なアップデート**
   ```bash
   sudo apt update && sudo apt upgrade -y
   ```
5. **アップロードファイルの検証**（既に実装済み）
6. **レート制限の実装**（オプション、Flask-Limiter使用）

---

## バックアップ

重要なデータのバックアップ：

```bash
# 設定ファイル
cp .env .env.backup

# データベース（将来的に追加する場合）
# mysqldump -u user -p database > backup.sql
```

---

## 監視とログ

### アプリケーションの監視

```bash
# サービスの状態
sudo systemctl status singlecell-converter

# リソース使用状況
htop  # または top

# ディスク使用状況
df -h
du -sh uploads/ outputs/
```

### ログローテーション

```bash
sudo nano /etc/logrotate.d/singlecell-converter
```

```
/var/log/singlecell-converter/*.log {
    daily
    rotate 14
    compress
    delaycompress
    notifempty
    missingok
    create 0640 webapp webapp
}
```

---

## サポート

問題が発生した場合は、GitHubのIssuesセクションで報告してください。
