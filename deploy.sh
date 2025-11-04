#!/bin/bash
# Quick deployment script for VPS/Ubuntu servers

set -e

echo "=========================================="
echo "Single-cell File Converter Deployment"
echo "=========================================="
echo ""

# Check if running on Ubuntu/Debian
if [ ! -f /etc/debian_version ]; then
    echo "This script is designed for Ubuntu/Debian systems."
    exit 1
fi

# Update system
echo "Updating system packages..."
sudo apt update && sudo apt upgrade -y

# Install required packages
echo "Installing required packages..."
sudo apt install -y python3 python3-pip python3-venv git nginx

# Get application directory
read -p "Enter installation directory [/opt/singlecell-converter]: " INSTALL_DIR
INSTALL_DIR=${INSTALL_DIR:-/opt/singlecell-converter}

# Create directory
echo "Creating installation directory: $INSTALL_DIR"
sudo mkdir -p $INSTALL_DIR
sudo chown $USER:$USER $INSTALL_DIR

# Clone or copy application
if [ ! -f "$INSTALL_DIR/app.py" ]; then
    echo "Copying application files..."
    cp -r . $INSTALL_DIR/
fi

cd $INSTALL_DIR

# Create virtual environment
echo "Creating virtual environment..."
python3 -m venv venv
source venv/bin/activate

# Install dependencies
echo "Installing Python dependencies..."
pip install --upgrade pip
pip install -r requirements.txt

# Create .env file
if [ ! -f .env ]; then
    echo "Creating .env file..."
    cp .env.example .env

    # Generate secret key
    SECRET_KEY=$(python3 -c "import secrets; print(secrets.token_hex(32))")
    sed -i "s/your-secret-key-here-change-this/$SECRET_KEY/" .env
    sed -i "s/FLASK_ENV=production/FLASK_ENV=production/" .env

    echo ".env file created. Please review and edit if needed."
fi

# Create necessary directories
mkdir -p uploads outputs
chmod 755 uploads outputs

# Create systemd service
echo "Creating systemd service..."
sudo tee /etc/systemd/system/singlecell-converter.service > /dev/null <<EOF
[Unit]
Description=Single-cell File Converter Web Application
After=network.target

[Service]
Type=simple
User=$USER
WorkingDirectory=$INSTALL_DIR
Environment="PATH=$INSTALL_DIR/venv/bin"
ExecStart=$INSTALL_DIR/venv/bin/gunicorn \\
    --bind 127.0.0.1:5000 \\
    --workers 2 \\
    --timeout 300 \\
    --access-logfile /var/log/singlecell-converter/access.log \\
    --error-logfile /var/log/singlecell-converter/error.log \\
    app:app
Restart=always

[Install]
WantedBy=multi-user.target
EOF

# Create log directory
sudo mkdir -p /var/log/singlecell-converter
sudo chown $USER:$USER /var/log/singlecell-converter

# Enable and start service
echo "Enabling and starting service..."
sudo systemctl daemon-reload
sudo systemctl enable singlecell-converter
sudo systemctl start singlecell-converter

# Configure Nginx
echo "Configuring Nginx..."
read -p "Enter your domain name (or IP address): " DOMAIN

sudo tee /etc/nginx/sites-available/singlecell-converter > /dev/null <<EOF
server {
    listen 80;
    server_name $DOMAIN;

    client_max_body_size 500M;

    location / {
        proxy_pass http://127.0.0.1:5000;
        proxy_set_header Host \$host;
        proxy_set_header X-Real-IP \$remote_addr;
        proxy_set_header X-Forwarded-For \$proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto \$scheme;

        # タイムアウト設定
        proxy_read_timeout 600s;
        proxy_connect_timeout 600s;
        proxy_send_timeout 600s;
    }
}
EOF

# Enable Nginx site
sudo ln -sf /etc/nginx/sites-available/singlecell-converter /etc/nginx/sites-enabled/
sudo rm -f /etc/nginx/sites-enabled/default

# Test Nginx configuration
sudo nginx -t

# Restart Nginx
sudo systemctl restart nginx

# Configure firewall
echo "Configuring firewall..."
sudo apt install -y ufw
sudo ufw --force enable
sudo ufw allow 22/tcp
sudo ufw allow 80/tcp
sudo ufw allow 443/tcp

echo ""
echo "=========================================="
echo "Deployment completed!"
echo "=========================================="
echo ""
echo "Application is now running at: http://$DOMAIN"
echo ""
echo "Useful commands:"
echo "  - Check status: sudo systemctl status singlecell-converter"
echo "  - View logs: sudo journalctl -u singlecell-converter -f"
echo "  - Restart: sudo systemctl restart singlecell-converter"
echo ""
echo "To setup SSL certificate (recommended):"
echo "  sudo apt install certbot python3-certbot-nginx"
echo "  sudo certbot --nginx -d $DOMAIN"
echo ""
echo "Configuration file: $INSTALL_DIR/.env"
echo ""
