"""
Single-cell file format converter web application
Converts CSV, RDS, and MEX files to h5mu format
"""
from flask import Flask, render_template, request, send_file, jsonify
import os
from pathlib import Path
from werkzeug.utils import secure_filename
import uuid
from datetime import datetime
from dotenv import load_dotenv

from converters.csv_converter import csv_to_h5mu, validate_csv
from converters.rds_converter import rds_to_h5mu, check_rpy2_available
from converters.mex_converter import mex_to_h5mu, validate_mex_directory
import shutil
import zipfile

# Load environment variables
load_dotenv()

app = Flask(__name__)

# Configuration from environment variables
app.config['SECRET_KEY'] = os.getenv('SECRET_KEY', 'dev-secret-key-change-in-production')
app.config['UPLOAD_FOLDER'] = os.getenv('UPLOAD_FOLDER', 'uploads')
app.config['OUTPUT_FOLDER'] = os.getenv('OUTPUT_FOLDER', 'outputs')
app.config['MAX_CONTENT_LENGTH'] = int(os.getenv('MAX_CONTENT_LENGTH', 500 * 1024 * 1024))

# Allowed file extensions
ALLOWED_EXTENSIONS = {'csv', 'rds', 'mtx', 'tsv', 'gz', 'zip'}


def allowed_file(filename):
    """Check if file extension is allowed"""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def generate_unique_filename(prefix='file'):
    """Generate unique filename with timestamp"""
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    unique_id = str(uuid.uuid4())[:8]
    return f"{prefix}_{timestamp}_{unique_id}"


def find_mex_directory(root_dir):
    """
    Find directory containing MEX files (matrix.mtx, features/genes.tsv, barcodes.tsv)

    Parameters:
    -----------
    root_dir : str
        Root directory to search

    Returns:
    --------
    str or None : Path to MEX directory or None if not found
    """
    from converters.mex_converter import find_file

    root_path = Path(root_dir)

    # Check root directory first
    if find_file(root_path, ['matrix.mtx', 'matrix.mtx.gz']):
        return str(root_path)

    # Search subdirectories
    for dirpath, dirnames, filenames in os.walk(root_dir):
        dir_path = Path(dirpath)
        if find_file(dir_path, ['matrix.mtx', 'matrix.mtx.gz']):
            return str(dir_path)

    return None


@app.route('/')
def index():
    """Main page"""
    return render_template('index.html')


@app.route('/convert', methods=['POST'])
def convert():
    """Handle file conversion"""
    try:
        # Check if file was uploaded
        if 'file' not in request.files:
            return jsonify({'error': 'ファイルが選択されていません'}), 400

        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'ファイルが選択されていません'}), 400

        if not allowed_file(file.filename):
            return jsonify({'error': '対応していないファイル形式です'}), 400

        # Get conversion options
        file_type = request.form.get('file_type', 'csv')
        transpose = request.form.get('transpose', 'false') == 'true'
        has_header = request.form.get('has_header', 'true') == 'true'
        has_index = request.form.get('has_index', 'true') == 'true'

        # Save uploaded file
        filename = secure_filename(file.filename)
        unique_name = generate_unique_filename('upload')
        file_extension = filename.rsplit('.', 1)[1].lower()
        upload_path = os.path.join(app.config['UPLOAD_FOLDER'], f"{unique_name}.{file_extension}")
        file.save(upload_path)

        # Generate output filename
        output_filename = f"{generate_unique_filename('converted')}.h5mu"
        output_path = os.path.join(app.config['OUTPUT_FOLDER'], output_filename)

        # Convert based on file type
        if file_type == 'csv':
            csv_to_h5mu(
                upload_path,
                output_path,
                transpose=transpose,
                has_header=has_header,
                has_index=has_index
            )
            # Clean up uploaded file
            os.remove(upload_path)

        elif file_type == 'rds':
            # Check if rpy2 is available
            if not check_rpy2_available():
                return jsonify({
                    'error': 'RDS変換にはrpy2が必要です。インストール方法については README.md を参照してください。'
                }), 501

            rds_to_h5mu(upload_path, output_path, object_type='auto')
            # Clean up uploaded file
            os.remove(upload_path)

        elif file_type == 'mex':
            # MEX files can be uploaded as a zip file or individual files
            # If it's a zip, extract it first
            if file_extension == 'gz':
                # This might be a single .mtx.gz or .tsv.gz file
                # For MEX, we need a directory with multiple files
                return jsonify({
                    'error': 'MEXファイルは、matrix.mtx、features.tsv、barcodes.tsvを含むZIPファイルとしてアップロードしてください。'
                }), 400
            elif file_extension == 'zip':
                # Extract zip file
                extract_dir = os.path.join(app.config['UPLOAD_FOLDER'], unique_name + '_extracted')
                os.makedirs(extract_dir, exist_ok=True)
                with zipfile.ZipFile(upload_path, 'r') as zip_ref:
                    zip_ref.extractall(extract_dir)

                # Find MEX directory (might be in a subdirectory)
                mex_dir = find_mex_directory(extract_dir)
                if not mex_dir:
                    shutil.rmtree(extract_dir)
                    os.remove(upload_path)
                    return jsonify({
                        'error': 'ZIPファイル内にMEXファイル（matrix.mtx、features.tsv、barcodes.tsv）が見つかりません。'
                    }), 400

                mex_to_h5mu(mex_dir, output_path)

                # Clean up
                shutil.rmtree(extract_dir)
                os.remove(upload_path)
            else:
                return jsonify({
                    'error': 'MEXファイルはZIP形式でアップロードしてください。'
                }), 400

        else:
            os.remove(upload_path)
            return jsonify({'error': '不明なファイル形式です'}), 400

        # Return success with download link
        return jsonify({
            'success': True,
            'message': '変換が完了しました',
            'download_url': f'/download/{output_filename}'
        })

    except Exception as e:
        return jsonify({'error': f'変換エラー: {str(e)}'}), 500


@app.route('/download/<filename>')
def download(filename):
    """Download converted file"""
    try:
        filepath = os.path.join(app.config['OUTPUT_FOLDER'], secure_filename(filename))
        if not os.path.exists(filepath):
            return jsonify({'error': 'ファイルが見つかりません'}), 404

        return send_file(
            filepath,
            as_attachment=True,
            download_name=filename,
            mimetype='application/octet-stream'
        )
    except Exception as e:
        return jsonify({'error': f'ダウンロードエラー: {str(e)}'}), 500


@app.route('/validate', methods=['POST'])
def validate():
    """Validate uploaded file"""
    try:
        if 'file' not in request.files:
            return jsonify({'error': 'ファイルが選択されていません'}), 400

        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'ファイルが選択されていません'}), 400

        # Save temporarily
        filename = secure_filename(file.filename)
        unique_name = generate_unique_filename('validate')
        file_extension = filename.rsplit('.', 1)[1].lower()
        temp_path = os.path.join(app.config['UPLOAD_FOLDER'], f"{unique_name}.{file_extension}")
        file.save(temp_path)

        # Validate based on file type
        file_type = request.form.get('file_type', 'csv')

        if file_type == 'csv':
            info = validate_csv(temp_path)
        else:
            info = {'message': 'Validation not yet implemented for this file type'}

        # Clean up
        os.remove(temp_path)

        return jsonify({'success': True, 'info': info})

    except Exception as e:
        return jsonify({'error': f'検証エラー: {str(e)}'}), 500


# Create necessary directories on startup
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['OUTPUT_FOLDER'], exist_ok=True)


if __name__ == '__main__':
    # Get configuration from environment
    port = int(os.getenv('PORT', 5000))
    host = os.getenv('HOST', '0.0.0.0')
    debug = os.getenv('FLASK_ENV', 'development') == 'development'

    # Run the app
    print("Starting Single-cell File Converter Web Application...")
    print(f"Access the application at: http://{host}:{port}")
    print(f"Debug mode: {debug}")
    app.run(debug=debug, host=host, port=port)
