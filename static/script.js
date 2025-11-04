// Single-cell File Converter JavaScript

let selectedFile = null;

document.addEventListener('DOMContentLoaded', function() {
    const fileInput = document.getElementById('fileInput');
    const uploadArea = document.getElementById('uploadArea');
    const convertBtn = document.getElementById('convertBtn');
    const fileInfo = document.getElementById('fileInfo');
    const fileName = document.getElementById('fileName');
    const fileSize = document.getElementById('fileSize');
    const statusMessage = document.getElementById('statusMessage');
    const progressBar = document.getElementById('progressBar');
    const csvOptions = document.getElementById('csvOptionsSection');
    const removeFileBtn = document.getElementById('removeFileBtn');
    const statusSection = document.getElementById('statusSection');
    const fileTypeRadios = document.querySelectorAll('input[name="file_type"]');

    // File type change handler
    fileTypeRadios.forEach(radio => {
        radio.addEventListener('change', function() {
            // Show/hide CSV options based on file type
            if (this.value === 'csv') {
                csvOptions.style.display = 'block';
            } else {
                csvOptions.style.display = 'none';
            }

            // Reset file selection
            resetFileSelection();
        });
    });

    // Upload area click handler
    uploadArea.addEventListener('click', function() {
        fileInput.click();
    });

    // File input change handler
    fileInput.addEventListener('change', function(e) {
        handleFileSelect(e.target.files[0]);
    });

    // Drag and drop handlers
    uploadArea.addEventListener('dragover', function(e) {
        e.preventDefault();
        e.stopPropagation();
        uploadArea.classList.add('drag-over');
    });

    uploadArea.addEventListener('dragleave', function(e) {
        e.preventDefault();
        e.stopPropagation();
        uploadArea.classList.remove('drag-over');
    });

    uploadArea.addEventListener('drop', function(e) {
        e.preventDefault();
        e.stopPropagation();
        uploadArea.classList.remove('drag-over');

        const files = e.dataTransfer.files;
        if (files.length > 0) {
            handleFileSelect(files[0]);
        }
    });

    // Convert button click handler
    convertBtn.addEventListener('click', function() {
        if (selectedFile) {
            convertFile();
        }
    });

    // Remove file button click handler
    if (removeFileBtn) {
        removeFileBtn.addEventListener('click', function(e) {
            e.stopPropagation();
            resetFileSelection();
        });
    }

    // Browse button click handler
    const btnBrowse = document.querySelector('.btn-browse');
    if (btnBrowse) {
        btnBrowse.addEventListener('click', function(e) {
            e.stopPropagation();
            fileInput.click();
        });
    }

    function handleFileSelect(file) {
        if (!file) return;

        selectedFile = file;

        // Display file information
        fileName.textContent = file.name;
        fileSize.textContent = formatFileSize(file.size);
        fileInfo.style.display = 'block';
        convertBtn.disabled = false;

        // Hide upload area
        uploadArea.style.display = 'none';

        hideStatus();
    }

    function resetFileSelection() {
        selectedFile = null;
        fileInput.value = '';
        fileInfo.style.display = 'none';
        uploadArea.style.display = 'block';
        convertBtn.disabled = true;
        hideStatus();
    }

    function convertFile() {
        const formData = new FormData();
        formData.append('file', selectedFile);

        // Get selected file type
        const fileType = document.querySelector('input[name="file_type"]:checked').value;
        formData.append('file_type', fileType);

        // Get CSV options if CSV is selected
        if (fileType === 'csv') {
            formData.append('transpose', document.getElementById('transpose').checked);
            formData.append('has_header', document.getElementById('hasHeader').checked);
            formData.append('has_index', document.getElementById('hasIndex').checked);
        }

        // Show progress
        showStatus('変換中...', 'info');
        progressBar.style.display = 'block';
        convertBtn.disabled = true;

        // Send request
        fetch('/convert', {
            method: 'POST',
            body: formData
        })
        .then(response => {
            if (!response.ok) {
                return response.json().then(data => {
                    throw new Error(data.error || '変換に失敗しました');
                });
            }
            return response.json();
        })
        .then(data => {
            progressBar.style.display = 'none';

            if (data.success) {
                showStatus('✅ ' + data.message, 'success');

                // Automatically download the file
                setTimeout(() => {
                    window.location.href = data.download_url;

                    // Reset after download
                    setTimeout(() => {
                        resetFileSelection();
                    }, 1000);
                }, 500);
            } else {
                showStatus('❌ ' + (data.error || '変換に失敗しました'), 'error');
                convertBtn.disabled = false;
            }
        })
        .catch(error => {
            progressBar.style.display = 'none';
            showStatus('❌ エラー: ' + error.message, 'error');
            convertBtn.disabled = false;
        });
    }

    function showStatus(message, type) {
        statusMessage.textContent = message;
        statusMessage.className = 'status-message ' + type;
        statusMessage.style.display = 'block';
        if (statusSection) {
            statusSection.style.display = 'block';
        }
    }

    function hideStatus() {
        statusMessage.style.display = 'none';
        progressBar.style.display = 'none';
        if (statusSection) {
            statusSection.style.display = 'none';
        }
    }

    function formatFileSize(bytes) {
        if (bytes === 0) return '0 Bytes';

        const k = 1024;
        const sizes = ['Bytes', 'KB', 'MB', 'GB'];
        const i = Math.floor(Math.log(bytes) / Math.log(k));

        return Math.round(bytes / Math.pow(k, i) * 100) / 100 + ' ' + sizes[i];
    }
});
