// Single-cell File Converter JavaScript

let selectedFiles = [];
let batchMode = false;

document.addEventListener('DOMContentLoaded', function () {
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
    const batchToggle = document.getElementById('batchMode');
    const fileList = document.getElementById('fileList');
    const fileListActions = document.getElementById('fileListActions');
    const addMoreBtn = document.getElementById('addMoreBtn');
    const clearAllBtn = document.getElementById('clearAllBtn');
    const uploadText = document.getElementById('uploadText');
    const convertNote = document.getElementById('convertNote');

    function getSelectedFileType() {
        const checked = document.querySelector('input[name="file_type"]:checked');
        return checked ? checked.value : 'auto';
    }

    function updateCsvOptionsVisibility() {
        const ftype = getSelectedFileType();
        const showCsvOptions = ftype === 'csv' || ftype === 'auto';
        csvOptions.style.display = showCsvOptions ? 'block' : 'none';
    }

    function updateUploadUiMode() {
        if (batchMode) {
            fileInput.setAttribute('multiple', 'multiple');
            uploadText.textContent = '複数のファイルをドラッグ&ドロップ';
            convertNote.innerHTML = '出力: 各ファイルが <code>convert_元のファイル名.h5mu</code> となり、まとめて ZIP でダウンロードされます';
        } else {
            fileInput.removeAttribute('multiple');
            uploadText.textContent = 'ファイルをドラッグ&ドロップ';
            convertNote.innerHTML = '出力: <code>convert_元のファイル名.h5mu</code>';
        }
    }

    fileTypeRadios.forEach((radio) => {
        radio.addEventListener('change', function () {
            updateCsvOptionsVisibility();
        });
    });

    batchToggle.addEventListener('change', function () {
        batchMode = this.checked;
        resetFileSelection();
        updateUploadUiMode();
    });

    uploadArea.addEventListener('click', function () {
        fileInput.click();
    });

    fileInput.addEventListener('change', function (e) {
        handleFilesSelected(Array.from(e.target.files));
        // Reset native input so selecting the same file again still triggers change
        fileInput.value = '';
    });

    uploadArea.addEventListener('dragover', function (e) {
        e.preventDefault();
        e.stopPropagation();
        uploadArea.classList.add('drag-over');
    });

    uploadArea.addEventListener('dragleave', function (e) {
        e.preventDefault();
        e.stopPropagation();
        uploadArea.classList.remove('drag-over');
    });

    uploadArea.addEventListener('drop', function (e) {
        e.preventDefault();
        e.stopPropagation();
        uploadArea.classList.remove('drag-over');

        const files = Array.from(e.dataTransfer.files);
        if (files.length > 0) {
            handleFilesSelected(files);
        }
    });

    convertBtn.addEventListener('click', function () {
        if (selectedFiles.length === 0) return;
        if (batchMode || selectedFiles.length > 1) {
            convertBatch();
        } else {
            convertSingle();
        }
    });

    if (removeFileBtn) {
        removeFileBtn.addEventListener('click', function (e) {
            e.stopPropagation();
            resetFileSelection();
        });
    }

    const btnBrowse = document.querySelector('.upload-area .btn-browse');
    if (btnBrowse) {
        btnBrowse.addEventListener('click', function (e) {
            e.stopPropagation();
            fileInput.click();
        });
    }

    if (addMoreBtn) {
        addMoreBtn.addEventListener('click', function (e) {
            e.stopPropagation();
            fileInput.click();
        });
    }

    if (clearAllBtn) {
        clearAllBtn.addEventListener('click', function (e) {
            e.stopPropagation();
            resetFileSelection();
        });
    }

    function handleFilesSelected(files) {
        if (!files || files.length === 0) return;

        if (batchMode) {
            // Append, dedupe by name+size
            const key = (f) => f.name + ':' + f.size;
            const existing = new Set(selectedFiles.map(key));
            for (const f of files) {
                if (!existing.has(key(f))) {
                    selectedFiles.push(f);
                    existing.add(key(f));
                }
            }
        } else {
            // Single mode: only keep the first file
            selectedFiles = [files[0]];
        }

        renderSelection();
        hideStatus();
    }

    function renderSelection() {
        if (selectedFiles.length === 0) {
            uploadArea.style.display = 'block';
            fileInfo.style.display = 'none';
            fileList.style.display = 'none';
            fileListActions.style.display = 'none';
            convertBtn.disabled = true;
            return;
        }

        uploadArea.style.display = 'none';
        convertBtn.disabled = false;

        if (selectedFiles.length === 1 && !batchMode) {
            // Single-file legacy view
            const f = selectedFiles[0];
            fileName.textContent = f.name;
            fileSize.textContent = formatFileSize(f.size);
            fileInfo.style.display = 'flex';
            fileList.style.display = 'none';
            fileListActions.style.display = 'none';
        } else {
            // Multi-file list view
            fileInfo.style.display = 'none';
            fileList.innerHTML = '';
            selectedFiles.forEach((f, idx) => {
                const li = document.createElement('li');

                const icon = document.createElement('span');
                icon.className = 'item-icon';
                icon.textContent = iconForFile(f.name);
                li.appendChild(icon);

                const name = document.createElement('span');
                name.className = 'item-name';
                name.textContent = f.name;
                li.appendChild(name);

                const typeBadge = document.createElement('span');
                typeBadge.className = 'item-type-badge';
                typeBadge.textContent = inferTypeLabel(f.name);
                li.appendChild(typeBadge);

                const size = document.createElement('span');
                size.className = 'item-size';
                size.textContent = formatFileSize(f.size);
                li.appendChild(size);

                const rm = document.createElement('button');
                rm.type = 'button';
                rm.className = 'item-remove';
                rm.textContent = '✕';
                rm.addEventListener('click', function (e) {
                    e.stopPropagation();
                    selectedFiles.splice(idx, 1);
                    renderSelection();
                });
                li.appendChild(rm);

                fileList.appendChild(li);
            });
            fileList.style.display = 'flex';
            fileListActions.style.display = 'flex';
        }
    }

    function iconForFile(filename) {
        const lower = filename.toLowerCase();
        if (lower.endsWith('.csv')) return '📊';
        if (lower.endsWith('.h5ad')) return '🧪';
        if (lower.endsWith('.rds')) return '📦';
        if (lower.endsWith('.zip') || lower.endsWith('.tar.gz') || lower.endsWith('.tar')) return '🧬';
        return '📄';
    }

    function inferTypeLabel(filename) {
        const lower = filename.toLowerCase();
        if (lower.endsWith('.csv')) return 'CSV';
        if (lower.endsWith('.h5ad')) return 'h5ad';
        if (lower.endsWith('.rds')) return 'RDS';
        if (lower.endsWith('.zip') || lower.endsWith('.tar.gz') || lower.endsWith('.tar')) return 'MEX';
        return '?';
    }

    function resetFileSelection() {
        selectedFiles = [];
        fileInput.value = '';
        renderSelection();
        hideStatus();
    }

    function buildCommonForm() {
        const fd = new FormData();
        fd.append('file_type', getSelectedFileType());
        // transpose: hidden input now holds 'auto' (default). When the backend
        // sees 'auto' it runs the gene-name heuristic.
        const transposeEl = document.getElementById('transpose');
        const transposeVal = transposeEl ? (transposeEl.value || 'auto') : 'auto';
        fd.append('transpose', transposeVal);
        fd.append('has_header', document.getElementById('hasHeader').checked);
        fd.append('has_index', document.getElementById('hasIndex').checked);
        const umapEl = document.getElementById('makeUmap');
        if (umapEl) fd.append('make_umap', umapEl.checked);
        return fd;
    }

    function convertSingle() {
        const formData = buildCommonForm();
        formData.append('file', selectedFiles[0]);
        // Backwards-compatible single-file endpoint expects 'file_type' as a
        // concrete type, but our backend now accepts 'auto' too.

        const wantUmap = document.getElementById('makeUmap') && document.getElementById('makeUmap').checked;
        showStatus(wantUmap ? '変換 + UMAP を実行中...（数十秒〜数分かかります）' : '変換中...', 'info');
        progressBar.style.display = 'block';
        convertBtn.disabled = true;

        fetch('/convert', { method: 'POST', body: formData })
            .then(parseJsonOrThrow)
            .then((data) => {
                progressBar.style.display = 'none';
                if (data.success) {
                    showStatus('✅ ' + data.message, 'success');
                    if (data.detection) renderDetection(data.detection);
                    if (data.umap) renderInlineUmap(data.umap);
                    const hasUmapJob = !!(data.umap_job && data.umap_job.job_id);
                    if (hasUmapJob) {
                        startUmapJobPolling(data.umap_job.job_id);
                    }
                    triggerDownload(data.download_url, { preserveStatus: hasUmapJob });
                } else {
                    showStatus('❌ ' + (data.error || '変換に失敗しました'), 'error');
                    convertBtn.disabled = false;
                }
            })
            .catch((error) => {
                progressBar.style.display = 'none';
                showStatus('❌ エラー: ' + error.message, 'error');
                convertBtn.disabled = false;
            });
    }

    function renderDetection(det) {
        if (!det || !det.auto_detected) return;
        const old = document.getElementById('detectionInfo');
        if (old) old.remove();

        const wrap = document.createElement('div');
        wrap.id = 'detectionInfo';
        wrap.className = 'batch-summary';
        const head = document.createElement('strong');
        head.textContent = det.transpose_applied
            ? '🔄 自動転置を適用しました'
            : '✓ 転置不要と判定しました（そのまま使用）';
        wrap.appendChild(head);

        if (det.reason) {
            const p = document.createElement('div');
            p.className = 'muted';
            p.style.marginTop = '0.25rem';
            p.textContent = '理由: ' + det.reason;
            wrap.appendChild(p);
        }
        if (det.scores) {
            const small = document.createElement('div');
            small.className = 'muted';
            small.style.fontSize = '0.75rem';
            small.style.marginTop = '0.25rem';
            const s = det.scores;
            small.textContent = `スコア — 行(gene): ${s.row_gene_score ?? '?'}, 列(gene): ${s.col_gene_score ?? '?'}, 行(barcode): ${s.row_barcode_score ?? '?'}, 列(barcode): ${s.col_barcode_score ?? '?'}`;
            wrap.appendChild(small);
        }

        statusSection.appendChild(wrap);
    }

    function startUmapJobPolling(jobId) {
        // Build / show a live status panel for the running job.
        const old = document.getElementById('umapJobPanel');
        if (old) old.remove();

        const panel = document.createElement('div');
        panel.id = 'umapJobPanel';
        panel.className = 'convert-umap-preview';
        panel.style.alignItems = 'flex-start';

        const head = document.createElement('strong');
        head.textContent = '🌐 UMAP 生成中...';
        panel.appendChild(head);

        const stage = document.createElement('div');
        stage.className = 'muted';
        stage.id = 'umapJobStage';
        stage.textContent = '開始中...';
        panel.appendChild(stage);

        const elapsed = document.createElement('div');
        elapsed.className = 'muted';
        elapsed.id = 'umapJobElapsed';
        elapsed.style.fontSize = '0.75rem';
        elapsed.textContent = '経過時間: 0s';
        panel.appendChild(elapsed);

        const bar = document.createElement('div');
        bar.className = 'progress-bar';
        bar.style.width = '100%';
        bar.style.marginTop = '0.5rem';
        const fill = document.createElement('div');
        fill.className = 'progress-fill';
        bar.appendChild(fill);
        panel.appendChild(bar);

        statusSection.appendChild(panel);

        const stageEl = panel.querySelector('#umapJobStage');
        const elapsedEl = panel.querySelector('#umapJobElapsed');

        let stopped = false;
        const poll = () => {
            if (stopped) return;
            fetch(`/umap_progress/${jobId}`)
                .then((r) => r.json())
                .then((data) => {
                    if (!data.success) {
                        stopped = true;
                        stageEl.textContent = '⚠️ ' + (data.error || 'ジョブを取得できません');
                        return;
                    }
                    if (data.stage) stageEl.textContent = '現在: ' + data.stage;
                    if (typeof data.elapsed_sec === 'number') {
                        elapsedEl.textContent = '経過時間: ' + Math.round(data.elapsed_sec) + 's';
                    }

                    if (data.status === 'done') {
                        stopped = true;
                        bar.remove();
                        head.textContent = '✅ UMAP が完成しました';
                        const r = data.result || {};
                        renderInlineUmap({
                            success: true,
                            preview_url: r.preview_url,
                            download_url: r.download_url,
                            modality: r.modality,
                            n_cells: r.n_cells,
                            n_genes: r.n_genes,
                        });
                        return;
                    }
                    if (data.status === 'error') {
                        stopped = true;
                        bar.remove();
                        head.textContent = '❌ UMAP に失敗しました';
                        stageEl.textContent = data.error || '不明なエラー';
                        return;
                    }
                    setTimeout(poll, 1500);
                })
                .catch((err) => {
                    setTimeout(poll, 3000);
                });
        };
        poll();
    }

    function renderInlineUmap(umap) {
        const oldPreview = document.getElementById('inlineUmapPreview');
        if (oldPreview) oldPreview.remove();

        const wrap = document.createElement('div');
        wrap.id = 'inlineUmapPreview';
        wrap.className = 'convert-umap-preview';

        if (umap.success && umap.preview_url) {
            const head = document.createElement('strong');
            head.textContent = `🌐 UMAP（${umap.modality || 'rna'} / ${umap.n_cells || '?'} cells × ${umap.n_genes || '?'} genes）`;
            wrap.appendChild(head);

            const img = document.createElement('img');
            img.src = umap.preview_url;
            img.alt = 'UMAP';
            wrap.appendChild(img);

            const dl = document.createElement('a');
            dl.href = umap.download_url;
            dl.className = 'btn-browse';
            dl.textContent = 'PNG をダウンロード';
            dl.setAttribute('download', '');
            wrap.appendChild(dl);
        } else {
            const err = document.createElement('span');
            err.style.color = 'var(--error-color)';
            err.textContent = '⚠️ ' + (umap.error || 'UMAP の生成に失敗しました');
            wrap.appendChild(err);
        }

        statusSection.appendChild(wrap);
    }

    function convertBatch() {
        const formData = buildCommonForm();
        selectedFiles.forEach((f) => formData.append('files', f));

        showStatus(`変換中... (${selectedFiles.length} 件)`, 'info');
        progressBar.style.display = 'block';
        convertBtn.disabled = true;

        fetch('/convert_batch', { method: 'POST', body: formData })
            .then(parseJsonOrThrow)
            .then((data) => {
                progressBar.style.display = 'none';
                if (data.success) {
                    showStatus('✅ ' + data.message, 'success');
                    renderBatchSummary(data.results);
                    triggerDownload(data.download_url);
                } else {
                    showStatus('❌ ' + (data.error || '変換に失敗しました'), 'error');
                    if (data.results) renderBatchSummary(data.results);
                    convertBtn.disabled = false;
                }
            })
            .catch((error) => {
                progressBar.style.display = 'none';
                showStatus('❌ エラー: ' + error.message, 'error');
                convertBtn.disabled = false;
            });
    }

    function renderBatchSummary(results) {
        if (!Array.isArray(results) || results.length === 0) return;

        const oldSummary = document.getElementById('batchSummary');
        if (oldSummary) oldSummary.remove();

        const wrapper = document.createElement('div');
        wrapper.id = 'batchSummary';
        wrapper.className = 'batch-summary';

        const succeeded = results.filter((r) => r.success).length;
        const failed = results.length - succeeded;
        const head = document.createElement('strong');
        head.textContent = `結果: 成功 ${succeeded} 件 / 失敗 ${failed} 件`;
        wrapper.appendChild(head);

        const ul = document.createElement('ul');
        results.forEach((r) => {
            const li = document.createElement('li');
            li.className = r.success ? 'succeeded' : 'failed';
            li.textContent = r.success
                ? `✓ ${r.filename}${r.file_type ? ' [' + r.file_type + ']' : ''}`
                : `✗ ${r.filename}${r.file_type ? ' [' + r.file_type + ']' : ''} — ${r.error || '不明なエラー'}`;
            ul.appendChild(li);
        });
        wrapper.appendChild(ul);

        statusSection.appendChild(wrapper);
    }

    function triggerDownload(url, options) {
        options = options || {};
        setTimeout(() => {
            window.location.href = url;
            // When a follow-up job (e.g. UMAP) is running we must keep the
            // status section visible so the user can see its progress and
            // the resulting image. Skip the auto reset in that case.
            if (options.preserveStatus) return;
            setTimeout(() => {
                resetFileSelection();
            }, 1000);
        }, 500);
    }

    function parseJsonOrThrow(response) {
        return response.json().then((data) => {
            if (!response.ok && !data.success) {
                const err = new Error(data.error || '変換に失敗しました');
                err.payload = data;
                throw err;
            }
            return data;
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
        const oldSummary = document.getElementById('batchSummary');
        if (oldSummary) oldSummary.remove();
        if (statusSection) {
            statusSection.style.display = 'none';
        }
    }

    function formatFileSize(bytes) {
        if (bytes === 0) return '0 Bytes';

        const k = 1024;
        const sizes = ['Bytes', 'KB', 'MB', 'GB'];
        const i = Math.floor(Math.log(bytes) / Math.log(k));

        return Math.round((bytes / Math.pow(k, i)) * 100) / 100 + ' ' + sizes[i];
    }

    // Initial state
    updateCsvOptionsVisibility();
    updateUploadUiMode();
});
