#!/usr/bin/env python3
"""
将HTML文件转换为PDF
使用Playwright自动化浏览器打印功能
"""
import subprocess
import sys
import os

def convert_html_to_pdf_via_chrome(html_path, pdf_path):
    """使用Chrome的headless模式将HTML转换为PDF"""
    chrome_paths = [
        '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
        '/Applications/Chromium.app/Contents/MacOS/Chromium',
        '/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge',
    ]

    chrome_bin = None
    for path in chrome_paths:
        if os.path.exists(path):
            chrome_bin = path
            break

    if not chrome_bin:
        print("错误：未找到Chrome、Chromium或Edge浏览器")
        print("请手动在浏览器中打开 socp_acopf_thesis.html 并使用打印功能保存为PDF")
        return False

    abs_html_path = os.path.abspath(html_path)
    abs_pdf_path = os.path.abspath(pdf_path)

    cmd = [
        chrome_bin,
        '--headless',
        '--disable-gpu',
        '--print-to-pdf=' + abs_pdf_path,
        '--no-pdf-header-footer',
        'file://' + abs_html_path
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        if os.path.exists(abs_pdf_path):
            print(f"✓ PDF已成功创建: {pdf_path}")
            return True
        else:
            print(f"转换失败: {result.stderr}")
            return False
    except Exception as e:
        print(f"错误: {e}")
        return False

if __name__ == '__main__':
    html_file = 'socp_acopf_thesis.html'
    pdf_file = 'socp_acopf_thesis.pdf'

    if not os.path.exists(html_file):
        print(f"错误：HTML文件不存在: {html_file}")
        sys.exit(1)

    success = convert_html_to_pdf_via_chrome(html_file, pdf_file)
    sys.exit(0 if success else 1)
