"""Optional Chromium layout tests; uses synthetic published evidence only."""
from pathlib import Path
import os
import shutil
import sys
import unittest

from playwright.sync_api import sync_playwright

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tests'))
from test_evidence_report import sample_data
from qpcr_pipeline.researcher_report_evidence import render_researcher_report_html


class EvidenceLayoutTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.playwright = sync_playwright().start()
        executable = os.environ.get('GEISON_CHROMIUM_EXECUTABLE') or shutil.which('chromium')
        options = {'executable_path': executable} if executable else {}
        cls.browser = cls.playwright.chromium.launch(headless=True, **options)

    @classmethod
    def tearDownClass(cls):
        cls.browser.close()
        cls.playwright.stop()

    def test_desktop_and_mobile_have_no_page_overflow_or_remote_requests(self):
        html = render_researcher_report_html(sample_data())
        for width in (390, 1440):
            with self.subTest(width=width):
                page = self.browser.new_page(viewport={'width': width, 'height': 900})
                requests, errors = [], []
                page.on('request', lambda r: requests.append(r.url))
                page.on('pageerror', lambda e: errors.append(str(e)))
                try:
                    page.set_content(html, wait_until='load')
                    self.assertEqual(page.locator('h1').inner_text(), 'Geison Evidence Report')
                    self.assertFalse(page.evaluate('document.documentElement.scrollWidth > innerWidth'))
                    self.assertEqual(requests, [])
                    self.assertEqual(errors, [])
                finally:
                    page.close()
