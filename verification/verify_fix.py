
from playwright.sync_api import sync_playwright
import os

def run_verification():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        cwd = os.getcwd()
        page.goto(f'file://{cwd}/index.html')

        page.click('#mode-toggle')

        # 1. Linear Matrix
        page.click('.app-card:has-text("Linear")')
        page.wait_for_selector('#tool-linear', state='visible')

        # Now click the back button specifically inside #tool-linear
        page.click('#tool-linear .back-btn')
        page.wait_for_selector('#app-grid', state='visible')

        # 2. Converter
        page.click('.app-card:has-text("Converter")')
        page.wait_for_selector('#tool-converter', state='visible')

        page.select_option('#conv-category', 'length')
        page.wait_for_timeout(200)
        page.select_option('#conv-unit-from', 'm')
        page.select_option('#conv-unit-to', 'ft')
        page.fill('#conv-val-from', '1')

        # Screenshot
        page.screenshot(path='verification/converter_fixed.png')

        val_to = page.input_value('#conv-val-to')
        print(f'1m = {val_to} ft')

        browser.close()

if __name__ == '__main__':
    run_verification()
