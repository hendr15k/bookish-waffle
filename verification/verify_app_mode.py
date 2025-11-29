
from playwright.sync_api import sync_playwright
import os

def run_verification():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        # Load the index.html directly
        cwd = os.getcwd()
        page.goto(f'file://{cwd}/index.html')

        # Enable App Mode
        page.click('#mode-toggle')
        page.wait_for_selector('.tool-view', state='hidden')

        # 1. Verify Calculator Clear Button
        page.click('.app-card:has-text("Calculator")')
        page.wait_for_selector('#tool-calculator', state='visible')
        page.fill('#calc-main-input', '1+1')
        page.click('#tool-calculator button:text("×")') # Click the X clear button
        val = page.input_value('#calc-main-input')
        if val != '':
             print(f'Calculator clear failed, value is {val}')

        page.click('button:text("←")') # Back

        # 2. Verify Matrix Builder
        page.click('.app-card:has-text("Linear")')
        page.wait_for_selector('#tool-linear', state='visible')
        page.fill('#mat-rows', '2')
        page.fill('#mat-cols', '2')
        page.click('button:text("Generate Grid")')
        page.wait_for_selector('.matrix-cell')

        # Fill matrix grid
        page.fill('.matrix-cell[data-row="0"][data-col="0"]', '1')
        page.fill('.matrix-cell[data-row="0"][data-col="1"]', '2')
        page.fill('.matrix-cell[data-row="1"][data-col="0"]', '3')
        page.fill('.matrix-cell[data-row="1"][data-col="1"]', '4')

        page.click('button:text("Insert Matrix")')
        lin_expr = page.input_value('#lin-expr')
        if lin_expr != '[[1,2],[3,4]]':
            print(f'Matrix insert failed: {lin_expr}')

        page.screenshot(path='verification/linear_matrix.png')

        # 3. Verify Unit Converter
        page.click('button:text("←")') # Back
        page.click('.app-card:has-text("Converter")')
        page.wait_for_selector('#tool-converter', state='visible')

        # Select Length, m to ft
        page.select_option('#conv-category', 'length')
        page.wait_for_timeout(100) # Wait for units to populate
        page.select_option('#conv-unit-from', 'm')
        page.select_option('#conv-unit-to', 'ft')

        page.fill('#conv-val-from', '1')
        # Expect ~3.28
        val_to = page.input_value('#conv-val-to')
        print(f'1m in ft is {val_to}')

        page.screenshot(path='verification/converter.png')

        browser.close()

if __name__ == '__main__':
    run_verification()
