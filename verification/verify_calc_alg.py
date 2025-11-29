
from playwright.sync_api import sync_playwright
import os

def run_verification():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        cwd = os.getcwd()
        page.goto(f'file://{cwd}/index.html')

        # Enable App Mode
        page.click('#mode-toggle')

        # 1. Verify Calculator Plot
        page.click('.app-card:has-text("Calculator")')
        page.wait_for_selector('#tool-calculator', state='visible')
        page.fill('#calc-main-input', 'x^2')
        page.click('#tool-calculator button:text("Plot")')

        page.wait_for_selector('#tool-calculator .tool-result canvas')
        print("Canvas found for Calculator Plot")
        page.screenshot(path='verification/calc_plot.png')

        # Back
        page.click('#tool-calculator .back-btn')

        # 2. Verify Algebra Plot
        page.click('.app-card:has-text("Algebra")')
        page.wait_for_selector('#tool-algebra', state='visible')
        page.fill('#alg-expr', 'sin(x)')
        page.click('#tool-algebra button:text("Plot")')

        page.wait_for_selector('#tool-algebra .tool-result canvas')
        print("Canvas found for Algebra Plot")
        page.screenshot(path='verification/algebra_plot.png')

        browser.close()

if __name__ == '__main__':
    run_verification()
