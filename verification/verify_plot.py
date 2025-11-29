
from playwright.sync_api import sync_playwright
import os

def run_verification():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        cwd = os.getcwd()
        page.goto(f'file://{cwd}/index.html')

        # Go to Calculus
        page.click('#mode-toggle')
        page.click('.app-card:has-text("Calculus")')
        page.wait_for_selector('#tool-calculus', state='visible')

        # Test Taylor Plot
        page.fill('#taylor-expr', 'sin(x)')
        page.click('button:text("Compare & Plot")')

        # Check result container has canvas
        page.wait_for_selector('.tool-result canvas')
        print("Canvas found for Taylor Plot")

        # Screenshot
        page.screenshot(path='verification/taylor_plot.png')

        browser.close()

if __name__ == '__main__':
    run_verification()
