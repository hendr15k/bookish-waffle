
from playwright.sync_api import sync_playwright
import os

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        # Load local index.html using file:// protocol
        path = os.path.abspath('index.html')
        page.goto(f'file://{path}')

        # Click 'App Mode' button first to show the app grid
        page.click('#mode-toggle')

        # Click 'Physics' app
        page.click('.app-card:has-text("Physics")')

        # Click 'Wechselstrom (AC)' tab - Force visible check?
        # The tab content should be visible.
        page.click('button:has-text("Wechselstrom (AC)")')

        # Inputs inside tabs might be hidden if tab logic didn't fire correctly
        # Wait for tab content to be visible
        page.wait_for_selector('#tab-phys-ac', state='visible')

        # Select Calculation 'Reactance' (default)
        # Fill inputs
        page.fill('#phys-f', '50')
        page.fill('#phys-L', '0.1')

        # Click Calculate
        page.click('#tool-physics .action-btn:has-text("Calculate")')

        # Wait for result text to be populated (not empty)
        # Find the visible result box inside the active tool
        page.wait_for_function('document.querySelector("#tool-physics .tool-result").innerText.length > 0')

        # Take screenshot
        page.screenshot(path='verification/ac_reactance.png')

        browser.close()

if __name__ == '__main__':
    run()
