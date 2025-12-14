
from playwright.sync_api import sync_playwright
import os

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        # Load local index.html using file:// protocol
        path = os.path.abspath('index.html')
        page.goto(f'file://{path}')

        # Use ID selector for robustness
        page.click('#mode-toggle')

        # Open Settings Tool
        page.wait_for_selector('#app-grid', state='visible')
        page.locator('.app-card').filter(has_text='Settings').click()

        # Toggle Dark Mode
        page.wait_for_selector('#tool-settings', state='visible')
        page.locator('label.switch').click()

        # Verify Dark Mode Class on Body
        body_class = page.evaluate('document.body.className')
        print(f'Body Class: {body_class}')
        assert 'dark-mode' in body_class

        # Take Screenshot of Dark Mode Settings
        page.screenshot(path='verification/dark_mode.png')

        # Back button is inside the specific tool view
        # We need to click the back button INSIDE #tool-settings
        page.click('#tool-settings .back-btn')

        # Wait for grid again
        page.wait_for_selector('#app-grid', state='visible')

        # Open Statistics Tool
        page.locator('.app-card').filter(has_text='Statistics').click()

        # Open Distributions Tab
        page.wait_for_selector('#tool-statistics', state='visible')
        page.click("//button[contains(text(), 'Distributions')]")

        # Select Chi-Squared
        page.select_option('#stat-dist-type', 'chisquare')

        # Verify Input visible
        expect_chi = page.locator('#dist-inputs-chisquare')
        if expect_chi.is_visible():
            print('Chi-Square inputs visible')
        else:
            print('Chi-Square inputs NOT visible')

        # Take Screenshot of Chi-Square Plot
        page.wait_for_timeout(500)
        page.screenshot(path='verification/chisquare.png')

        browser.close()

if __name__ == '__main__':
    run()
