
from playwright.sync_api import sync_playwright
import os

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        # Open local index.html
        url = 'file://' + os.path.abspath('index.html')
        page.goto(url)

        # 1. Switch to App Mode
        page.click('#mode-toggle')

        # 2. Open Statistics Tool
        page.click('.app-card:has-text("Statistics")')

        # 3. Open Distributions Tab
        page.click('button:has-text("Distributions")')

        # 4. Select Binomial
        page.select_option('#stat-dist-type', 'binomial')

        # 5. Take screenshot of Binomial Dist
        page.screenshot(path='verification/statistics_binomial.png')

        # 6. Go back - Be specific which back button (Statistics Tool Back Button)
        # The back button is inside .tool-header in #tool-statistics
        page.click('#tool-statistics .back-btn')

        # Wait for grid to be visible
        page.wait_for_selector('#app-grid', state='visible')

        # 7. Open Physics Tool
        page.click('.app-card:has-text("Physics")')

        # 8. Click Electricity Tab
        page.click('button:has-text("Electricity")')

        # 9. Calculate Ohm's Law (V=IR)
        # Ensure we are targeting the right select in physics tool
        page.select_option('#phys-elec-type', 'V')
        page.fill('#phys-I', '2')
        page.fill('#phys-R', '5')
        # Click Calculate inside the active tab/tool
        page.click('#tab-phys-elec button:has-text("Calculate")')

        # 10. Take screenshot of Physics Elec
        page.screenshot(path='verification/physics_electricity.png')

        browser.close()

if __name__ == '__main__':
    run()
