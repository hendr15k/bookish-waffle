
from playwright.sync_api import sync_playwright
import os

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        cwd = os.getcwd()
        page.goto(f"file://{cwd}/index.html")

        # 1. Click App Mode
        page.click("#mode-toggle")

        # 2. Click Statistics Tool
        page.click(".app-card:has-text('Statistics')")

        # 3. Click Regression Tab
        page.click("button:has-text('Regression')")

        # 4. Select Quadratic
        page.select_option("#reg-type", "quadratic")

        # 5. Enter Data
        page.fill("#reg-data", "[[-1,1],[0,0],[1,1],[2,4]]")

        # 6. Click Calculate (Target the one inside regression tab)
        page.click("#tab-stat-reg button:has-text('Calculate & Plot')")

        # Wait for result
        page.wait_for_selector("#tool-statistics .tool-result canvas")
        page.wait_for_timeout(500)

        page.screenshot(path="verification/stats_regression.png")

        # 7. Click Back button inside Statistics Tool
        page.click("#tool-statistics .back-btn")

        # 8. Click Physics Tool
        page.click(".app-card:has-text('Physics')")

        # 9. Click Constants Tab
        page.click("button:has-text('Constants')")

        page.wait_for_timeout(200)
        page.screenshot(path="verification/physics_constants.png")

        browser.close()

if __name__ == "__main__":
    run()
