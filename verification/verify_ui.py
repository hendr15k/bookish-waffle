
from playwright.sync_api import sync_playwright
import os

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        # Load local index.html
        cwd = os.getcwd()
        page.goto(f"file://{cwd}/index.html")

        # 1. Test Dark Mode
        # Click "App Mode" to see tools
        page.click("#mode-toggle")

        # Scroll to Settings
        page.locator(".app-card").filter(has_text="Settings").click()

        # Toggle Dark Mode
        page.click("label[for=dark-mode-toggle]")

        # Take Screenshot of Dark Mode
        page.screenshot(path="verification/dark_mode.png")

        # 2. Test Chi-Squared
        # Go back (close tool)
        page.click(".back-btn")

        # Open Statistics
        page.locator(".app-card").filter(has_text="Statistics").click()

        # Click Distributions Tab
        page.click("button:has-text(\"Distributions\")")

        # Select Chi-Squared
        page.select_option("#stat-dist-type", "chisquare")

        # Wait for canvas update (approx)
        page.wait_for_timeout(1000)

        # Take Screenshot of Chi-Squared
        page.screenshot(path="verification/chisquare.png")

        browser.close()

if __name__ == "__main__":
    run()
