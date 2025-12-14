
from playwright.sync_api import sync_playwright
import os

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        # Load local index.html
        cwd = os.getcwd()
        page.goto(f"file://{cwd}/index.html")

        # 1. Test Light Mode (Default)
        page.screenshot(path="verification/light_mode.png")

        # 2. Test Dark Mode
        page.click("#mode-toggle")
        page.locator(".app-card").filter(has_text="Settings").click()
        page.click("label[for=dark-mode-toggle]")
        page.wait_for_timeout(500)
        page.screenshot(path="verification/dark_mode_fixed.png")

        browser.close()

if __name__ == "__main__":
    run()
