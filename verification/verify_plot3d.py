
from playwright.sync_api import sync_playwright
import os

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        # Load local index.html
        cwd = os.getcwd()
        page.goto(f"file://{cwd}/index.html")

        # We are in Desktop mode by default.
        # Enter command in main input.
        page.fill("#command-input", "plot3d(sin(sqrt(x^2+y^2)), x, y, -5, 5, -5, 5)")
        page.click("#submit-btn")

        # Wait for history entry
        page.wait_for_selector(".result canvas")

        # Wait for rendering
        page.wait_for_timeout(1000)

        page.screenshot(path="verification/plot3d_verification.png")
        browser.close()

if __name__ == "__main__":
    run()
