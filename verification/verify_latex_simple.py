
from playwright.sync_api import sync_playwright
import os

def run():
    cwd = os.getcwd()
    file_url = f"file://{cwd}/index.html"

    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page()
        page.goto(file_url)

        # Test asin(x) -> \arcsin(x)
        page.fill("#command-input", "asin(x)")
        page.click("#submit-btn")
        page.wait_for_timeout(1000)

        content = page.content()
        if "\\arcsin" in content:
             print("Verified: arcsin present in DOM")
        else:
             print("FAILED: arcsin NOT found in DOM")

        # Test abs(x) -> |x|
        page.fill("#command-input", "abs(x)")
        page.click("#submit-btn")
        page.wait_for_timeout(1000)
        content = page.content()
        if "\\left|x\\right|" in content:
             print("Verified: abs present in DOM")
        else:
             print("FAILED: abs NOT found in DOM")

        browser.close()

if __name__ == "__main__":
    run()
