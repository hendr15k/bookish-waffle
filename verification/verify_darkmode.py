from playwright.sync_api import sync_playwright
import os

def run(playwright):
    browser = playwright.chromium.launch()
    page = browser.new_page()

    cwd = os.getcwd()
    page.goto(f"file://{cwd}/index.html")

    # Enable Mobile Mode
    page.click("#mode-toggle")

    # Open Settings
    page.click("text=Settings")

    # Toggle Dark Mode
    # Use selector for checkbox id="dark-mode-toggle"
    # Click the label or input
    page.check("#dark-mode-toggle")

    # Wait for body class
    page.wait_for_selector("body.dark-mode")

    # Screenshot Settings Page in Dark Mode
    page.screenshot(path="verification/darkmode_settings.png")

    # Go back and screenshot the main app grid in Dark Mode
    page.click(".back-btn")

    # Wait for animation/transition?
    page.wait_for_timeout(500)
    page.screenshot(path="verification/darkmode_home.png")

    browser.close()

with sync_playwright() as playwright:
    run(playwright)
