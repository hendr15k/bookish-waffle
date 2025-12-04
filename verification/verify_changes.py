from playwright.sync_api import sync_playwright
import os

def run(playwright):
    browser = playwright.chromium.launch(headless=True)
    page = browser.new_page()

    # Load index.html locally
    filepath = os.path.abspath("index.html")
    page.goto(f"file://{filepath}")

    # Check Physics Tool
    print("Checking Physics Tool...")
    page.click("text=App Mode")
    page.click("#app-grid > div:nth-child(13)") # Physics card (needs better selector maybe, but order is fixed)

    # Check Electricity Tab
    print("Checking Electricity Tab...")
    page.click("button:has-text('Electricity')")
    page.wait_for_selector("#tab-phys-elec.active")

    # Calculate Ohm's Law
    page.select_option("#phys-elec-type", "ohm")
    page.fill("#phys-I", "2")
    page.fill("#phys-R", "10")
    page.click("#tab-phys-elec button:has-text('Calculate')")

    # Verify result
    result = page.text_content("#tool-physics .tool-result")
    print(f"Ohm's Result: {result}")
    assert "20.0000 V" in result

    # Check Waves Tab
    print("Checking Waves Tab...")
    page.click("button:has-text('Waves')")
    page.wait_for_selector("#tab-phys-wave.active")

    # Calculate Wave Speed
    page.select_option("#phys-wave-type", "speed")
    page.fill("#phys-f", "50")
    page.fill("#phys-l", "2")
    page.click("#tab-phys-wave button:has-text('Calculate')")

    result = page.text_content("#tool-physics .tool-result")
    print(f"Wave Result: {result}")
    assert "100.0000 m/s" in result

    # Take screenshot of Physics
    page.screenshot(path="verification/physics_check.png")

    # Go back and check Geometry
    # Explicitly target the back button of the physics tool
    page.click("#tool-physics .back-btn")

    print("Checking Geometry Tool...")
    # Wait for animation/display change
    page.wait_for_selector("#app-grid", state="visible")

    page.click("#app-grid > div:nth-child(11)") # Geometry card
    page.click("button:has-text('3D Objects')")

    # Check Cone
    print("Checking Cone...")
    page.select_option("#geo-3d-shape", "cone")
    page.fill("#geo-r3", "3")
    page.fill("#geo-h3", "4")
    page.click("#tab-geo-3d button:has-text('Calculate')")

    result = page.text_content("#tool-geometry .tool-result")
    print(f"Cone Result: {result}")
    # Vol = 1/3 * pi * 9 * 4 = 12 pi approx 37.6991
    assert "Volume = 37.6991" in result

    # Check Pyramid
    print("Checking Pyramid...")
    page.select_option("#geo-3d-shape", "pyramid")
    page.fill("#geo-a3", "10")
    page.fill("#geo-h3", "12")
    page.click("#tab-geo-3d button:has-text('Calculate')")

    result = page.text_content("#tool-geometry .tool-result")
    print(f"Pyramid Result: {result}")
    # Vol = 100 * 12 / 3 = 400
    assert "Volume = 400.0000" in result

    # Take screenshot of Geometry
    page.screenshot(path="verification/geometry_check.png")

    browser.close()

with sync_playwright() as playwright:
    run(playwright)
