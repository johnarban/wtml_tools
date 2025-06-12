# wtml_tools
Tools to help deal with WTML, AVM, and WCS (stuff probably already done in other WWT packages)

## Basic Functionality

This codebase provides tools to generate `.wtml` files for use with the WorldWide Telescope (WWT). The primary functionality is to take an image file and its corresponding WCS (World Coordinate System) file and create a `.wtml` file. Below is an example of how to use the `create_wtml_from_image` function to achieve this:

```python
from wtml_tools import wtml_tools

# Input files
image_path = "example_image.jpg"
wcsfile = "example_image.wcs"

# Modifications to apply
mods = {
    'scale': {'x': 1.5, 'y': 1.5},  # Scale the image by a factor of 1.5 in both directions
}

# Generate WTML
wtml_tools.create_wtml_from_image(
    image_path=image_path,
    wcsfile=wcsfile,
    name="Example Image",
    image_url="https://example.com/example_image.jpg",
    thumb_url="https://example.com/example_thumbnail.jpg",
    description="This is an example image with scale modifications.",
    credits="NASA, ESA",
    credits_url="https://example.com/credits",
    mods=mods,
)
```

## Advanced Usage

### Modifying WTML Files

You can add or update metadata in an existing `.wtml` file using the `add_description_credits_to_wtml` function:

```python
from wtml_tools import wtml_tools

wtml_file = "example.wtml"
wtml_tools.add_description_credits_to_wtml(
    wtml_file=wtml_file,
    description="Updated description",
    credits="Updated credits",
    credits_url="https://updated-credits.com",
    place_name="Updated Place Name",
    imageset_name="Updated ImageSet Name",
)
print(f"Updated WTML file: {wtml_file}")
```

### Applying Modifications to WCS Headers

The `create_wtml_from_image` function allows you to modify the WCS header before generating the `.wtml` file. For example, you can scale the image using the `mods` parameter:

```python
from wtml_tools import wtml_tools

# Input files
image_path = "example_image.jpg"
wcsfile = "example_image.wcs"

# Modifications to apply
mods = {
    'scale': {'x': 2.0, 'y': 2.0},  # Scale the image by a factor of 2 in both directions
}

# Generate WTML with modifications
tree, imageset, out, image_header = wtml_tools.create_wtml_from_image(
    image_path=image_path,
    wcsfile=wcsfile,
    name="Scaled Image Example",
    image_url="https://example.com/scaled_image.jpg",
    thumb_url="https://example.com/scaled_thumbnail.jpg",
    description="This is a scaled image example.",
    credits="NASA, ESA",
    credits_url="https://example.com/credits",
    mods=mods,
)
print(f"WTML file created with modifications: {out}")
```
