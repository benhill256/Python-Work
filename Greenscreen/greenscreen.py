from PIL import Image

will_repeat = True
while will_repeat:


    # These functions are to define which of the available pictures will be spliced.
    # The "while" loops here are to account for unrecognized inputs.
    print("""\nThe following are the options for the foreground:
    boat
    cactus
    cat
    cat_small
    harvester
    hiker
    penguin
    spaceshuttle""")
    green = input("\nWhich image in the folder would you like to superimpose? ")
    while green.lower() not in ["boat", "cactus", "cat", "cat_small", "harvester", "hiker", "penguin", "spaceshuttle"]:
        print("Image not found. Please try again.")
        green = input("Which image in the folder would you like to superimpose? ")

    print("""\nThe following are the options for the background:
    beach
    blossoms
    building
    clouds
    desert
    earth
    field
    forest
    mirror
    sails
    sky
    snowscape
    sunrise
    sunset
    tiles
    tree""")
    base = input("\nWhich image in the folder would you like to be the background? ")
    while base.lower() not in ["beach", "blossoms", "building", "clouds", "desert", "earth", "field", "forest", "mirror", "sails", "sky", "snowscape", "sunrise", "sunset", "tiles", "tree"]:
        print("Image not found. Please try again.")
        base = input("Which image in the folder would you like to be the background? ")

    image_green = Image.open(f"{green}.jpg")
    width_green, height_green = image_green.size
    pixels_green = image_green.load()

    image_base = Image.open(f"{base}.jpg")
    width_base, height_base = image_base.size
    pixels_base = image_base.load()

    # These are scaling constants, for the case that the size of
    # the background is different from the size of the overlay.
    k1 = int(width_base / width_green)
    k2 = int(height_base / height_green)

    # This does the actual splicing of the two images, with the
    # scaling constants acting on the background.
    for i in range(width_green):
        for j in range(height_green):
            r, g, b = pixels_green[i, j]
            if r < 160 and g > 160 and b < 160:
                pixels_green[i, j] = pixels_base[k1 * i, k2 * j]

    image_green.show()

    # This saves the image with a custom name, dependent on the
    # names of the two images being combined.
    save_image = input("\nWould you like to save the image (YES / NO)? ").lower()
    while save_image not in ["yes", "no"]:
        save_image = input("Input not recognized. Please try again: ").lower()
    if save_image == "yes":
        filename = input("\nWould you like to input a custom filename, or use an auto-generated one (CUSTOM / AUTO)? ").lower()
        while filename not in ["custom", "auto"]:
            filename = input("Input not recognized. Please try again: ").lower()
        if filename == "custom":
            filename = input("\nPlease enter the desired filename (not including the file type extension): ")
        elif filename == "auto":
            filename = f"{green}_{base}"
        image_green.save(f"{filename}.jpg")
        print(f"\nCombined image saved under the filename '{filename}.jpg'.")


    # This asks if the user would like to use the program again.
    repeat = input("\nWould you like to combine more images (YES / NO)? ").lower()
    while repeat not in ["yes", "no"]:
        repeat = input("Input not recognized. Please try again: ").lower()
    if repeat == "no":
        will_repeat = False

print()


# r, g, b = pixels_green[50, 100]
# r, g, b = 255, 0, 255
# pixels_green[50, 100] = (r, g, b)