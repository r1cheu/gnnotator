def container_image(image_name):
    """
    Helper function to change mirror prefixes for container images.
    """
    return f"docker://{config["container_mirror"]}/{image_name}"
