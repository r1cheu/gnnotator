def container_image(image_name):
    """
    Helper function to change mirror prefixes for container images.
    """
    if config["container_mirror"] == "":
        return f"docker://{image_name}"
    return f"docker://{config["container_mirror"]}/{image_name}"
