import customtkinter as ctk
from PIL import Image, ImageTk

# Lazy Image Loader Class
class LazyImageLoader:
    """
    A class for lazily loading images. Images are only loaded when they are accessed,
    not when the instance of LazyImageLoader is created.
    """

    def __init__(self, image_path):
        """
        Initializes the LazyImageLoader with a specific image path.

        Args:
            image_path (str): The path to the image file.
        """
        self.image_path = image_path
        self._image = None  # Placeholder for the image, will be loaded on demand

    @property
    def image(self):
        """
        A property that loads and returns the image only when accessed. If the image is already loaded,
        it returns the existing image.

        Returns:
            CTkImage: The loaded image.
        """
        # Load the image only if it hasn't been loaded already
        if self._image is None:
            self._image = ctk.CTkImage(
                light_image=Image.open(self.image_path),
                dark_image=Image.open(self.image_path)
            )
        # Special case for a specific image path, adjusts size
        if self.image_path == ".\\assets\\logo_IPHC.png":
            return ctk.CTkImage(
                light_image=Image.open(self.image_path),
                size=(116, 80),
                dark_image=Image.open(self.image_path))
        elif self.image_path == ".\\assets\\juno.png":
            return ctk.CTkImage(
                light_image=Image.open(self.image_path),
                size=(96, 80),
                dark_image=Image.open(self.image_path))
        else: 
            # Return the loaded image
            return self._image