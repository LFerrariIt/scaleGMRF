# Install required packages (if not already installed)
#install.packages(c("hexSticker", "showtext"))

# Load libraries
library(hexSticker)
library(showtext)


# Path to your white scale icon (must be a white icon on transparent background)
# Replace with your actual file path or name
img_path <- "scale_icon_white.png"

# Create the hex sticker          # position of image
  sticker(
subplot = img_path,           # icon image
package = "scaleGMRF",        # package name
p_size = 130,  # Font size for the title
p_y = 0.6,  # Vertical position of the title
p_family = "mono",  # Typewriter-like font
p_color = "white",  # Title color
p_fontface = "bold",  # Make the text bold
h_color = "#16684e",
h_fill="#23A77E",# Set the border color to white
s_x = 1,  # Horizontal position of the image
s_y = 1.25,  # Vertical position of the image
s_width = 0.5,  # Size of the image
filename = "sticker.png",  # Output file path
url="github.com/LFerrariIt/scaleGMRF",
u_size=25,
u_color="#16684e",
dpi=2000
)

