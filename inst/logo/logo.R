img <- image_read("pills.png")

fonts <- font_add_google("Zilla Slab")

font <- "Zilla Slab"

x <- sticker(
  subplot = img,
  package = "BayesDissolution",
  s_width = 2.5,
  s_height = 3.5,
  s_x = 1,
  s_y = 0.85,
  p_size = 5,
  h_fill = "#3F8CCC",
  h_color = "darkgray",
  h_size = 1.5,
  #spotlight = T
  #l_alpha = 1,   ### l_ parameters are for the spotlight
  p_color = "white",
  #l_y = 1,
  #l_x = 1,
  #l_width = 3,
  #l_height = 3,
  p_family = font
)
