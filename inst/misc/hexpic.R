library(hexSticker)
library(ggplot2)
library(ggimage)

p <- ggplot(data = tibble::tibble(x = 0, y = 0), aes(x, y)) +
  geom_bgimage("clarabel_logo.png") + theme_void()

## sticker(p, package = "clarabel", p_color = "#b7410e", p_family = "sans", p_fontface = "italic",
##         p_size = 16, p_y = 1.55, s_x = 1.0, s_y = .90, s_width = 1, s_height = 1,
##         h_fill = "white", h_size = 0.75, h_color = "#b7410e", filename = "clarabel.png")

sticker(p, package = "clarabel", p_color = "#b7410e", p_family = "sans", p_fontface = "bold",
        p_size = 16, p_y = 1.55, s_x = 1.0, s_y = .90, s_width = 1, s_height = 1,
        h_fill = "white", h_size = 1.5, h_color = "#b7410e", filename = "clarabel.png")






