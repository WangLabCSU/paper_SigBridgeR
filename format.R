flir::fix(path = rstudioapi::getActiveDocumentContext()$path)
system2(
  command = "air",
  args = c("format", rstudioapi::getActiveDocumentContext()$path)
)
