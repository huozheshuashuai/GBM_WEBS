p = r"E:/KIMI/GBM_WEBS_Project/113_glass_second_axis.R"
s = open(p, encoding="utf-8").read()
old = '''match_col <- function(bc) {
  hit <- grep(paste0("^", gsub("([\\-])", "\\\\\\1", bc)), colnames(m), value = TRUE)
  if (!length(hit)) hit <- grep(substr(bc, 1, 19), colnames(m), value = TRUE)
  if (length(hit) >= 1) hit[1] else NA
}'''
new = '''match_col <- function(bc) {
  bc2 <- gsub("-", ".", bc, fixed = TRUE)
  hit <- grep(paste0("^", bc2), colnames(m), value = TRUE)
  if (!length(hit)) {
    bc3 <- gsub("-", ".", substr(bc, 1, 19), fixed = TRUE)
    hit <- grep(paste0("^", bc3), colnames(m), value = TRUE)
  }
  if (length(hit) >= 1) hit[1] else NA
}'''
assert old in s, "pattern not found"
s = s.replace(old, new)
open(p, "w", encoding="utf-8").write(s)
print("patched")
