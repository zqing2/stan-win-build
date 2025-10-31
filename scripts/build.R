# scripts/build.R
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("cmdstanr", quietly = TRUE)) install.packages("cmdstanr")
library(cmdstanr)

# 安装 CmdStan（Windows Runner 上很快）
install_cmdstan(cores = 2, quiet = TRUE)

# 编译（开启线程）
stan_file <- file.path("model", "hier_nb_mm_group_parallel_full.stan")
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

# 目标可执行文件路径
exe_path <- mod$exe_file()  # e.g., model/hier_nb_mm_group_parallel_full.exe

# ====== 收集运行时 DLL（TBB 与 MinGW）======
# TBB DLL（oneTBB）
tbb_dir <- file.path(cmdstan_path(), "stan", "lib", "stan_math", "lib", "tbb")
tbb_dlls <- if (dir.exists(tbb_dir)) list.files(tbb_dir, "\\.dll$", TRUE, TRUE) else character(0)

# MinGW 运行时 DLL（libstdc++/libgcc/libwinpthread）
rtools_home <- Sys.getenv("RTOOLS43_HOME", Sys.getenv("RTOOLS40_HOME", ""))
cand <- c(file.path(rtools_home, "ucrt64", "bin"), file.path(rtools_home, "mingw64", "bin"))
mingw_bin <- cand[dir.exists(cand)][1]
mingw_dlls <- character(0)
if (!is.na(mingw_bin)) {
  pats <- c("libstdc\\+\\+.*\\.dll$", "libgcc_s_.*\\.dll$", "libwinpthread-.*\\.dll$")
  mingw_dlls <- unlist(lapply(pats, function(p) list.files(mingw_bin, p, TRUE)))
}

# 打包输出
out_dir <- "dist"; if (!dir.exists(out_dir)) dir.create(out_dir)
file.copy(exe_path, file.path(out_dir, basename(exe_path)), overwrite = TRUE)
file.copy(tbb_dlls, out_dir, overwrite = TRUE)
file.copy(mingw_dlls, out_dir, overwrite = TRUE)

# 写运行说明
readme <- c(
  "# Windows 运行说明",
  "",
  "把整个压缩包解压到一个文件夹，确保以下文件与 .exe 位于同一目录：",
  "- tbb*.dll",
  "- libstdc++-6.dll / libgcc_s_seh-1.dll / libwinpthread-1.dll",
  "",
  "R 里使用：",
  "```r",
  "library(cmdstanr)",
  sprintf("mod <- cmdstan_model(exe_file = \"%s\")", basename(exe_path)),
  "Sys.setenv(STAN_NUM_THREADS=8)  # 如需并行",
  "fit <- mod$variational(data = stan_data, seed = 123)",
  "```",
  "",
  "命令行使用：",
  sprintf("%s method=variational data file=your_data.json", basename(exe_path)),
  "",
  "自检：",
  sprintf(".\\%s --help", basename(exe_path))
)
writeLines(readme, file.path(out_dir, "README_windows.md"))

# 压缩为 artifact
zipfile <- "stan_windows_bundle.zip"
if (file.exists(zipfile)) file.remove(zipfile)
old <- getwd(); setwd(out_dir)
zip(zipfile = file.path("..", zipfile), files = list.files("."))
setwd(old)

cat("\n== Done ==\nBundle:", zipfile, "\n")
