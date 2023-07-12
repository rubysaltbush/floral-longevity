# function to cache pre-prepared R data. If csv already in cache will read data
# from cache, if not will run defined script to produce and then cache data

cache_csv <- function(path, fetch_function, read_function = readr::read_csv, 
                      save_function = readr::write_csv) {
  if (file.exists(path)) {
    read_function(path, guess_max = 10000)
  } else {
    data <- fetch_function()
    save_function(data, path)
    data
  }
}
