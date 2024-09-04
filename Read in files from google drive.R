## ----------------------------------- ##
# Housekeeping ----
## ----------------------------------- ##

# Identify URLs of Drive folders with needed data
main_url <- googledrive::as_id("https://drive.google.com/drive/u/1/folders/0AIPkWhVuXjqFUk9PVA") #main WG folder, contains Site_Reference_Table
c_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/1dTENIB5W2ClgW0z-8NbjqARiaGO2_A7W") #Master_Chemistry, contains masterdata_chem
q_url <- googledrive::as_id("https://drive.google.com/drive/u/0/folders/1hbkUsTdo4WAEUnlPReOUuXdeeXm92mg-") #Master_Dataset folder, contains master Q

# Identify needed data in the Drive
wanted_files <- googledrive::drive_ls(path = main_url) %>%
  #dplyr::bind_rows(googledrive::drive_ls(path = q_url)) %>%
  dplyr::bind_rows(googledrive::drive_ls(path = c_url)) %>%
  # Filter to only needed files
  dplyr::filter(name %in% c("Site_Reference_Table", #saved as Google sheet, no file extension needed?
                            "20240201_masterdata_discharge.csv",
                            "20240312_masterdata_chem.csv"))

# Check those files
wanted_files
#Site Reference Table? Need to manually download

# Create a folder to download data into
dir.create(path = file.path("raw_data"), showWarnings = F)

# Download that data
purrr::walk2(.x = wanted_files$name, .y = wanted_files$id,
             .f = ~ googledrive::drive_download(file = .y, overwrite = T,
                                                path = file.path("raw_data", .x)))
