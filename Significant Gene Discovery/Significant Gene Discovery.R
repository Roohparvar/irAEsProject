install.packages("writexl")
install.packages("dplyr")

library(dplyr)
library(writexl)

# خواندن فایل CSV
data <- read.csv("AllMarkers.csv")  # اینجا اسم فایل CSV خودت رو بذار

# بررسی مقادیر یکتای ستون cluster
clusters <- unique(data$cluster)

# ساختن یک لیست برای نگهداری داده‌های هر کلاستر
sheets <- list()

# پر کردن لیست با داده‌های جداشده بر اساس cluster
for (cl in clusters) {
  sheets[[paste0("Cluster_", cl)]] <- filter(data, cluster == cl)
}

# نوشتن خروجی در یک فایل اکسل
write_xlsx(sheets, path = "clustered_output.xlsx")

# ....................................................................................

library(openxlsx)

# فایل اکسل
input_file <- "clustered_output.xlsx"
output_file <- "filtered_output.xlsx"

# لیست شیت‌های موجود
sheet_names <- getSheetNames(input_file)

# باز کردن یک فایل جدید اکسل برای ذخیره نتایج فیلتر شده
wb <- createWorkbook()

# برای هر شیت، فیلتر کردن داده‌ها و ذخیره در شیت جدید
for (sheet_name in sheet_names) {
  # خواندن داده‌ها از شیت فعلی
  data <- read.xlsx(input_file, sheet = sheet_name)
  
  # فیلتر کردن داده‌ها بر اساس شرایط مشخص شده
  filtered_data <- subset(data, 
                          p_val_adj < 0.05 & 
                            avg_log2FC > 0.1 & 
                            pct.1 > 0.1 & 
                            (pct.1 - pct.2) > 0.1)
  
  # اضافه کردن شیت جدید با داده‌های فیلتر شده
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = filtered_data)
}

# ذخیره کردن فایل اکسل جدید
saveWorkbook(wb, output_file, overwrite = TRUE)

# ....................................................................................